/*******************************************************************************************
 *
 *  filters overlaps by various criteria
 *
 *  Author :  MARVEL Team
 *
 *  Date   :  May 2015
 *
 *******************************************************************************************/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/pass.h"
#include "lib/read_loader.h"
#include "lib/tracks.h"
#include "lib/trim.h"
#include "lib/utils.h"

#include "dalign/align.h"
#include "db/DB.h"

#define DEF_ARG_R TRACK_REPEATS
#define DEF_ARG_T TRACK_TRIM
#define DEF_ARG_L TRACK_DUST
#define DEF_ARG_F 0
#define DEF_ARG_C 0
#define DEF_ARG_U 2000
#define DEF_ARG_M 2400
#define DEF_ARG_W 600
#define DEF_ARG_D 28
#define DEF_ARG_O 4000
#define DEF_ARG_Z 20000

#define VERBOSE

#undef DEBUG_CHAIN
#undef CHAIN_DEBUG
#undef DEBUG_GAPS
#undef DEBUG_FILTER

#define ANCHOR_INVALID 	(1 << 0)
#define ANCHOR_TRIM 		(1 << 1)
#define ANCHOR_LOWCOMP 	(1 << 2)

typedef struct
{
	Overlap **ovls;
	int novl;
	int maxOvl;

} Chain;

typedef struct
{

	int beg;
	int end;
	int flag;
} anchorItv;

typedef struct
{
	// stats counters
	int statsFiltRepeat;
	int statsFiltContained;
	int statsFiltInvalidChain;
	int statsLowCovALn;
	int statsGapAlns;

	// settings
	int nMinNonRepeatBases;
	int nFuzzBases;
	int nContPerc;
	int nVerbose;
	int keepIdentity;

	int maxOverallDiff;		// add up all diffs from all chains divided by number of overlapping bases in A- and B-read

	// repeat settings
	int mergeRepDist;
	int repeatWindowLookBack;
	int rp_mergeTips;
	int minChainLen;

	// stitch settings
	int stitchChain;
	int stitchMaxTipFuzzy;
	int stitchNonLowCompAnchorBases;
	int stitchMaxGapSize;
	int stitchMaxGapSizeInLowCompl;
	int stitchMinChainLen;
	int stitchMaxChainLASs;

	int gapMinSpanners;
	int minTipCoverage;

	HITS_DB* db;
	HITS_TRACK* trackRepeat;
	HITS_TRACK* trackTrim;
	HITS_TRACK* trackLowCompl;

	FILE* fileOutDiscardedReads;
	FILE* fileOutFullyDiscardedAreads;
	int   minLenOfFullyDiscardedAreads;

	ovl_header_twidth twidth;

	int nkeptChains; // 0 ... only best, otherwise keep all
	Chain *ovlChains;
	int curChains;
	int maxChains;

	int maxUniqAIntervals;
	int curUniqAIntervals;
	anchorItv *uniqAIntervals;

	int maxUniqBIntervals;
	int curUniqBIntervals;
	anchorItv *uniqBIntervals;

	TRIM* trim;
} FilterContext;

extern char* optarg;
extern int optind, opterr, optopt;

static int getRepeatBases(FilterContext *ctx, Overlap *ovl, int read)
{
	if (ctx->trackRepeat == NULL)
	{
		return 0;
	}

	assert(ovl->aread == read || ovl->bread == read);

//	int aLen = ovl->path.aepos - ovl->path.abpos;
	int bLen = ovl->path.bepos - ovl->path.bbpos;

	// get repeats track
	track_anno* rep_anno = ctx->trackRepeat->anno;
	track_data* rep_data = ctx->trackRepeat->data;

	track_anno rb, re;
	int nrep = 0;

	rb = rep_anno[read] / sizeof(track_data);
	re = rep_anno[read + 1] / sizeof(track_data);

	// loop through all repeats in a
	int rBeg, rEnd;
	while (rb < re)
	{
		rBeg = rep_data[rb];
		rEnd = rep_data[rb + 1];

		if (ovl->aread == read)
		{
			nrep += intersect(ovl->path.abpos, ovl->path.aepos, rBeg, rEnd);
		}
		else
		{
			if (ovl->flags & OVL_COMP)
			{
				nrep += intersect(bLen - ovl->path.bepos, bLen - ovl->path.bbpos, rBeg, rEnd);
			}
			else
			{
				nrep += intersect(ovl->path.bbpos, ovl->path.bepos, rBeg, rEnd);
			}

		}
		rb += 2;
	}
	return nrep;
}

static int cmp_chain_len(const void *a, const void *b)
{
	Chain* c1 = (Chain *) a;
	Chain* c2 = (Chain *) b;

	int i;

	int olen1 = c1->ovls[0]->path.aepos - c1->ovls[0]->path.abpos;

	for (i = 1; i < c1->novl; i++)
	{
		olen1 += c1->ovls[i]->path.aepos - c1->ovls[i]->path.abpos;

		if (c1->ovls[i - 1]->path.aepos > c1->ovls[i]->path.abpos)
			olen1 -= c1->ovls[i - 1]->path.aepos > c1->ovls[i]->path.abpos;
	}

	int olen2 = c2->ovls[0]->path.aepos - c2->ovls[0]->path.abpos;

	for (i = 1; i < c2->novl; i++)
	{
		olen2 += c2->ovls[i]->path.aepos - c2->ovls[i]->path.abpos;

		if (c2->ovls[i - 1]->path.aepos > c2->ovls[i]->path.abpos)
			olen2 -= c2->ovls[i - 1]->path.aepos > c2->ovls[i]->path.abpos;
	}

	return (olen2 - olen1);
}

static int cmp_ovls_abeg(const void* a, const void* b)
{
	Overlap* o1 = *(Overlap**) a;
	Overlap* o2 = *(Overlap**) b;

	int cmp = o1->path.abpos - o2->path.abpos;

	if (!cmp)
	{
		cmp = (o1->path.aepos - o1->path.abpos) - (o2->path.aepos - o2->path.abpos);
	}

	return cmp;
}

#ifdef DEBUG_CHAIN
static void printChain(Chain *chain)
{
	printf("CHAIN: nvols %d, %7d vs %7d %s\n", chain->novl, chain->ovls[0]->aread, chain->ovls[0]->bread, (chain->ovls[0]->flags & OVL_COMP) ? "COMP" : "NORM");
	int i;
	for (i = 0; i < chain->novl; ++i)
	{
		printf("%3d in A [%8d,%8d] in B [%8d,%8d] %5.2f\n", i, chain->ovls[i]->path.abpos, chain->ovls[i]->path.aepos, chain->ovls[i]->path.bbpos,
				chain->ovls[i]->path.bepos, chain->ovls[i]->path.diffs * 100.0 / (chain->ovls[i]->path.aepos - chain->ovls[i]->path.abpos));
	}
}
#endif

static int contained(int ab, int ae, int bb, int be)
{
	if (ab >= bb && ae <= be)
	{
		return 1;
	}

	return 0;
}

static void chain(FilterContext *ctx, Overlap *ovls, int n)
{
	/// TODO hard coded
	int MIN_OVL_LOOKAHEAD = 2000;
	int MAX_OVL_LOOKAHEAD = 30000;
	int STRIDE_OVL_LOOKAHEAD = 2000;

#ifdef DEBUG_CHAIN
	printf("chain(%d,%d,%d) CHAIN: n%d m%d\n", ovls->aread, ovls->bread, n, ctx->curChains, ctx->maxChains);
#endif
	if (n < 2)
	{

		if (ctx->ovlChains->ovls == NULL)
		{
			ctx->ovlChains->novl = 0;
			ctx->ovlChains->maxOvl = 10;
			ctx->ovlChains->ovls = (Overlap**) malloc(sizeof(Overlap*) * ctx->ovlChains->maxOvl);
		}

		ctx->ovlChains->ovls[0] = ovls;
		ovls->flags |= OVL_TEMP;
		ctx->ovlChains->novl++;
		ctx->curChains++;

#ifdef DEBUG_CHAIN
		printChain(ctx->ovlChains);
#endif

		return;
	}

	{
		int i, j;
		// get rid of contained overlaps
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = ovls + i;

			if (ovl_i->flags & (OVL_CONT))
				continue;

			for (j = i + 1; j < n; j++)
			{
				Overlap *ovl_j = ovls + j;

				if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos) && contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
				{
					ovl_j->flags |= (OVL_CONT);
					ctx->statsFiltContained++;
				}
			}
		}
	}

	int aread, bread;
	int alen, blen;
	int i;

	aread = ovls->aread;
	bread = ovls->bread;

	alen = DB_READ_LEN(ctx->db, aread);
	blen = DB_READ_LEN(ctx->db, bread);

	int nremain = n;

#ifdef DEBUG_CHAIN
	printf("find detect already excluded overlaps\n");
#endif
	{
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = ovls + i;

			if (ovl_i->flags & (OVL_CONT | OVL_TRIM))
			{
				nremain--;
			}
		}
	}
#ifdef DEBUG_CHAIN
	printf("nremain %d\n", nremain);
#endif

	if (nremain == 0)
		return;

	int round = 0;

	while (nremain > 0) // use nremain to find putative chain starting ovls
	{
		round++;

		int longestUniqOvlBases = -1;
		int longestUniqOvlIdx = -1;
		int longestOvlBases = -1;
		int longestOvlIdx = -1;

		// find longest overlap based on number of unique bases
		for (i = 0; i < n; i++)
		{
			Overlap *ovl = ovls + i;

			if (ovl->flags & (OVL_CONT | OVL_TEMP | OVL_REPEAT | OVL_TRIM))
			{
				continue;
			}

			int aLen = ovl->path.aepos - ovl->path.abpos;
			int bLen = ovl->path.bepos - ovl->path.bbpos;

			int aRep = getRepeatBases(ctx, ovl, ovl->aread);
			int bRep = getRepeatBases(ctx, ovl, ovl->bread);

#ifdef DEBUG_CHAIN
			printf("%d - %d [%d, %d] [%d, %d], aR %d/%d, bR %d/%d\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, aLen, aRep,
					bLen, bRep);
#endif
			int tmpBases = MAX(aLen - aRep, bLen - bRep);
			if (tmpBases > longestUniqOvlBases)
			{
				longestUniqOvlBases = tmpBases;
				longestUniqOvlIdx = i;
			}

			tmpBases = MAX(aLen, bLen);
			if (tmpBases > longestOvlBases)
			{
				longestOvlBases = tmpBases;
				longestOvlIdx = i;
			}
		}

		if (longestOvlIdx < 0)
		{
			//printf("Unable to find new chain start ovl! nremain %d, curChains %d\n", nremain, ctx->curChains);
			return;
		}

		if (longestUniqOvlBases < ctx->nMinNonRepeatBases + 1)
		{
#ifdef DEBUG_CHAIN
			printf("Number of unique bases to low. Use longest overlap.\n");
#endif
			if (round > 1)
			{
#ifdef DEBUG_CHAIN
				printf("Break out of chain. Cannot find unique anchor alignment to start a chain.\n");
#endif
				return;
			}

			longestUniqOvlBases = longestOvlBases;
			longestUniqOvlIdx = longestOvlIdx;
		}

#ifdef DEBUG_CHAIN
		printf("longest overlap:\n");
		printf("idx: %d --> uB %d, %d - %d [%d, %d] [%d, %d]\n", longestUniqOvlIdx, longestUniqOvlBases, ovls[longestUniqOvlIdx].aread,
				ovls[longestUniqOvlIdx].bread, ovls[longestUniqOvlIdx].path.abpos, ovls[longestUniqOvlIdx].path.aepos, ovls[longestUniqOvlIdx].path.bbpos,
				ovls[longestUniqOvlIdx].path.bepos);
#endif
// try to "elongate" longest overlap
// 1st on the right
// 2nd on the left side

		if (ctx->curChains == ctx->maxChains)
		{
			ctx->maxChains = ctx->maxChains * 1.2 + 5;
			ctx->ovlChains = (Chain*) realloc(ctx->ovlChains, sizeof(Chain) * ctx->maxChains);
			bzero(ctx->ovlChains + ctx->curChains, sizeof(Chain) * (ctx->maxChains - ctx->curChains));
		}

		Chain *chain = ctx->ovlChains + ctx->curChains;
#ifdef DEBUG_CHAIN
		printf("chain: nOvl: %d, maxOvl %d, nremain: %d\n", chain->novl, chain->maxOvl, nremain);
#endif

		if (chain->ovls == NULL)
		{
			chain->novl = 0;
			chain->maxOvl = 10;
			chain->ovls = (Overlap**) malloc(sizeof(Overlap*) * chain->maxOvl);

#ifdef DEBUG_CHAIN
			printf("chain: nOvl: %d, maxOvl %d\n", chain->novl, chain->maxOvl);
#endif
		}

		chain->ovls[0] = ovls + longestUniqOvlIdx;
		chain->ovls[0]->flags |= OVL_TEMP;
		chain->novl++;
		nremain--;
#ifdef DEBUG_CHAIN
		printf("chain: nOvl: %d, maxOvl %d, nremain %d\n", chain->novl, chain->maxOvl, nremain);
#endif

		int ab1, ae1;
		int bb1, be1;

		int ab2, ae2;
		int bb2, be2;

		if (nremain && longestUniqOvlIdx + 1 < n)
		{
			ab1 = ovls[longestUniqOvlIdx].path.abpos;
			ae1 = ovls[longestUniqOvlIdx].path.aepos;

			bb1 = ovls[longestUniqOvlIdx].path.bbpos;
			be1 = ovls[longestUniqOvlIdx].path.bepos;

#ifdef DEBUG_CHAIN
			printf("extend longest overlap in right direction\n");
#endif
// 1st right
			int cont = 1;
			int curBestUniqOffset = 1;
			int curBestUniqBases = -1;
			int curBestBases = -1;
			int curBestOffset = 1;
			int curBestIntersection = MAX(alen, blen);

			while (cont)
			{
				int stepSize;
				for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestUniqBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
				{
#ifdef DEBUG_CHAIN
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif
					for (i = longestUniqOvlIdx + curBestUniqOffset; i < n; i++)
					{
						Overlap * ovl = ovls + i;

						ab2 = ovl->path.abpos;
						ae2 = ovl->path.aepos;

						bb2 = ovl->path.bbpos;
						be2 = ovl->path.bepos;

						if ((ovl->flags & OVL_COMP) != (ovls[longestUniqOvlIdx].flags & OVL_COMP))
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> different orientations\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ovl->flags & OVL_CONT)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> really contained repeat\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ovl->flags & OVL_TEMP)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> is part of another chain\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ovl->flags & OVL_DISCARD)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> discarded\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						// todo mark those as discard ????
						if (contained(ab2, ae2, ab1, ae1))
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in A-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						// todo mark those as discard ????
						if (contained(bb2, be2, bb1, be1))
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in B-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ae2 < ae1 || be2 < be1) // also overlap must follow right direction
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> improper right extension direction (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2,
									bb2, be2, stepSize);
#endif
							continue;
						}

						if (MAX(ab2 - ae1, bb2 - be1) > stepSize)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> gap size too large (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
									stepSize);
#endif
							continue;
						}

						if (ae1 - ab2 > ae2 - ae1 || be1 - bb2 > be2 - be1) // at least 50% overhang
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> overhang too short (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
									stepSize);
#endif
							continue;
						}

						// check if current overlap is better (longer/more unique bases ?) then curBest

						int curUniqBasesInAIvl = (ae2 - ab2) - getRepeatBases(ctx, ovl, aread);
						int curUniqBasesInBIvl = (be2 - bb2) - getRepeatBases(ctx, ovl, bread);

						if (curBestIntersection > MAX(intersect(ab1, ae1, ab2, ae2), intersect(bb1, be1, bb2, be2)) && curBestBases < MIN(ae2 - ab2, be2 - bb2))
						{
							curBestBases = MIN(ae2 - ab2, be2 - bb2);
							curBestOffset = i - longestUniqOvlIdx;
							curBestIntersection = MAX(intersect(ab1, ae1, ab2, ae2), intersect(bb1, be1, bb2, be2));
						}

						if (curBestUniqBases < MIN(curUniqBasesInAIvl, curUniqBasesInBIvl))
						{
#ifdef DEBUG_CHAIN
							printf("found right current best overlap %d %d %d: [%d, %d] [%d, %d] right side\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							curBestUniqOffset = i - longestUniqOvlIdx;
							curBestUniqBases = MIN(curUniqBasesInAIvl, curUniqBasesInBIvl);
						}
						else if (curBestUniqBases == -1 && stepSize + STRIDE_OVL_LOOKAHEAD > MAX_OVL_LOOKAHEAD) // for repetitive genomes
						{
							Overlap *tmpOvl = ovls + (longestUniqOvlIdx + curBestOffset);

							if ((intersect(ab1, ae1, tmpOvl->path.abpos, tmpOvl->path.aepos) < ae1 - tmpOvl->path.abpos) && (intersect(bb1, be1, tmpOvl->path.bbpos, tmpOvl->path.bepos) < be1 - tmpOvl->path.bbpos))
							{
								curBestUniqOffset = curBestOffset;
								curBestUniqBases = 1;
							}
						}
						else
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> cannot be anchored (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
									stepSize);
#endif
						}
					}
				}

				// check if left best overlap can be used to extend overlap group on the right side
				if (curBestUniqBases < 0) // i.e. there was no good overlap at right side
				{
#ifdef DEBUG_CHAIN
					printf("could not extend ovlgroup on right side with proper overlap (with stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
					break;
				}

				/// todo further sanity check necessary ???
				ab2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.abpos;
				ae2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.aepos;

				bb2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.bbpos;
				be2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.bepos;
#ifdef DEBUG_CHAIN
				printf("extend ovlgroup with (right): %d %d %d: [%d, %d] [%d, %d] stepSize %d\n", longestUniqOvlIdx + curBestUniqOffset,
						ovls[longestUniqOvlIdx + curBestUniqOffset].aread, ovls[longestUniqOvlIdx + curBestUniqOffset].bread, ab1, ae1, ab2, ae2,
						stepSize - STRIDE_OVL_LOOKAHEAD);
#endif

				if (chain->novl == chain->maxOvl)
				{
					chain->maxOvl = chain->maxOvl * 1.2 + 5;
					chain->ovls = (Overlap**) realloc(chain->ovls, sizeof(Overlap*) * chain->maxOvl);
				}

				// append left side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
				chain->ovls[chain->novl] = ovls + (longestUniqOvlIdx + curBestUniqOffset);
				chain->ovls[chain->novl]->flags |= OVL_TEMP;
				chain->novl++;
				nremain--;
#ifdef DEBUG_CHAIN
				printf("chain: nOvl: %d, maxOvl %d nremain %d\n", chain->novl, chain->maxOvl, nremain);
#endif

				ab1 = ab2;
				ae1 = ae2;
				bb1 = bb2;
				be1 = be2;

				curBestUniqOffset++;
				curBestOffset = curBestUniqOffset;

				curBestUniqBases = -1;
				curBestBases = -1;

				curBestIntersection = MAX(alen, blen);

				if (longestUniqOvlIdx + curBestUniqOffset >= n)
				{
					cont = 0;
				}
			}
		}

		if (nremain && longestUniqOvlIdx > 0)
		{
			ab1 = ovls[longestUniqOvlIdx].path.abpos;
			ae1 = ovls[longestUniqOvlIdx].path.aepos;

			bb1 = ovls[longestUniqOvlIdx].path.bbpos;
			be1 = ovls[longestUniqOvlIdx].path.bepos;

#ifdef DEBUG_CHAIN
			printf("extend longest overlap in left direction\n");
#endif
// 2nd left side
			int cont = 1;
			int curBestUniqOffset = 1;
			int curBestUniqBases = -1;
			int curBestBases = -1;
			int curBestOffset = 1;
			int curBestIntersection = MAX(alen, blen);

			while (cont)
			{
				int stepSize;
				for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestUniqBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
				{
#ifdef DEBUG_CHAIN
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif

					// try to find next best overlap with lookahead of stepSize bases
					for (i = longestUniqOvlIdx - curBestUniqOffset; i >= 0; --i)
					{
						Overlap * ovl = ovls + i;

#ifdef DEBUG_CHAIN
						printf("LEFT: Check ovl: a[%d, %d] b[%d,%d]\n", ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos);
#endif

						ab2 = ovl->path.abpos;
						ae2 = ovl->path.aepos;

						bb2 = ovl->path.bbpos;
						be2 = ovl->path.bepos;

						if ((ovl->flags & OVL_COMP) != (ovls[longestUniqOvlIdx].flags & OVL_COMP))
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> different orientations\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ovl->flags & OVL_CONT)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> really contained repeat\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ovl->flags & OVL_TEMP)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> is part of another chain\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ovl->flags & OVL_DISCARD)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> discarded\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						// todo mark those as discard ????
						if (contained(ab2, ae2, ab1, ae1))
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in A-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						// todo mark those as discard ????
						if (contained(bb2, be2, bb1, be1))
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in B-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							continue;
						}

						if (ab2 > ab1 || bb2 > bb1) // also overlap must follow left direction
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> improper left extension direction (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2,
									bb2, be2, stepSize);
#endif
							continue;
						}

						if (MAX(ab1 - ae2, bb1 - be2) > stepSize)
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> gap size too large (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
									stepSize);
#endif
							continue;
						}

						if (ae2 - ab1 > ab1 - ab2 || be2 - bb1 > bb1 - bb2) // at least 50% overhang
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> overhang too short (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
									stepSize);
#endif
							continue;
						}

						// check if current overlap is better (longer/more unique bases ?) then curLeftBest

						int curUniqBasesInAIvl = (ae2 - ab2) - getRepeatBases(ctx, ovl, aread);
						int curUniqBasesInBIvl = (be2 - bb2) - getRepeatBases(ctx, ovl, bread);

						if (curBestIntersection > MAX(intersect(ab2, ae2, ab1, ae1), intersect(bb2, be2, bb1, be1)) && curBestBases < MIN(ae2 - ab2, be2 - bb2))
						{
							curBestBases = MIN(ae2 - ab2, be2 - bb2);
							curBestOffset = longestUniqOvlIdx - i;
							curBestIntersection = MAX(intersect(ab2, ae2, ab1, ae1), intersect(bb2, be2, bb1, be1));
						}

						if (curBestUniqBases < MIN(curUniqBasesInAIvl, curUniqBasesInBIvl))
						{
#ifdef DEBUG_CHAIN
							printf("found left current best overlap %d %d %d: [ab2 %d, ae2 %d] [bb2 %d, be2 %d] left side\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							curBestUniqOffset = longestUniqOvlIdx - i;
							curBestUniqBases = curUniqBasesInAIvl + curUniqBasesInBIvl;
						}
						else if (curBestUniqBases == -1 && stepSize + STRIDE_OVL_LOOKAHEAD > MAX_OVL_LOOKAHEAD) // for repetitive genomes
						{
							Overlap *tmpOvl = ovls + (longestUniqOvlIdx - curBestOffset);

							if ((intersect(tmpOvl->path.abpos, tmpOvl->path.aepos, ab1, ae1) < ae1 - tmpOvl->path.abpos) && (intersect(tmpOvl->path.bbpos, tmpOvl->path.bepos, bb1, be1) < be1 - tmpOvl->path.bbpos))
							{
								curBestUniqOffset = curBestOffset;
								curBestUniqBases = 1;
							}
						}
						else
						{
#ifdef DEBUG_CHAIN
							printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> cannot be anchored (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
									stepSize);
#endif
						}
					}
				}

				// check if left best overlap can be used to extend overlap group on the left side
				if (curBestUniqBases < 0) // i.e. there was no good overlap at left side
				{
#ifdef DEBUG_CHAIN
					printf("could not extend ovlgroup on left side with proper overlap (stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
					break;
				}

				/// todo further sanity check necessary ???
				ab2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.abpos;
				ae2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.aepos;

				bb2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.bbpos;
				be2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.bepos;

#ifdef DEBUG_CHAIN
				printf("extend ovlgroup with (left): %d %d %d: [%d, %d] [%d, %d] with stepSize %d\n", longestUniqOvlIdx - curBestUniqOffset,
						ovls[longestUniqOvlIdx - curBestUniqOffset].aread, ovls[longestUniqOvlIdx - curBestUniqOffset].bread, ab1, ae1, ab2, ae2,
						stepSize - STRIDE_OVL_LOOKAHEAD);
#endif

				if (ctx->curChains == ctx->maxChains)
				{
					ctx->maxChains = ctx->maxChains * 1.2 + 5;
					ctx->ovlChains = (Chain*) realloc(ctx->ovlChains, sizeof(Chain) * ctx->maxChains);
					bzero(ctx->ovlChains + ctx->curChains, sizeof(Chain) * (ctx->maxChains - ctx->curChains));
				}

				if (chain->novl == chain->maxOvl)
				{
					chain->maxOvl = chain->maxOvl * 1.2 + 5;
					chain->ovls = (Overlap**) realloc(chain->ovls, sizeof(Overlap*) * chain->maxOvl);
				}

				// append left side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
				chain->ovls[chain->novl] = ovls + (longestUniqOvlIdx - curBestUniqOffset);
				chain->ovls[chain->novl]->flags |= OVL_TEMP;
				chain->novl++;
				nremain--;
#ifdef DEBUG_CHAIN
				printf("chain: nOvl: %d, maxOvl %d nremain %d\n", chain->novl, chain->maxOvl, nremain);
#endif

				ab1 = ab2;
				ae1 = ae2;
				bb1 = bb2;
				be1 = be2;

				curBestUniqOffset++;
				curBestOffset = curBestUniqOffset;

				curBestUniqBases = -1;
				curBestBases = -1;

				curBestIntersection = MAX(alen, blen);

				if (longestUniqOvlIdx - curBestUniqOffset < 0)
				{
					cont = 0;
				}
			}

			if (chain->novl > 1)
			{	// sort chain
				qsort(chain->ovls, chain->novl, sizeof(Overlap*), cmp_ovls_abeg);
			}
		}

#ifdef DEBUG_CHAIN
		printf("chain: nOvl: %d, maxOvl %d nremain: %d\n", chain->novl, chain->maxOvl, nremain);
#endif

		// find possible ovls that could be added to chain (i.e. fill gaps)

		if (chain->novl > 1 /*&& nremain > 0*/)   // allow extension with previously repeat-marked overlaps
		{
#ifdef DEBUG_CHAIN
			printf("find possible ovls that could be added to chain (i.e. fill gaps)\n");
#endif
			int chainIdx = 0;
			int chainLastIdx = chain->novl - 1;
			int j;
			for (i = 0; i < n; i++)
			{
				Overlap *ovl = ovls + i;
				if ((ovl->flags & (OVL_TEMP | OVL_CONT | OVL_TRIM)) || ((ovl->flags & OVL_COMP) != (chain->ovls[chainIdx]->flags & OVL_COMP)))
					continue;

				if (ovl->path.abpos < chain->ovls[chainIdx]->path.abpos)
					continue;

				if (ovl->path.abpos > chain->ovls[chainLastIdx]->path.abpos)
					break;

				int ovhBases = 100;

				for (j = chainIdx; j < chainLastIdx; j++)
				{
					if (chain->ovls[j]->path.aepos - ovhBases < ovl->path.abpos && chain->ovls[j + 1]->path.abpos + ovhBases > ovl->path.aepos && chain->ovls[j]->path.bepos - ovhBases < ovl->path.bbpos && chain->ovls[j + 1]->path.bbpos + ovhBases > ovl->path.bepos)
					{
						Overlap *lastAddedOvl = chain->ovls[chain->novl - 1];

						if (intersect(ovl->path.abpos, ovl->path.aepos, lastAddedOvl->path.abpos, lastAddedOvl->path.aepos) > ovhBases || intersect(ovl->path.bbpos, ovl->path.bepos, lastAddedOvl->path.bbpos, lastAddedOvl->path.bepos) > ovhBases)
							break;

						if (chain->novl == chain->maxOvl)
						{
							chain->maxOvl = chain->maxOvl * 1.2 + 5;
							chain->ovls = (Overlap**) realloc(chain->ovls, sizeof(Overlap*) * chain->maxOvl);
						}

						// append left side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
						ovl->flags &= ~(OVL_DISCARD); //
						ovl->flags |= OVL_TEMP;
						chain->ovls[chain->novl] = ovl;
						chain->novl++;
						//nremain--;
#ifdef DEBUG_CHAIN
						printf("chain: nOvl: %d, maxOvl %d nremain %d\n", chain->novl, chain->maxOvl, nremain);
#endif
					}

					if (ovl->path.abpos > chain->ovls[j + 1]->path.abpos)
						chainIdx++;
				}

			}

			if (chainLastIdx < chain->novl - 1)
			{
				qsort(chain->ovls, chain->novl, sizeof(Overlap*), cmp_ovls_abeg);
			}
		}

#ifdef DEBUG_CHAIN
		printChain(chain);
#endif

		// sanity check // there should be no intersection with other chains (with same orientation) possible
		int valid = 1;
		if (ctx->curChains)
		{
#ifdef DEBUG_CHAIN
			printf("DO SANITY CHECK\n");
#endif
			int j;
			for (i = 0; i < ctx->curChains && valid; i++)
			{
				if ((chain->ovls[0]->flags & OVL_COMP) == (ctx->ovlChains[i].ovls[0]->flags & OVL_COMP))
				{
					for (j = 0; j < chain->novl; j++)
					{
						if ((chain->ovls[j]->path.abpos > ctx->ovlChains[i].ovls[0]->path.abpos && chain->ovls[j]->path.aepos < ctx->ovlChains[i].ovls[ctx->ovlChains[i].novl - 1]->path.aepos)
								|| (chain->ovls[j]->path.bbpos > ctx->ovlChains[i].ovls[0]->path.bbpos && chain->ovls[j]->path.bepos < ctx->ovlChains[i].ovls[ctx->ovlChains[i].novl - 1]->path.bepos))
						{
#ifdef DEBUG_CHAIN
							printf("CHAIN is invalid - DISCARD\n");
#endif
							valid = 0;
							break;
						}
					}
				}
			}
		}

		if (valid)
			ctx->curChains++;
		else
		{
			for (i = 0; i < chain->novl; i++)
			{
				chain->ovls[i]->flags |= OVL_DISCARD;
				ctx->statsFiltInvalidChain++;
			}

			chain->novl = 0;

		}

#ifdef DEBUG_CHAIN
		printf("curChain: %d, remain unchained OVls: %d\n", ctx->curChains, nremain);
#endif
	}

	// sort chains according to alignment lengths
	if (ctx->curChains > 1)
	{
#ifdef DEBUG_CHAIN
		printf("SORT CHAINS (longest first):\n");
#endif
		qsort(ctx->ovlChains, ctx->curChains, sizeof(Chain), cmp_chain_len);
	}
}

static int filter(FilterContext* ctx, Overlap* ovl)
{
	int ret = 0;
	
	int ovlBLen = DB_READ_LEN(ctx->db, ovl->bread);

	if (ctx->nMinNonRepeatBases != -1 && ovl->aread != ovl->bread)
	{
		int b, e, rb, re, ovllen, repeat;

		track_anno* repeats_anno = ctx->trackRepeat->anno;
		track_data* repeats_data = ctx->trackRepeat->data;

		ovllen = ovl->path.aepos - ovl->path.abpos;
		repeat = 0;

		b = repeats_anno[ovl->aread] / sizeof(track_data);
		e = repeats_anno[ovl->aread + 1] / sizeof(track_data);

		for (; b < e; b += 2)
		{
			rb = repeats_data[b];
			re = repeats_data[b + 1];

			if (rb > ovl->path.aepos)
				break;
			else if (re < ovl->path.abpos)
				continue;

			repeat += intersect(ovl->path.abpos, ovl->path.aepos, rb, re);
		}

		if (repeat > 0 && ovllen - repeat < ctx->nMinNonRepeatBases)
		{
			if (ctx->nVerbose)
			{
				printf("overlap %d -> %d: drop due to repeat in a\n", ovl->aread, ovl->bread);
			}

			ctx->statsFiltRepeat++;
			ret |= OVL_REPEAT;
		}

		if (!(ret & OVL_REPEAT))
		{
			ovllen = ovl->path.bepos - ovl->path.bbpos;

			int bbpos, bepos;

			if (ovl->flags & OVL_COMP)
			{
				bbpos = ovlBLen - ovl->path.bepos;
				bepos = ovlBLen - ovl->path.bbpos;
			}
			else
			{
				bbpos = ovl->path.bbpos;
				bepos = ovl->path.bepos;
			}

			repeat = 0;

			b = repeats_anno[ovl->bread] / sizeof(track_data);
			e = repeats_anno[ovl->bread + 1] / sizeof(track_data);
			for (; b < e; b += 2)
			{
				rb = repeats_data[b];
				re = repeats_data[b + 1];

				if (rb > bepos)
					break;
				else if (re < bbpos)
					continue;

				repeat += intersect(bbpos, bepos, rb, re);
			}

			if (repeat > 0 && ovllen - repeat < ctx->nMinNonRepeatBases)
			{
				if (ctx->nVerbose)
				{
					printf("overlap %d -> %d: drop due to repeat in b\n", ovl->aread, ovl->bread);
				}

				ctx->statsFiltRepeat++;

				ret |= OVL_REPEAT;
			}
		}
	}

	return ret;
}

static void filter_pre(PassContext* pctx, FilterContext* fctx)
{
#ifdef VERBOSE
	printf( ANSI_COLOR_GREEN "PASS filtering\n" ANSI_COLOR_RESET);

	printf( ANSI_COLOR_RED "OPTIONS\n" ANSI_COLOR_RESET);
	printf( ANSI_COLOR_RED "  keepIdentity %d\n" ANSI_COLOR_RESET, fctx->keepIdentity);
	printf( ANSI_COLOR_RED "  mergeRepDist %d\n" ANSI_COLOR_RESET, fctx->mergeRepDist);
	printf( ANSI_COLOR_RED "  nContPerc %d\n" ANSI_COLOR_RESET, fctx->nContPerc);
	printf( ANSI_COLOR_RED "  nFuzzBases %d\n" ANSI_COLOR_RESET, fctx->nFuzzBases);
	printf( ANSI_COLOR_RED "  nMinNonRepeatBases %d\n" ANSI_COLOR_RESET, fctx->nMinNonRepeatBases);
	printf( ANSI_COLOR_RED "  repeatWindowLookBack %d\n" ANSI_COLOR_RESET, fctx->repeatWindowLookBack);
      
        if(fctx->trackRepeat)
		printf( ANSI_COLOR_RED "  RepeatTrack %s\n" ANSI_COLOR_RESET, fctx->trackRepeat->name);
	if(fctx->trackTrim)
                printf( ANSI_COLOR_RED "  TrimTrack %s\n" ANSI_COLOR_RESET, fctx->trackTrim->name);
	if(fctx->trackLowCompl)
                printf( ANSI_COLOR_RED "  LowComplexityTrack %s\n" ANSI_COLOR_RESET, fctx->trackLowCompl->name);
#endif

	fctx->twidth = pctx->twidth;

	fctx->curChains = 0;
	fctx->maxChains = 5;
	fctx->ovlChains = (Chain*) malloc(sizeof(Chain) * MAX(fctx->maxChains, fctx->nkeptChains));
	bzero(fctx->ovlChains, sizeof(Chain) * MAX(fctx->maxChains, fctx->nkeptChains));

	fctx->maxUniqAIntervals = 20;
	fctx->uniqAIntervals = malloc(sizeof(anchorItv) * fctx->maxUniqAIntervals);
	bzero(fctx->uniqAIntervals, sizeof(anchorItv) * fctx->maxUniqAIntervals);
	fctx->curUniqAIntervals = 0;

	fctx->maxUniqBIntervals = 20;
	fctx->uniqBIntervals = malloc(sizeof(anchorItv) * fctx->maxUniqBIntervals);
	bzero(fctx->uniqBIntervals, sizeof(anchorItv) * fctx->maxUniqBIntervals);
	fctx->curUniqBIntervals = 0;

	if(fctx->trackTrim)
		fctx->trim = trim_init(fctx->db, pctx->twidth, fctx->trackTrim, NULL);
}

static void filter_post(FilterContext* ctx)
{
#ifdef VERBOSE

	if (ctx->statsFiltRepeat > 0)
	{
		printf("min non-repeat bases of %4d discarded     %10d\n", ctx->nMinNonRepeatBases, ctx->statsFiltRepeat);
	}

	if (ctx->statsFiltContained > 0)
	{
		printf("contained LAS discarded     %10d\n", ctx->statsFiltContained);
	}

	if (ctx->statsFiltInvalidChain > 0)
	{
		printf("invalid LAS chains discarded     %10d\n", ctx->statsFiltInvalidChain);
	}

	if (ctx->statsGapAlns > 0)
	{
		printf("gap LAS discarded     %10d\n", ctx->statsGapAlns);
	}

	if (ctx->statsLowCovALn > 0)
	{
		printf("low coverage LAS discarded     %10d\n", ctx->statsLowCovALn);
	}

#endif

	int i;
	for (i = 0; i < ctx->maxChains; i++)
	{
		Chain *chain = ctx->ovlChains + i;
		if (chain)
			free(chain->ovls);
		else
			break;
	}
	free(ctx->ovlChains);
	free(ctx->uniqAIntervals);
	free(ctx->uniqBIntervals);

	if(ctx->trackTrim)
		trim_close(ctx->trim);
}

static void getRepeatBasesFromInterval(HITS_TRACK* repeat, int readID, int beg, int end, int *cumBases, int *largest)
{
	track_anno* rep_anno = repeat->anno;
	track_data* rep_data = repeat->data;

	track_anno rb, re;

	*cumBases = 0;
	*largest = 0;
	int rBeg, rEnd, tmp;

	// repeat bases in a-read
	rb = rep_anno[readID] / sizeof(track_data);
	re = rep_anno[readID + 1] / sizeof(track_data);

	while (rb < re)
	{
		rBeg = rep_data[rb];
		rEnd = rep_data[rb + 1];

		tmp = intersect(beg, end, rBeg, rEnd);

		if (tmp)
		{
			*cumBases += tmp;
			if (*largest < (rEnd - rBeg))
				*largest = (rEnd - rBeg);
		}

		rb += 2;
	}
}

static int cmp_aIvl(const void *a, const void *b)
{
	anchorItv * a1 = (anchorItv*) a;
	anchorItv * a2 = (anchorItv*) b;

	if (a1->flag & ANCHOR_INVALID)
	{
		return 1;
	}
	if (a2->flag & ANCHOR_INVALID)
	{
		return -1;
	}
	return a1->beg - a2->beg;
}

static void createUniqueMask(FilterContext *ctx, int read, int isAread)
{
	int *numIntervals;
	int *curItv;
	anchorItv *uniqIntervals;

	if (isAread)
	{
		numIntervals = &(ctx->maxUniqAIntervals);
		curItv = &(ctx->curUniqAIntervals);
		uniqIntervals = ctx->uniqAIntervals;
	}
	else
	{
		numIntervals = &(ctx->maxUniqBIntervals);
		curItv = &(ctx->curUniqBIntervals);
		uniqIntervals = ctx->uniqBIntervals;
	}

	int trim_beg, trim_end;
	int rlen = DB_READ_LEN(ctx->db, read);

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, read, &trim_beg, &trim_end);
	}
	else
	{
		trim_beg = 0;
		trim_end = rlen;
	}

#ifdef DEBUG_CHAIN
	printf("trim_beg %d, trim_end %d\n", trim_beg, trim_end);
#endif

	if(ctx->trackRepeat == NULL)
	{
		ctx->curUniqAIntervals=1;
		ctx->uniqAIntervals[0].beg = trim_beg;
		ctx->uniqAIntervals[0].end = trim_end;
		return;
	}

	int MINANCHOR = 10;
	int WINDOW = ctx->repeatWindowLookBack;
	int MAXMERGE = ctx->mergeRepDist;
	int i, b, e;

	track_anno* repeats_anno = ctx->trackRepeat->anno;
	track_data* repeats_data = ctx->trackRepeat->data;

	b = repeats_anno[read] / sizeof(track_data);
	e = repeats_anno[read + 1] / sizeof(track_data);

	if (*numIntervals < (e - b + 1) + 4)
	{
		*numIntervals = (e - b + 1) + 4;
		if (isAread)
		{
			ctx->uniqAIntervals = (anchorItv*) realloc(ctx->uniqAIntervals, *numIntervals * sizeof(anchorItv));
			uniqIntervals = ctx->uniqAIntervals;
		}
		else
		{
			ctx->uniqBIntervals = (anchorItv*) realloc(ctx->uniqBIntervals, *numIntervals * sizeof(anchorItv));
			uniqIntervals = ctx->uniqBIntervals;
		}
	}

	// reset current anchor interval index
	*curItv = 0;
	bzero(uniqIntervals, sizeof(anchorItv) * (*numIntervals));
	int anchorbases = 0;

	if (b < e)
	{
		int rb1, rb2;
		int re1, re2;

		rb1 = repeats_data[b];
		re1 = repeats_data[b + 1];

		if (rb1 > 0 && rb1 > MINANCHOR)
		{
			uniqIntervals[*curItv].beg = 0;
			uniqIntervals[*curItv].end = rb1;
			(*curItv)++;
		}

		b += 2;
		while (b < e)
		{
			rb2 = repeats_data[b];
			re2 = repeats_data[b + 1];

			if (rb2 - re1 > MINANCHOR)
			{
				uniqIntervals[*curItv].beg = re1;
				uniqIntervals[*curItv].end = rb2;
				(*curItv)++;
			}

			rb1 = rb2;
			re1 = re2;
			b += 2;
		}

		if (re1 < rlen && rlen - re1 > MINANCHOR)
		{
			uniqIntervals[*curItv].beg = re1;
			uniqIntervals[*curItv].end = rlen;
			(*curItv)++;
		}

		anchorbases = 0;
		for (i = 0; i < *curItv; i++)
		{
			anchorItv *a = uniqIntervals + i;
			if (a->flag & ANCHOR_INVALID)
				continue;

			anchorbases += a->end - a->beg;
		}

		// update unique intervals based on trim track
		if (trim_beg > 0 || trim_end < rlen)
		{
			for (i = 0; i < *curItv; i++)
			{
				anchorItv *a = uniqIntervals + i;

				if (trim_beg >= a->end)
				{
					a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
				}
				else if (trim_beg > a->beg)
				{
					a->beg = trim_beg;
				}

				if (a->beg >= trim_end)
				{
					a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
				}
				else if (a->end > trim_end)
				{
					a->end = trim_end;
				}
			}
		}

		anchorbases = 0;
		for (i = 0; i < *curItv; i++)
		{
			anchorItv *a = uniqIntervals + i;
			if (a->flag & ANCHOR_INVALID)
				continue;

			anchorbases += a->end - a->beg;
		}

		// update unique intervals based on low complexity and tandem repeat
		// todo hardcoded values !!!
		int predust, dust, postdust, longestDust, longestDustl, longestDustr;
		for (i = 0; i < *curItv; i++)
		{
			anchorItv *a = uniqIntervals + i;

			if (a->flag & ANCHOR_INVALID)
				continue;

			if (a->end - a->beg > MAXMERGE)
				continue;

			getRepeatBasesFromInterval(ctx->trackLowCompl, read, a->beg, a->end, &dust, &longestDust);

			if (dust * 100.0 / (a->end - a->beg) > 50.0)
			{
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
			}
			else if (dust * 100.0 / (a->end - a->beg) > 15.0 && longestDust > 100)
			{
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
			}
			else if ((a->end - a->beg) < 100 && longestDust > 20)
			{
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
			}
			else // check if neighboring repeats end in low complexity interval
			{
				int checkFlanks = MIN(100, WINDOW);

				while (checkFlanks <= WINDOW)
				{
					getRepeatBasesFromInterval(ctx->trackLowCompl, read, MAX(0, a->beg - checkFlanks), a->beg, &predust, &longestDustl);
					getRepeatBasesFromInterval(ctx->trackLowCompl, read, a->end, MIN(a->end + checkFlanks, rlen), &postdust, &longestDustr);

					if ((predust * 100.0 / checkFlanks > 20.0 && longestDustl > 30) || (postdust * 100 / checkFlanks > 20.0 && longestDustr > 30))
					{
						a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
						break;
					}

					if (checkFlanks < WINDOW && checkFlanks + 100 > WINDOW)
						checkFlanks = WINDOW;
					else
						checkFlanks += 100;
				}
			}
#ifdef DEBUG_CHAIN
			printf("#LC %d %d %d f%d PRE %d %d %.2f DUST %d %d %.2f post %d %d %.2f SUM %d %d %.2f\n", read, a->beg, a->end, a->flag, predust,
					a->beg - MAX(0, a->beg - WINDOW), predust * 100.0 / (a->beg - MAX(0, a->beg - WINDOW)), dust, a->end - a->beg, dust * 100.0 / (a->end - a->beg),
					postdust, MIN(a->end + WINDOW, rlen) - a->end, postdust * 100.0 / (MIN(a->end + WINDOW, rlen) - a->end), predust + dust + postdust,
					(a->beg - MAX(0, a->beg - WINDOW)) + (a->end - a->beg) + (MIN(a->end + WINDOW, rlen) - a->end),
					(predust + dust + postdust) * 100.0 / ((a->beg - MAX(0, a->beg - WINDOW)) + (a->end - a->beg) + (MIN(a->end + WINDOW, rlen) - a->end)));
#endif
		}
	}
	else // add full read interval as unique range
	{
		uniqIntervals[0].beg = trim_beg;
		uniqIntervals[0].end = trim_end;
		(*curItv)++;
	}

	anchorbases = 0;
	for (i = 0; i < *curItv; i++)
	{
		anchorItv *a = uniqIntervals + i;
		if (a->flag & ANCHOR_INVALID)
			continue;

		anchorbases += a->end - a->beg;
	}

	repeats_anno = ctx->trackLowCompl->anno;
	repeats_data = ctx->trackLowCompl->data;

	b = repeats_anno[read] / sizeof(track_data);
	e = repeats_anno[read + 1] / sizeof(track_data);

	int rb, re;

	// update unique anchors with all low complexity intervals !!!
	int c = *curItv;
	for (i = 0; i < c; i++)
	{
		if (uniqIntervals[i].flag & ANCHOR_INVALID)
			continue;

		while (b < e)
		{
			rb = repeats_data[b];
			re = repeats_data[b + 1];

			if (rb > uniqIntervals[i].end)
			{
				break;
			}

			if (re < uniqIntervals[i].beg)
			{
				b += 2;
				continue;
			}

			// dust fully covers unique part
			if (rb <= uniqIntervals[i].beg && re >= uniqIntervals[i].end)
			{
				uniqIntervals[i].flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
				break;
			}

			// dust aligns left with unique part
			if (rb <= uniqIntervals[i].beg)
			{
				uniqIntervals[i].beg = re;
			}
			// dust aligns with right unique part
			if (re >= uniqIntervals[i].end)
			{
				uniqIntervals[i].end = rb;
			}
			// dust splits uniq part, i.e. make unique part invalid an append splits to the end of uniqueIntervals
			if (*curItv >= *numIntervals)
			{
				*numIntervals = 1.2 * (*numIntervals) + 10;
				if (isAread)
				{
					ctx->uniqAIntervals = (anchorItv*) realloc(ctx->uniqAIntervals, *numIntervals * sizeof(anchorItv));
					bzero(ctx->uniqAIntervals + (*curItv), sizeof(anchorItv) * (*numIntervals - *curItv));
					uniqIntervals = ctx->uniqAIntervals;
				}
				else
				{
					ctx->uniqBIntervals = (anchorItv*) realloc(ctx->uniqBIntervals, *numIntervals * sizeof(anchorItv));
					bzero(ctx->uniqBIntervals + (*curItv), sizeof(anchorItv) * (*numIntervals - *curItv));
					uniqIntervals = ctx->uniqBIntervals;
				}
			}

			uniqIntervals[*curItv].beg = uniqIntervals[i].beg;
			uniqIntervals[*curItv].end = rb;

			uniqIntervals[i].beg = re;
			(*curItv)++;

			b += 2;
		}
	}

	qsort(uniqIntervals, *curItv, sizeof(anchorItv), cmp_aIvl);
	for (i = 0; i < *curItv; i++)
	{
		anchorItv *a = uniqIntervals + i;
		if (a->flag & ANCHOR_INVALID)
			break;
	}
	*curItv = i;

	// merge tips if required, i.e. if there is any repeat annotation within the first/last 2k?! sequence
	if (ctx->rp_mergeTips && *curItv > 0)
	{
		int resort = 0;
		if (uniqIntervals[0].beg > trim_beg || uniqIntervals[0].end < trim_beg + ctx->rp_mergeTips)
		{
			for (i = 0; i < *curItv; i++)
			{
				anchorItv *a = uniqIntervals + i;

				if (a->flag & ANCHOR_INVALID)
					continue;

				if (a->end < trim_beg + ctx->rp_mergeTips)
				{
					a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
					resort = 1;
				}
				else if (a->beg < trim_beg + ctx->rp_mergeTips)
				{
					a->beg = trim_beg + ctx->rp_mergeTips;
				}
				else
				{
					break;
				}
			}
		}

		if (uniqIntervals[*curItv - 1].end < trim_end || uniqIntervals[*curItv - 1].beg > trim_end - ctx->rp_mergeTips)
		{
			for (i = *curItv - 1; i >= 0; --i)
			{
				anchorItv *a = uniqIntervals + i;

				if (a->flag & ANCHOR_INVALID)
					continue;

				if (a->beg > trim_end - ctx->rp_mergeTips)
				{
					a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
					resort = 1;
				}
				else if (a->end > trim_end - ctx->rp_mergeTips)
				{
					a->end = trim_end - ctx->rp_mergeTips;
				}
				else
				{
					break;
				}
			}
		}

		if (resort)
		{
			qsort(uniqIntervals, *curItv, sizeof(anchorItv), cmp_aIvl);
			for (i = 0; i < *curItv; i++)
			{
				anchorItv *a = uniqIntervals + i;
				if (a->flag & ANCHOR_INVALID)
					break;
			}
			*curItv = i;
		}
	}

#ifdef DEBUG_CHAIN
	// report final unique anchors:
	anchorbases = 0;
	printf("#final anchors %d", read);
	for (i = 0; i < *curItv; i++)
	{
		anchorItv *a = uniqIntervals + i;
		if (a->flag & ANCHOR_INVALID)
		break;

		printf(" %d-%d-%d", a->beg, a->end, a->flag);
		anchorbases += a->end - a->beg;
	}
	printf(" sum n%d b%d\n", i, anchorbases);
#endif
}

static int gapIsLowComplexity(FilterContext *ctx, int read, int beg, int end, float fraction)
{
	if (!ctx->trackLowCompl)
		return 0;

	int longest;
	int cumLC;

	getRepeatBasesFromInterval(ctx->trackLowCompl, read, beg, end, &cumLC, &longest);

	if ((cumLC * 1.0) / (end - beg) < fraction)
		return 0;

	return 1;
}

static int stitchChain(FilterContext *ctx, Chain *chain)
{
	int stitched = 0;

	if (chain->novl < 2)
	{
		return stitched;
	}

	int i;
	int ab2, ae1, ae2;
	int bb2, be1, be2;

	int aread = chain->ovls[0]->aread;
	int bread = chain->ovls[0]->bread;

	// check if any overlap is contained in low complexity interval, if so then don't stitch at all
	if (ctx->trackLowCompl)
	{
		for (i = 0; i < chain->novl; i++)
		{
			Overlap *ovl = chain->ovls[i];
			int repeatBases;
			int longestRep;

			// check LowComp of A-read
			getRepeatBasesFromInterval(ctx->trackLowCompl, aread, ovl->path.abpos, ovl->path.aepos, &repeatBases, &longestRep);
			if ((ovl->path.aepos - ovl->path.abpos) - repeatBases < ctx->stitchNonLowCompAnchorBases)
			{
				return stitched;
			}

			// check LowComp of B-read
			if (ovl->flags & OVL_COMP)
			{
				getRepeatBasesFromInterval(ctx->trackLowCompl, bread, DB_READ_LEN(ctx->db, bread) - ovl->path.bepos, DB_READ_LEN(ctx->db, bread) - ovl->path.bbpos, &repeatBases, &longestRep);
			}
			else
			{
				getRepeatBasesFromInterval(ctx->trackLowCompl, bread, ovl->path.bbpos, ovl->path.bepos, &repeatBases, &longestRep);
			}
			if ((ovl->path.bepos - ovl->path.bbpos) - repeatBases < ctx->stitchNonLowCompAnchorBases)
			{
				return 0;
			}
		}
	}

	Overlap* ovli = chain->ovls[0];
	ae1 = ovli->path.aepos;
	be1 = ovli->path.bepos;

	for (i = 1; i < chain->novl; i++)
	{
		Overlap* ovlk = chain->ovls[i];
		assert((ovli->flags & OVL_COMP) == (ovlk->flags & OVL_COMP));

		ab2 = ovlk->path.abpos;
		ae2 = ovlk->path.aepos;

		bb2 = ovlk->path.bbpos;
		be2 = ovlk->path.bepos;

		int deltaa = abs(ae1 - ab2);
		int deltab = abs(be1 - bb2);

		float LOWCOMPFRACTION = 0.8;

		int isLowCompl = 0;
		int cumLowCompBases = 0;
		int longestLowCompBases = 0;
		if(ae1 <= ab2)
		{
			getRepeatBasesFromInterval(ctx->trackLowCompl, ovli->aread, ae1, ab2, &cumLowCompBases, &longestLowCompBases);
		}
		else
		{
			getRepeatBasesFromInterval(ctx->trackLowCompl, ovli->aread, ab2, ae1, &cumLowCompBases, &longestLowCompBases);
		}

		if(deltaa*LOWCOMPFRACTION <= cumLowCompBases)
		{
			isLowCompl = 1;
		}
		else	// check b-read
		{
			cumLowCompBases = 0;
			longestLowCompBases = 0;
			if(ovli->flags & OVL_COMP)
			{
				assert(ovlk->flags & OVL_COMP);
				int cbe1 = DB_READ_LEN(ctx->db, ovli->bread) - be1;
				int cbb2 = DB_READ_LEN(ctx->db, ovli->bread) - bb2;

				if(cbe1 <= cbb2)
				{
					getRepeatBasesFromInterval(ctx->trackLowCompl, ovli->bread, cbe1, cbb2, &cumLowCompBases, &longestLowCompBases);
				}
				else
				{
					getRepeatBasesFromInterval(ctx->trackLowCompl, ovli->bread, cbb2, cbe1, &cumLowCompBases, &longestLowCompBases);
				}
			}
			else
			{
				if(be1 <= bb2)
				{
					getRepeatBasesFromInterval(ctx->trackLowCompl, ovli->bread, be1, bb2, &cumLowCompBases, &longestLowCompBases);
				}
				else
				{
					getRepeatBasesFromInterval(ctx->trackLowCompl, ovli->bread, bb2, be1, &cumLowCompBases, &longestLowCompBases);
				}
			}

			if(deltab*LOWCOMPFRACTION <= cumLowCompBases)
			{
				isLowCompl = 1;
			}
		}

		if (ctx->stitchMaxGapSize < 0 || (deltaa < ctx->stitchMaxGapSize && deltab < ctx->stitchMaxGapSize)
				|| (isLowCompl && (deltaa < ctx->stitchMaxGapSizeInLowCompl && deltab < ctx->stitchMaxGapSizeInLowCompl)))
		{
#ifdef VERBOSE_STITCH
			int ab1 = ovli->path.abpos;
			int bb1 = ovli->path.bbpos;

			printf("STITCH %8d @ %5d..%5d -> %8d @ %5d..%5d %c\n"
					"                  %5d..%5d -> %8d @ %5d..%5d %c\n",
					ovli->aread,
					ab1, ae1, ovli->bread, bb1, be1, OVL_STRAND(ovli),
					ab2, ae2, ovlk->bread, bb2, be2, OVL_STRAND(ovlk)));
#endif

			ovli->path.aepos = ae2;
			ovli->path.bepos = be2;
			ovli->path.diffs += ovlk->path.diffs;
			ovli->path.tlen = 0;

			ae1 = ae2;
			be1 = be2;

			assert(ovli->bread == ovlk->bread);

			ovli->flags &= ~( OVL_DISCARD | OVL_LOCAL); // force a re-evaluation of the OVL_LOCAL flags

			ovlk->flags |= OVL_DISCARD | OVL_STITCH;

			stitched += 1;

#ifdef VERBOSE_STITCH
			printf( "    -> %8d @ %5d..%5d -> %8d @ %5d..%5d %c  delta a %3d b %3d\n",
					ovli->aread,
					ovli->path.abpos, ovli->path.aepos,
					ovli->bread,
					ovli->path.abpos, ovli->path.aepos,
					OVL_STRAND( ovli ),
					deltaa, deltab );
#endif
		}
		else
		{
			ovli = ovlk;
		}
	}

	return stitched;
}

static void findGaps(FilterContext *ctx, Overlap *ovl, int novl)
{

	// todo hard coded: use a dynamic window dependent on repeat track ?
	int SWINDOW = 1000;
	int ANCHOR = MAX(ctx->minChainLen / 2, 1 + SWINDOW / 2);
	int MINSPANNER = ctx->gapMinSpanners;

	int trim_abeg, trim_aend;

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trim_abeg, &trim_aend);
	}
	else
	{
		trim_abeg = 0;
		trim_aend = DB_READ_LEN(ctx->db, ovl->aread);
	}

	int foundGap = 0;
	int i, j, k, l;
#ifdef DEBUG_GAPS
	printf("find gaps: ");
#endif
	int count = 0;
	for (i = 0; i < novl; i++)
	{
		if (!(ovl[i].flags & OVL_DISCARD))
		{
			count++;
		}
	}
	if(count == 0)
		return;

#ifdef DEBUG_GAPS
	printf("ovl: %d/%d\n", count, novl);
#endif

	for (i = trim_abeg; i < trim_aend && foundGap == 0; i += SWINDOW)
	{
		int nspanner = 0;

		j = k = 0;
		while (j < novl)
		{
			while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
			{
				k++;
			}

			if (ovl[j].aread != ovl[j].bread)
			{
				int abpos = DB_READ_LEN(ctx->db, ovl[j].aread);
				int aepos = 0;

				for (l = j; l <= k; l++)
				{
					Overlap * o = ovl + l;
					if (o->flags & OVL_DISCARD)
						continue;

					if (o->path.abpos < abpos)
						abpos = o->path.abpos;

					if (o->path.aepos > aepos)
						aepos = o->path.aepos;
#ifdef DEBUG_GAPS
					printf("t|%d,%d| i|%d,%d| [%d, %d] {%d %d} (%d %d)\n", trim_abeg, trim_aend, i - ANCHOR, i + ANCHOR, o->aread, o->bread, abpos, aepos, o->path.abpos, o->path.aepos);
#endif
				}

				if (abpos < aepos)
				{
					if ((abpos <= MAX(trim_abeg, i - ANCHOR)) && (aepos >= MIN(i + ANCHOR, trim_aend)))
					{
						nspanner++;
					}
#ifdef DEBUG_GAPS
					printf(" ---> nspanner: %d\n", nspanner);
#endif
				}

				if(nspanner >= MINSPANNER)
					break;
			}
			k++;
			j = k;
		}

		if (nspanner < MINSPANNER)
		{

			printf("FOUND GAP in READ %d in range [%d, %d]\n", ovl->aread, MAX(trim_abeg, i - ANCHOR), MIN(i + ANCHOR, trim_aend));

			for (j = 0; j < novl; j++)
			{
				ovl[j].flags |= (OVL_DISCARD | OVL_GAP);
				ctx->statsGapAlns++;
			}
			fprintf(ctx->fileOutDiscardedReads, "%d GAP\n", ovl->aread);
			return;
		}
	}
}

static void checkTipCoverage(FilterContext *ctx, Overlap *ovl, int novl)
{

	int trimBeg, trimEnd;
	trimBeg = 0;
	trimEnd = DB_READ_LEN(ctx->db, ovl->aread);

	if (ctx->trackTrim)
		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trimBeg, &trimEnd);

	int i;
	int count = 0;
	for (i = 0; i < novl; i++)
	{
		if (!(ovl[i].flags & OVL_DISCARD))
		{
			count++;
		}
	}

	if(count == 0)
		return;

	if (trimEnd - trimBeg)
	{
		int entercov = 0;
		int leavecov = 0;
		int bases = 0;

		char * cov_read_active = malloc(DB_READ_MAXLEN(ctx->db));
		bzero(cov_read_active, DB_READ_MAXLEN(ctx->db));

		int j;
		for (j = 0; j < novl; j++)
		{
			Overlap* ovl_j = ovl + j;

			if (ovl_j->flags & OVL_DISCARD)
				continue;

			if (ovl_j->path.abpos <= trimBeg)
				entercov++;

			if (ovl_j->path.aepos >= trimEnd)
				leavecov++;

			bases += ovl_j->path.aepos - ovl_j->path.abpos;
			memset(cov_read_active + ovl_j->path.abpos, 1, ovl_j->path.aepos - ovl_j->path.abpos);
		}

		int active = 0;
		for (j = trimBeg; j < trimEnd; j++)
		{
			active += cov_read_active[j];
		}

		if (bases / (trimEnd - trimBeg) <= ctx->minTipCoverage || (leavecov < ctx->minTipCoverage || entercov < ctx->minTipCoverage) || ((trimEnd - trimBeg) - active > ctx->minTipCoverage))
		{
			for (j = 0; j < novl; j++)
			{
				Overlap* ovl_j = ovl + j;

				if (ovl_j->flags & OVL_DISCARD)
					continue;

				ovl_j->flags |= OVL_DISCARD;
				ctx->statsLowCovALn++;
			}
			printf("DROP LOWCOV AREAD %d (b: %d, e: %d, avgCov %d, gapBases: %d)\n", ovl->aread, entercov, leavecov, bases / (trimEnd - trimBeg), (trimEnd - trimBeg) - active);
			fprintf(ctx->fileOutDiscardedReads, "%d LCOV\n", ovl->aread);
		}
		free(cov_read_active);
	}
}

static int filter_handler(void* _ctx, Overlap* ovl, int novl)
{
	FilterContext* ctx = (FilterContext*) _ctx;
	int i, j, k;

	int trim_abeg, trim_aend;
	int trim_bbeg, trim_bend;

	if (ctx->trackTrim)
	{
		for (j = 0; j < novl; j++)
		{
			trim_overlap(ctx->trim, ovl + j);
//			if(ovl[j].flags & OVL_TRIM)
//			printf("TRIMMED %d %d [%d, %d] [%d,%d]\n", ovl[j].aread, ovl[j].bread, ovl[j].path.abpos, ovl[j].path.aepos, ovl[j].path.bbpos, ovl[j].path.bepos);
		}
//		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trim_abeg, &trim_aend);
//		printf("T[%d, %d]\n", trim_abeg, trim_aend);
	}

// set filter flags
	for (j = 0; j < novl; j++)
	{
		// get rid of all previous flags
		ovl[j].flags &= ~(OVL_CONT | OVL_TEMP | OVL_REPEAT);
		ovl[j].flags |= filter(ctx, ovl + j);
	}

	j = k = 0;

// create unique mask for a-read
	createUniqueMask(ctx, ovl->aread, 1);

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trim_abeg, &trim_aend);
	}
	else
	{
		trim_abeg = 0;
		trim_aend = DB_READ_LEN(ctx->db, ovl->aread);
	}

	while (j < novl)
	{

		//printf("%d vs %d\n", ovl[j].aread, ovl[j].bread); 
		assert(k == j);

		int nAnchorOvls = (ovl[j].flags & (OVL_REPEAT | OVL_TRIM)) ? 0 : 1;
		//printf("---< nAnchorOvls: %d a(%d, %d) b(%d,%d) R? %d, T? %d\n", nAnchorOvls,ovl[j].path.abpos,ovl[j].path.aepos,ovl[j].path.bbpos,ovl[j].path.bepos, (ovl[j].flags & (OVL_REPEAT)), (ovl[j].flags & (OVL_TRIM)));
		while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
		{
			k++;
			nAnchorOvls += (ovl[k].flags & (OVL_REPEAT | OVL_TRIM)) ? 0 : 1;
			//printf("---< nAnchorOvls: %d a(%d, %d) b(%d,%d) R? %d, T? %d\n", nAnchorOvls, ovl[k].path.abpos,ovl[k].path.aepos,ovl[k].path.bbpos,ovl[k].path.bepos, (ovl[j].flags & (OVL_REPEAT)), (ovl[j].flags & (OVL_TRIM)));
		}
#ifdef CHAIN_DEBUG
		printf("AnchorOvls: %d k: %d vs %d j: %d vs %d REP: %d, TRIM: %d\n", nAnchorOvls, ovl[j].aread, ovl[j].bread, ovl[k].aread, ovl[k].bread, ovl[j].flags & (OVL_REPEAT), ovl[j].flags & (OVL_TRIM));
#endif
		// ignore all self alignments if those are present
		if (ovl[j].aread == ovl[j].bread)
			nAnchorOvls = 0;

		if (nAnchorOvls)
		{
#ifdef CHAIN_DEBUG
			printf("read: %8d len(%10d) | read: %8d len(%10d) novl: %10d, anchorOvl: %5d\n", ovl[j].aread, DB_READ_LEN(ctx->db, ovl[j].aread), ovl[j].bread, DB_READ_LEN(ctx->db, ovl[j].bread), k - j + 1, nAnchorOvls);
			fflush(stdout);
#endif
			chain(ctx, ovl + j, k - j + 1);
#ifdef CHAIN_DEBUG
			printf("FINAL CHAINS: %d %7d vs %7d\n", ctx->curChains, ctx->ovlChains[0].ovls[0]->aread, ctx->ovlChains[0].ovls[0]->bread);

			for (i = 0; i < ctx->curChains; i++)
			{
				printf(" CHAIN %d/%d: #novl %d\n", i + 1, ctx->curChains, ctx->ovlChains[i].novl);

				int j;
				for (j = 0; j < ctx->ovlChains[i].novl; j++)
				{
					printf("  OVL %d/%d: a[%7d, %7d] b[%7d, %7d] %s\n", j + 1, ctx->ovlChains[i].novl, ctx->ovlChains[i].ovls[j]->path.abpos, ctx->ovlChains[i].ovls[j]->path.aepos, ctx->ovlChains[i].ovls[j]->path.bbpos, ctx->ovlChains[i].ovls[j]->path.bepos,
							(ctx->ovlChains[i].ovls[j]->flags & OVL_COMP) ? "COMP" : "NORM");
				}
			}
#endif
			int a, b;

			// discard all overlaps, that are not part of a valid chain
			for (i = j; i <= k; i++)
			{
				if (!(ovl[i].flags & OVL_TEMP))
				{
					ovl[i].flags |= OVL_DISCARD;
					ctx->statsFiltInvalidChain++;
				}
			}

			if (ctx->nkeptChains == 0)
			{
				for (a = 1; a < ctx->curChains; a++)
					for (b = 0; b < ctx->ovlChains[a].novl; b++)
					{
						ctx->ovlChains[a].ovls[b]->flags |= OVL_DISCARD;
						ctx->statsFiltInvalidChain++;
					}
			}

			{
				// create unique mask for b-read
				// printf("-----> create unique mask for b_read %d\n", ovl[j].bread);
				createUniqueMask(ctx, ovl[j].bread, 0);
				//printf("ctx->curUniqBIntervals: %d\n", ctx->curUniqBIntervals);

				if (ctx->trackTrim)
				{
					get_trim(ctx->db, ctx->trackTrim, ovl[j].bread, &trim_bbeg, &trim_bend);
				}
				else
				{
					trim_bbeg = 0;
					trim_bend = DB_READ_LEN(ctx->db, ovl[j].bread);
				}

				// check for proper chain
				for (a = 0; a < ctx->curChains; a++)
				{
					if (ctx->nkeptChains == 0 && a > 0) // do not check, they are already discarded
						break;

					Chain *chain = ctx->ovlChains + a;

					int properBegA = 0;
					int properEndA = 0;
					int properBegB = 0;
					int properEndB = 0;
					int coveredBasesInAread = 0;
					int coveredBasesInBread = 0;
					int properGapLen = 1;

					int alignedBasesInAread;
					int alignedBasesInBread;
					int numdiffs = 0;

					int fuzzy = ctx->nFuzzBases;

					// check for proper begin
					if (chain->ovls[0]->path.abpos <= trim_abeg + fuzzy)
						properBegA = 1;

					if (chain->ovls[0]->flags & OVL_COMP)
					{
						if (DB_READ_LEN(ctx->db,chain->ovls[0]->bread) - chain->ovls[0]->path.bbpos + fuzzy >= trim_bend)
							properBegB = 1;
					}
					else
					{
						if (chain->ovls[0]->path.bbpos <= trim_bbeg + fuzzy)
							properBegB = 1;
					}

					// check for proper end
					if (chain->ovls[chain->novl - 1]->path.aepos + fuzzy >= trim_aend)
						properEndA = 1;

					if (chain->ovls[0]->flags & OVL_COMP)
					{
						if (DB_READ_LEN(ctx->db,chain->ovls[0]->bread) - chain->ovls[chain->novl - 1]->path.bepos - fuzzy <= trim_bbeg)
							properEndB = 1;
					}
					else
					{
						if (chain->ovls[chain->novl - 1]->path.bepos + fuzzy >= trim_bend)
							properEndB = 1;
					}

					coveredBasesInAread = chain->ovls[0]->path.aepos - chain->ovls[0]->path.abpos;
					coveredBasesInBread = chain->ovls[0]->path.bepos - chain->ovls[0]->path.bbpos;

					alignedBasesInAread = chain->ovls[0]->path.aepos - chain->ovls[0]->path.abpos;
					alignedBasesInBread = chain->ovls[0]->path.bepos - chain->ovls[0]->path.bbpos;

					numdiffs = chain->ovls[0]->path.diffs;

					for (b = 1; b < chain->novl; b++)
					{
						coveredBasesInAread += chain->ovls[b]->path.aepos - chain->ovls[b]->path.abpos;
						coveredBasesInBread += chain->ovls[b]->path.bepos - chain->ovls[b]->path.bbpos;

						alignedBasesInAread += chain->ovls[b]->path.aepos - chain->ovls[b]->path.abpos;
						alignedBasesInBread += chain->ovls[b]->path.bepos - chain->ovls[b]->path.bbpos;

						numdiffs += chain->ovls[b]->path.diffs;

						// check for intersection in A
						if (chain->ovls[b]->path.abpos < chain->ovls[b - 1]->path.aepos)
						{
							coveredBasesInAread -= chain->ovls[b - 1]->path.aepos - chain->ovls[b]->path.abpos;
						}
						// check for gap in A
						else
						{
							if (chain->ovls[b]->path.abpos - chain->ovls[b - 1]->path.aepos > fuzzy)
							{
								// check gap
								//printf("gap [%d, %d] in %d is low complexity: %d\n", chain->ovls[b - 1]->path.aepos, chain->ovls[b]->path.abpos,chain->ovls[b]->aread, gapIsLowComplexity(ctx, chain->ovls[b]->aread, chain->ovls[b - 1]->path.aepos, chain->ovls[b]->path.abpos, 0.8),);
								if (!gapIsLowComplexity(ctx, chain->ovls[b]->aread, chain->ovls[b - 1]->path.aepos, chain->ovls[b]->path.abpos, 0.8))
								{
									properGapLen = 0;
									break;
								}
							}
						}
						// check for intersection in B
						if (chain->ovls[b]->path.bbpos < chain->ovls[b - 1]->path.bepos)
						{
							coveredBasesInBread -= chain->ovls[b - 1]->path.bepos - chain->ovls[b]->path.bbpos;
						}
						// check for gap in B
						else
						{
							if (chain->ovls[b]->path.bbpos - chain->ovls[b - 1]->path.bepos > fuzzy)
							{
								int bbpos = chain->ovls[b - 1]->path.bepos;
								int bepos = chain->ovls[b]->path.bbpos;

								if (chain->ovls[b]->flags & OVL_COMP)
								{
									int tmp = bbpos;
									bbpos = DB_READ_LEN(ctx->db, chain->ovls[b]->bread) - bepos;
									bepos = DB_READ_LEN(ctx->db, chain->ovls[b]->bread) - tmp;
								}
								//printf("gap [%d, %d] in %d is low complexity: %d\n", bbpos, bepos, chain->ovls[b]->bread, gapIsLowComplexity(ctx, chain->ovls[b]->bread, bbpos, bepos, 0.8));
								if (!gapIsLowComplexity(ctx, chain->ovls[b]->bread, bbpos, bepos, 0.8))
								{
									properGapLen = 0;
									break;
								}
							}
						}
					}

					int validContainment = 0;
					int validDiff = 0;
					int validMinLen = 0;

					if (ctx->nFuzzBases || ctx->nContPerc)
					{
						if (properGapLen && ((properBegA || properBegB) && (properEndA || properEndB)))
						{
							if (MAX(coveredBasesInAread, coveredBasesInBread) >= (int) (ctx->nContPerc / 100.0 * MIN(DB_READ_LEN(ctx->db, chain->ovls[0]->aread), DB_READ_LEN(ctx->db, chain->ovls[0]->bread))))
							{
								validContainment = 1;
							}
						}
					}
					else // containments are irrelevant
					{
						validContainment = 1;
					}

					if (numdiffs * 100.0 / alignedBasesInAread <= ctx->maxOverallDiff && numdiffs * 100.0 / alignedBasesInBread <= ctx->maxOverallDiff)
						validDiff = 1;

					if (coveredBasesInAread > ctx->minChainLen && coveredBasesInBread > ctx->minChainLen)
						validMinLen = 1;

#ifdef CHAIN_DEBUG
					printf("properBegA %d properBegB %d properEndA %d properEndB %d properGapLen %d validContainment %d validDiff %d validMinLen: %d\n", properBegA, properBegB, properEndA, properEndB, properGapLen, validContainment, validDiff, validMinLen );
#endif
					if ((!properBegA && !properBegB) || (!properEndA && !properEndB) || !properGapLen || !validContainment || !validDiff || !validMinLen)
					{
#ifdef CHAIN_DEBUG
						printf("   *** DISCARD chain ****\n");
#endif
						for (b = 0; b < chain->novl; b++)
						{
							chain->ovls[b]->flags |= OVL_DISCARD;
							ctx->statsFiltInvalidChain++;
						}
					}
					else
					{
						int count = 0;
						for (b = 0; b < chain->novl; b++)
						{
							if (!(chain->ovls[b]->flags & OVL_DISCARD))
								count++;
						}
#ifdef DEBUG_FILTER
						printf(" validLAS %d\n", count);
#endif

						if (ctx->stitchChain)
						{
							int nStitch = stitchChain(ctx, chain);

							if (chain->novl - nStitch > ctx->stitchMaxChainLASs)
							{
								for (b = 0; b < chain->novl; b++)
								{
									chain->ovls[b]->flags |= OVL_DISCARD;
									ctx->statsFiltInvalidChain++;
								}
							}
						}
						count = 0;
#ifdef DEBUG_FILTER


						for (b = 0; b < chain->novl; b++)
						{
							if (!(chain->ovls[b]->flags & OVL_DISCARD))
							count++;
						}
						printf(" validLAS %d\n", count);
#endif
					}
				}
			}
			// reset chain and ovl counter
			for (a = 0; a < ctx->curChains; a++)
				ctx->ovlChains[a].novl = 0;
			ctx->curChains = 0;
		}
		else
		{
		   //printf("DISCARD all overlaps\n");
			if (ovl[j].aread != ovl[j].bread || (ctx->keepIdentity == 0))
			{					// discard all overlaps
				for (i = j; i <= k; i++)
				{
					//printf("DISCARD: %d vs %d a[%d, %d] b[%d, %d] \n",ovl[i].aread,ovl[i].bread,ovl[i].path.abpos,ovl[i].path.aepos,ovl[i].path.bbpos,ovl[i].path.bepos);
					ovl[i].flags |= OVL_DISCARD;
					ctx->statsFiltInvalidChain++;
				}
			}
		}

		k++;
		j = k;
	}

	if (ctx->fileOutDiscardedReads)
	{
		findGaps(ctx, ovl, novl);
	}

	if (ctx->minTipCoverage)
	{
		checkTipCoverage(ctx, ovl, novl);
	}

	if(ctx->fileOutFullyDiscardedAreads)
	{
		if(trim_aend - trim_abeg > ctx->minLenOfFullyDiscardedAreads)
		{
			int count=0;
			for (i = 0; i < novl; i++)
				{
					if (ovl[i].flags & OVL_DISCARD)
					{
						count++;
					}
				}
			if(count == novl)
			{
				fprintf(ctx->fileOutFullyDiscardedAreads, "%d %d %d\n", ovl->aread, trim_abeg, trim_aend);
			}
		}

	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "[-vpiS] [-nkfcmwdyoULGOCMZ <int>] [-BR <file>][-rlt <track>] <db> <overlaps_in> <overlaps_out>\n");

	fprintf(stderr, "options: -v        verbose\n");
	fprintf(stderr, "         -i        keep identity overlaps\n");
	fprintf(stderr, "         -p        purge discarded overlaps\n");
	fprintf(stderr, "         -r <trc>  repeat track name (%s)\n", DEF_ARG_R);
	fprintf(stderr, "         -l <trc>  low complexity track (e.g. tan, dust, tan_dust, default: %s)\n", DEF_ARG_L);
	fprintf(stderr, "         -t <trc>  trim-track (default: %s)\n", DEF_ARG_T);
	fprintf(stderr, "         -k <int>  keep valid overlap chains: 0 ... best, 1 ... all\n");

	fprintf(stderr, "\n 1. Chain overlaps\n");
	fprintf(stderr, "         -n <int>  at least one alignment of a valid chain must have n non-repetitive bases\n");
	fprintf(stderr, "         -m <int>  max merge distance of neighboring repeats (default: %d)\n", DEF_ARG_M);
	fprintf(stderr, "         -w <int>  window size in bases. Merge repeats that are closer then -V bases and have a decent number of low complexity bases in between both repeats\n");
	fprintf(stderr, "                   or at -W bases at the tips of the neighboring repeat. Those can cause a fragmented repeat mask. (default: %d)\n", DEF_ARG_W);
	fprintf(stderr, "         -d <int>  max overall divergence allowed [0,100] (default: %d)\n", DEF_ARG_D);
	fprintf(stderr, "         -y <int>  merge repeats with start/end position of read if repeat interval starts/ends with fewer then -Y\n");
	fprintf(stderr, "         -o <int>  minimum chain length (default: %d)\n", DEF_ARG_O);
	fprintf(stderr, "         -f <int>  allow maximum of -f fuzzy bases \"of structural variations\" between two neighboring overlaps of a chain, (default %d)\n", DEF_ARG_F);
	fprintf(stderr, "         -c <int>  chain alignment must cover at least -p percent of the shorter read. p=[1,100], (default: %d)\n", DEF_ARG_C);

	fprintf(stderr, "\n 2. Stitch chains (optional)\n");
	fprintf(stderr, "         -S        stitch LASchains\n");
	fprintf(stderr, "         -U <int>  maximum unaligned bases for first and last overlap of LAchain (default: 0)\n");
	fprintf(stderr, "         -L <int>  do not merge LAS that cover low complexity regions and have less than -L \"unique\" anchor bases (default: 1000)\n");
	fprintf(stderr, "         -G <int>  maximum merge distance (default: -1)\n");
	fprintf(stderr, "         -M <int>  maximum merge distance in low complexity regions (default: -1)\n");
	fprintf(stderr, "         -O <int>  minimum chain length (default: -1)\n");
	fprintf(stderr, "         -C <int>  max number of LAS in a LASchain (default: 1)\n");

	fprintf(stderr, "\n 3. find Gaps and exclude all those reads (optional)\n");
	fprintf(stderr, "         -B <file> find breaks and gaps and write corresponding read ID into <file>\n");
	fprintf(stderr, "         -N <int>  min number of spanning alignments, to be not a gap (default: 1)\n");

	fprintf(stderr, "\n 4. further Filter parameter (optional)\n");
	fprintf(stderr, "         -E <int>  minimum leaving/entering coverage (default: 0). If coverage is less then -E, than all overlaps are discarded.\n");
	fprintf(stderr, "         -R <file> Write out Aread ids, that are at least Z-bases <int> long and were fully filtered with LAfilterChains.\n");
	fprintf(stderr, "         -Z <int>  minimum A-read length for repeatitive reads (default: %d)\n", DEF_ARG_Z);

}

int main(int argc, char* argv[])
{
	HITS_DB db;
	FilterContext fctx;
	PassContext* pctx;
	FILE* fileOvlIn;
	FILE* fileOvlOut;

	bzero(&fctx, sizeof(FilterContext));

	fctx.db = &db;

// args

	char* pcTrackRepeats = DEF_ARG_R;
	char* pcTrackTrim = DEF_ARG_T;
	char* pcTrackLowCompl = DEF_ARG_L;
	char* pathOutDiscardReads = NULL;
	char* pathOutFullyDiscardAReads = NULL;

	int arg_purge = 0;

	fctx.nMinNonRepeatBases = -1;
	fctx.nVerbose = 0;
	fctx.nkeptChains = 0;
	fctx.nFuzzBases = DEF_ARG_F;
	fctx.nContPerc = DEF_ARG_C;
	fctx.keepIdentity = 0;
	fctx.maxOverallDiff = DEF_ARG_D;
	fctx.mergeRepDist = DEF_ARG_M;
	fctx.repeatWindowLookBack = DEF_ARG_W;
	fctx.rp_mergeTips = 0;
	fctx.minChainLen = DEF_ARG_O;
	fctx.stitchChain = 0;
	fctx.stitchNonLowCompAnchorBases = 1000;
	fctx.stitchMaxGapSize = -1;  // by default stitch all valid chains
	fctx.stitchMaxGapSizeInLowCompl = -1;  // by default stitch all valid chains
	fctx.stitchMaxTipFuzzy = 0;
	fctx.stitchMinChainLen = fctx.minChainLen;
	fctx.stitchMaxChainLASs = 1;
	fctx.minTipCoverage = 0;
	fctx.fileOutDiscardedReads = NULL;
	fctx.gapMinSpanners = 1;
	fctx.fileOutFullyDiscardedAreads = NULL;
	fctx.minLenOfFullyDiscardedAreads = DEF_ARG_Z;
	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "vpn:k:r:f:c:l:t:d:n:m:w:y:io:U:L:G:O:SC:B:E:M:N:R:Z:")) != -1)
	{
		switch (c)
		{
			case 'v':
				fctx.nVerbose = 1;
				break;

			case 'p':
				arg_purge = 1;
				break;

			case 'S':
				fctx.stitchChain = 1;
				break;

			case 'i':
				fctx.keepIdentity = 1;
				break;

			case 'n':
				fctx.nMinNonRepeatBases = atoi(optarg);
				break;

			case 'B':
				pathOutDiscardReads = optarg;
				break;

			case 'R':
				pathOutFullyDiscardAReads = optarg;
				break;

			case 'Z':
				fctx.minLenOfFullyDiscardedAreads = atoi(optarg);
				break;

			case 'E':
				fctx.minTipCoverage = atoi(optarg);
				break;

			case 'C':
				fctx.stitchMaxChainLASs = atoi(optarg);
				break;

			case 'U':
				fctx.stitchMaxTipFuzzy = atoi(optarg);
				break;

			case 'L':
				fctx.stitchNonLowCompAnchorBases = atoi(optarg);
				break;

			case 'G':
				fctx.stitchMaxGapSize = atoi(optarg);
				break;

			case 'M':
				fctx.stitchMaxGapSizeInLowCompl = atoi(optarg);
				break;

			case 'N':
				fctx.gapMinSpanners = atoi(optarg);
				break;

			case 'O':
				fctx.stitchMinChainLen = atoi(optarg);
				break;

			case 'o':
				fctx.minChainLen = atoi(optarg);
				break;

			case 'y':
				fctx.rp_mergeTips = atoi(optarg);
				break;

			case 'f':
				fctx.nFuzzBases = atoi(optarg);
				break;

			case 'c':
				fctx.nContPerc = atoi(optarg);
				break;

			case 'r':
				pcTrackRepeats = optarg;
				break;

			case 'l':
				pcTrackLowCompl = optarg;
				break;

			case 't':
				pcTrackTrim = optarg;
				break;

			case 'd':
				fctx.maxOverallDiff = atoi(optarg);
				break;

			case 'k':
				fctx.nkeptChains = atoi(optarg);
				break;

			case 'm':
				fctx.mergeRepDist = atoi(optarg);
				break;

			case 'w':
				fctx.repeatWindowLookBack = atoi(optarg);
				break;

			default:
				fprintf(stderr, "unknown option %c\n", optopt);
				usage();
				exit(1);
		}
	}

	if (argc - optind != 3)
	{
		usage();
		exit(1);
	}

	char* pcPathReadsIn = argv[optind++];
	char* pcPathOverlapsIn = argv[optind++];
	char* pcPathOverlapsOut = argv[optind++];

	if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL)
	{
		fprintf(stderr, "could not open %s\n", pcPathOverlapsIn);
		exit(1);
	}

	if ((fileOvlOut = fopen(pcPathOverlapsOut, "w")) == NULL)
	{
		fprintf(stderr, "could not open %s\n", pcPathOverlapsOut);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "could not open %s\n", pcPathReadsIn);
		exit(1);
	}

	if (fctx.nMinNonRepeatBases != -1)
	{
		fctx.trackRepeat = track_load(&db, pcTrackRepeats);

		if (!fctx.trackRepeat)
		{
			fprintf(stderr, "could not load track %s\n", pcTrackRepeats);
			exit(1);
		}
	}

	// try to load further non-mandatory tracks
	fctx.trackTrim = track_load(&db, pcTrackTrim);
	if (!fctx.trackTrim)
		fprintf(stderr, "[WARNING] - could not load track %s\n", pcTrackTrim);

	fctx.trackLowCompl = track_load(&db, pcTrackLowCompl);
	if (!fctx.trackLowCompl)
		fprintf(stderr, "[WARNING] - could not load track %s\n", pcTrackLowCompl);

	if (fctx.nContPerc < 0 || fctx.nContPerc > 100)
	{
		fprintf(stderr, "[ERROR] Invalid range for minimum percent of chain alignments %d. Must be in [1,100]\n", fctx.nContPerc);
		exit(1);
	}

	if (fctx.nFuzzBases < 0)
	{
		fprintf(stderr, "[ERROR] -c fuzzy SV bases must be positive! (%d)\n", fctx.nFuzzBases);
		exit(1);
	}

	if (pathOutDiscardReads)
	{
		FILE* fileOut = fopen(pathOutDiscardReads, "w");

		if (fileOut == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathOutDiscardReads);
			exit(1);
		}
		fctx.fileOutDiscardedReads = fileOut;
	}

	if (pathOutFullyDiscardAReads)
	{
		FILE* fileOut = fopen(pathOutFullyDiscardAReads, "w");

		if (fileOut == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathOutFullyDiscardAReads);
			exit(1);
		}
		fctx.fileOutFullyDiscardedAreads = fileOut;
	}


// passes

	pctx = pass_init(fileOvlIn, fileOvlOut);

	pctx->split_b = 0;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &fctx;
	pctx->write_overlaps = 1;
	pctx->purge_discarded = arg_purge;

	filter_pre(pctx, &fctx);

	pass(pctx, filter_handler);

	filter_post(&fctx);

	pass_free(pctx);

// cleanup

	Close_DB(&db);
	if (fctx.fileOutDiscardedReads)
	{
		fclose(fctx.fileOutDiscardedReads);
	}
	if (fctx.fileOutFullyDiscardedAreads)
	{
		fclose(fctx.fileOutFullyDiscardedAreads);
	}
	fclose(fileOvlOut);
	fclose(fileOvlIn);

	return 0;
}
