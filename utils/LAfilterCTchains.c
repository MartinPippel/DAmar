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

#define READ_NONE 0x0
#define READ_DISCARD ( 0x1 << 0 )

#define VERBOSE

#undef DEBUG_REPEAT_EXTENSION
#undef DEBUG_CHAIN

typedef struct
{
	Overlap **ovls;
	int novl;
	int maxOvl;

} Chain;

typedef struct
{
	// stats counters
	int nFilteredRepeat;
	int nDiscardedOvls;
	int nKeptOvls;

	// settings
	int nMinNonRepeatBases;
	int nVerbose;

	int rp_mergeTips;

	HITS_DB* db;
	HITS_TRACK* trackRepeat;

	ovl_header_twidth twidth;

	FILE* fileOutDiscardedOverlaps;

	int nMaxProperChains;
	Chain *ovlChains;
	int curChains;
	int maxChains;
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
	int MAX_OVL_LOOKAHEAD = 10000;
	int STRIDE_OVL_LOOKAHEAD = 2000;

	int trim_bb, trim_be;

	trim_bb = 0;
	trim_be = DB_READ_LEN(ctx->db, ovls->bread);

#ifdef DEBUG_CHAIN
	printf("chain(%d,%d,%d) CHAIN: n%d m%d trim [%d, %d]\n", ovls->aread, ovls->bread, n, ctx->curChains, ctx->maxChains, trim_ab, trim_ae);
#endif
	if (n < 2)
	{
		return;
	}

	int aread, bread;
	int alen, blen;
	int i;

	// reset OVL_TEMP flag
	for (i = 0; i < n; i++)
	{
		if (ovls[i].flags & OVL_TEMP)
			ovls[i].flags &= ~OVL_TEMP;
	}

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

			if (ovl_i->flags & (OVL_CONT))
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

// mark contained overlaps
#ifdef DEBUG_CHAIN
	printf("mark contained overlaps\n");
#endif
	{
		int j;
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = ovls + i;

			if (ovl_i->flags & (OVL_CONT | OVL_TEMP))
				continue;

			for (j = i + 1; j < n; j++)
			{
				Overlap *ovl_j = ovls + j;

				if (ovl_j->flags & (OVL_CONT | OVL_TEMP))
					continue;

				if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos)
						&& contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
				{
					nremain--;
					ovl_j->flags |= (OVL_TEMP);
				}
			}

		}
	}

#ifdef DEBUG_CHAIN
	printf("nremain %d\n", nremain);
#endif

	assert(nremain >= 1);

	while (nremain > 0)
	{
		int longestUniqOvlBases = -1;
		int longestUniqOvlIdx = -1;
		int longestOvlBases = -1;
		int longestOvlIdx = -1;

		// find longest overlap based on number of unique bases
		for (i = 0; i < n; i++)
		{
			Overlap *ovl = ovls + i;

			if (ovl->flags & (OVL_CONT | OVL_TEMP))
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

		assert(longestOvlIdx >= 0);

		if (longestUniqOvlBases < ctx->twidth && longestOvlBases > longestUniqOvlBases)
		{
#ifdef DEBUG_CHAIN
			printf("Number of unique bases to low. Use longest overlap.\n");
#endif
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

							if ((intersect(ab1, ae1, tmpOvl->path.abpos, tmpOvl->path.aepos) < ae1 - tmpOvl->path.abpos)
									&& (intersect(bb1, be1, tmpOvl->path.bbpos, tmpOvl->path.bepos) < be1 - tmpOvl->path.bbpos))
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

							if ((intersect(tmpOvl->path.abpos, tmpOvl->path.aepos, ab1, ae1) < ae1 - tmpOvl->path.abpos)
									&& (intersect(tmpOvl->path.bbpos, tmpOvl->path.bepos, bb1, be1) < be1 - tmpOvl->path.bbpos))
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

		if (chain->novl > 1 && nremain > 0)
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
				if ((ovl->flags & (OVL_TEMP | OVL_CONT)) || ((ovl->flags & OVL_COMP) != (chain->ovls[chainIdx]->flags & OVL_COMP)))
					continue;

				if (ovl->path.abpos < chain->ovls[chainIdx]->path.abpos)
					continue;

				if (ovl->path.abpos > chain->ovls[chainLastIdx]->path.abpos)
					break;

				for (j = chainIdx; j < chainLastIdx; j++)
				{
					if (chain->ovls[j]->path.aepos <= ovl->path.abpos && chain->ovls[j + 1]->path.abpos >= ovl->path.aepos
							&& chain->ovls[j]->path.bepos <= ovl->path.bbpos && chain->ovls[j + 1]->path.bbpos >= ovl->path.bepos)
					{
						Overlap *lastAddedOvl = chain->ovls[chain->novl - 1];

						if (intersect(ovl->path.abpos, ovl->path.aepos, lastAddedOvl->path.abpos, lastAddedOvl->path.aepos)
								|| intersect(ovl->path.bbpos, ovl->path.bepos, lastAddedOvl->path.bbpos, lastAddedOvl->path.bepos))
							break;

						if (chain->novl == chain->maxOvl)
						{
							chain->maxOvl = chain->maxOvl * 1.2 + 5;
							chain->ovls = (Overlap**) realloc(chain->ovls, sizeof(Overlap*) * chain->maxOvl);
						}

						// append left side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
						ovl->flags |= OVL_TEMP;
						chain->ovls[chain->novl] = ovl;
						chain->novl++;
						nremain--;
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

		if (nremain)
		{
			// mark remaining ovls as OVL_TEMP if the overlap with a chain overlap !!
#ifdef DEBUG_CHAIN
			printf("// mark remaining ovls as DISC if they overlap with a chain overlap !!\n");
#endif
			int chainIdx = 0;
			int chainLastIdx = chain->novl - 1;
			int j;
			for (i = 0; i < n; i++)
			{
				Overlap *ovl = ovls + i;
				if ((ovl->flags & (OVL_TEMP | OVL_CONT)) || ((ovl->flags & OVL_COMP) != (chain->ovls[chainIdx]->flags & OVL_COMP)))
					continue;

				for (j = chainIdx; j <= chainLastIdx; j++)
				{
					if (intersect(chain->ovls[j]->path.abpos, chain->ovls[j]->path.aepos, ovl->path.abpos, ovl->path.aepos)
							|| intersect(chain->ovls[j]->path.bbpos, chain->ovls[j]->path.bepos, ovl->path.bbpos, ovl->path.bepos))
					{
						ovl->flags |= OVL_TEMP;
						nremain--;
#ifdef DEBUG_CHAIN
						printf("OVL_TEMP [%d, %d] [%d, %d] nremain %d\n", ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, nremain);
#endif
					}

					if (j + 1 < chain->novl && ovl->path.abpos > chain->ovls[j + 1]->path.abpos)
						chainIdx++;
				}
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
				if ((chain->ovls[0]->flags & OVL_COMP) == (ctx->ovlChains[i].ovls[0]->flags && OVL_COMP))
				{
					for (j = 0; j < chain->novl; j++)
					{
						if ((chain->ovls[j]->path.abpos > ctx->ovlChains[i].ovls[0]->path.abpos
								&& chain->ovls[j]->path.aepos < ctx->ovlChains[i].ovls[ctx->ovlChains[i].novl - 1]->path.aepos)
								|| (chain->ovls[j]->path.bbpos > ctx->ovlChains[i].ovls[0]->path.bbpos
										&& chain->ovls[j]->path.bepos < ctx->ovlChains[i].ovls[ctx->ovlChains[i].novl - 1]->path.bepos))
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
			chain->novl = 0;
		}

#ifdef DEBUG_CHAIN
		printf("curChain: %d, remain unchained OVls: %d\n", ctx->curChains, nremain);
#endif
	}

#ifdef DEBUG_CHAIN
	printf("FINAL CHAINS: %d %7d vs %7d\n", ctx->curChains, ctx->ovlChains[0].ovls[0]->aread, ctx->ovlChains[0].ovls[0]->bread);
	for (i = 0; i < ctx->curChains; i++)
	{
		printf(" CHAIN %d/%d: #novl %d\n", i + 1, ctx->curChains, ctx->ovlChains[0].novl);
		int j;
		for (j = 0; j < ctx->ovlChains[i].novl; j++)
		{
			printf("  OVL %d/%d: a[%7d, %7d] b[%7d, %7d] %s\n", j + 1, ctx->ovlChains[i].novl, ctx->ovlChains[i].ovls[j]->path.abpos,
					ctx->ovlChains[i].ovls[j]->path.aepos, ctx->ovlChains[i].ovls[j]->path.bbpos, ctx->ovlChains[i].ovls[j]->path.bepos,
					(ctx->ovlChains[i].ovls[j]->flags & OVL_COMP) ? "COMP" : "NORM");
		}
	}
#endif

	// sort chains according to alignment lengths
	if (ctx->curChains > 1)
	{
#ifdef DEBUG_CHAIN
		printf("SORT CHAINS (longest first):\n");
#endif
		qsort(ctx->ovlChains, ctx->curChains, sizeof(Chain), cmp_chain_len);

#ifdef DEBUG_CHAIN
		printf("FINAL CHAINS: %d %7d vs %7d\n", ctx->curChains, ctx->ovlChains[0].ovls[0]->aread, ctx->ovlChains[0].ovls[0]->bread);
		for (i = 0; i < ctx->curChains; i++)
		{
			printf(" CHAIN %d/%d: #novl %d\n", i + 1, ctx->curChains, ctx->ovlChains[0].novl);
			int j;
			for (j = 0; j < ctx->ovlChains[i].novl; j++)
			{
				printf("  OVL %d/%d: a[%7d, %7d] b[%7d, %7d] %s\n", j + 1, ctx->ovlChains[i].novl, ctx->ovlChains[i].ovls[j]->path.abpos,
						ctx->ovlChains[i].ovls[j]->path.aepos, ctx->ovlChains[i].ovls[j]->path.bbpos, ctx->ovlChains[i].ovls[j]->path.bepos,
						(ctx->ovlChains[i].ovls[j]->flags & OVL_COMP) ? "COMP" : "NORM");
			}
		}
#endif
	}
}

static int filter(FilterContext* ctx, Overlap* ovl)
{
	int ret = 0;

	int trim_ab, trim_ae, trim_bb, trim_be;
	int trim_alen, trim_blen;

	int ovlALen = DB_READ_LEN(ctx->db, ovl->aread);
	int ovlBLen = DB_READ_LEN(ctx->db, ovl->bread);

	trim_ab = 0;
	trim_ae = ovlALen;

	trim_bb = 0;
	trim_be = ovlBLen;

	trim_alen = ovlALen;
	trim_blen = ovlBLen;

	if (ctx->nMinNonRepeatBases != -1)
	{
		int b, e, rb, re, ovllen, repeat;

		track_anno* repeats_anno = ctx->trackRepeat->anno;
		track_data* repeats_data = ctx->trackRepeat->data;

		int rp_mergeTip_ab = trim_ab;
		int rp_mergeTip_ae = trim_ae;

		if (ctx->rp_mergeTips)
		{
			int cumRep = 0;

			b = repeats_anno[ovl->aread] / sizeof(track_data);
			e = repeats_anno[ovl->aread + 1] / sizeof(track_data);

			while (b < e)
			{
				rb = repeats_data[b];
				re = repeats_data[b + 1];

				// ignore repeats starting behind rp_mergeTip offset
				if (rb > trim_ab + ctx->rp_mergeTips)
					break;

				// ignore repeat in front of trim interval
				if (re < trim_ab)
				{
					b += 2;
					continue;
				}

				if (rb < trim_ab)
				{
					cumRep += re - trim_ab;
				}
				else
				{
					cumRep += re - rb;
				}

				if (re > trim_ab + ctx->rp_mergeTips)
					break;

				b += 2;
			}

			if (cumRep > 1 + ctx->rp_mergeTips / 3)
			{
				rp_mergeTip_ab = trim_ab + ctx->rp_mergeTips;
			}

			// check end of the Aread
			cumRep = 0;

			b = repeats_anno[ovl->aread] / sizeof(track_data);
			e = repeats_anno[ovl->aread + 1] / sizeof(track_data);

			while (b < e)
			{
				rb = repeats_data[b];
				re = repeats_data[b + 1];

				// ignore repeats that end before trim_ae - rp_mergeTip offset
				if (re < trim_ae - ctx->rp_mergeTips)
				{
					b += 2;
					continue;
				}

				// ignore repeat behind end of trim interval
				if (rb > trim_ae)
				{
					break;
				}

				if (re > trim_ae)
				{
					cumRep += trim_ae - rb;
				}
				else
				{
					cumRep += re - rb;
				}

				b += 2;
			}

			if (cumRep > 1 + ctx->rp_mergeTips / 3)
			{
				rp_mergeTip_ae = trim_ae - ctx->rp_mergeTips;
			}
		}

		if (rp_mergeTip_ae < rp_mergeTip_ab)
		{
			if (ctx->nVerbose)
			{
				printf("overlap %d -> %d: drop due to repeat in a\n", ovl->aread, ovl->bread);
			}

			ctx->nFilteredRepeat++;
			ret |= OVL_REPEAT;
		}

		if (!(ret & OVL_REPEAT))
		{
			// Check A-read !!!!

			b = repeats_anno[ovl->aread] / sizeof(track_data);
			e = repeats_anno[ovl->aread + 1] / sizeof(track_data);
			ovllen = MIN(ovl->path.aepos, rp_mergeTip_ae) - MAX(ovl->path.abpos, rp_mergeTip_ab);
			repeat = 0;

			if (ovllen < ctx->nMinNonRepeatBases)
			{
				if (ctx->nVerbose)
				{
					printf("overlap %d -> %d: drop due to repeat in a\n", ovl->aread, ovl->bread);
				}

				ctx->nFilteredRepeat++;
				ret |= OVL_REPEAT;
			}
		}

		if (!(ret & OVL_REPEAT))
		{
			while (b < e)
			{
				rb = repeats_data[b];
				re = repeats_data[b + 1];

				if (re < rp_mergeTip_ab)
				{
					b += 2;
					continue;
				}
				else if (rb > rp_mergeTip_ae)
				{
					break;
				}

				if (rb < rp_mergeTip_ab)
				{
					rb = rp_mergeTip_ab;
				}

				if (re > rp_mergeTip_ae)
				{
					re = rp_mergeTip_ae;
				}

				repeat += intersect(ovl->path.abpos, ovl->path.aepos, rb, re);

				b += 2;
			}

			if (repeat > 0 && ovllen - repeat < ctx->nMinNonRepeatBases)
			{
				if (ctx->nVerbose)
				{
					printf("overlap %d -> %d: drop due to repeat in a\n", ovl->aread, ovl->bread);
				}

				ctx->nFilteredRepeat++;
				ret |= OVL_REPEAT;
			}
		}

		if (!(ret & OVL_REPEAT))
		{
			// check B-Read only if we don't know yet if overlap is discarded by repeat

			b = repeats_anno[ovl->bread] / sizeof(track_data);
			e = repeats_anno[ovl->bread + 1] / sizeof(track_data);

			// roughly map rp_mergeTip positions to B-read
			ovllen = ovl->path.bepos - ovl->path.bbpos;

			if (ovl->path.aepos > rp_mergeTip_ae)
				ovllen -= (ovl->path.aepos - rp_mergeTip_ae);

			if (ovl->path.abpos < rp_mergeTip_ab)
				ovllen -= (rp_mergeTip_ab - ovl->path.abpos);

			if (ovllen < ctx->nMinNonRepeatBases)
			{
				if (ctx->nVerbose)
				{
					printf("overlap %d -> %d: drop due to repeat in a\n", ovl->aread, ovl->bread);
				}

				ctx->nFilteredRepeat++;
				ret |= OVL_REPEAT;
			}
		}

		if (!(ret & OVL_REPEAT))
		{
			int bbpos, bepos;

			if (ovl->flags & OVL_COMP)
			{
				if (ovl->path.aepos > rp_mergeTip_ae)
					bbpos = ovlBLen - (ovl->path.bepos - (ovl->path.aepos - rp_mergeTip_ae));
				else
					bbpos = ovlBLen - ovl->path.bepos;
				if (ovl->path.abpos < rp_mergeTip_ab)
					bepos = ovlBLen - (ovl->path.bbpos - (rp_mergeTip_ab - ovl->path.abpos));
				else
					bepos = ovlBLen - ovl->path.bbpos;
			}
			else
			{
				bbpos = ovl->path.bbpos;
				if (ovl->path.abpos < rp_mergeTip_ab)
					bbpos += rp_mergeTip_ab - ovl->path.abpos;
				bepos = ovl->path.bepos;
				if (ovl->path.aepos > rp_mergeTip_ae)
					bepos -= ovl->path.aepos - rp_mergeTip_ae;
			}

			repeat = 0;

			while (b < e)
			{
				rb = repeats_data[b];
				re = repeats_data[b + 1];

				repeat += intersect(bbpos, bepos, rb, re);

				b += 2;
			}

			if (repeat > 0 && ovllen - repeat < ctx->nMinNonRepeatBases)
			{
				if (ctx->nVerbose)
				{
					printf("overlap %d -> %d: drop due to repeat in b\n", ovl->aread, ovl->bread);
				}

				ctx->nFilteredRepeat++;

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
#endif

	fctx->twidth = pctx->twidth;

	fctx->curChains = 0;
	fctx->maxChains = 5;
	fctx->ovlChains = (Chain*) malloc(sizeof(Chain) * MAX(fctx->maxChains, fctx->nMaxProperChains));
	bzero(fctx->ovlChains, sizeof(Chain) * MAX(fctx->maxChains, fctx->nMaxProperChains));

}

static void filter_post(FilterContext* ctx)
{
#ifdef VERBOSE

	if (ctx->nFilteredRepeat > 0)
	{
		printf("min non-repeat bases of %4d discarded     %10d\n", ctx->nMinNonRepeatBases, ctx->nFilteredRepeat);
	}

	if (ctx->nKeptOvls > 0)
	{
		printf("kept %d overlaps\n", ctx->nKeptOvls);
	}

	if (ctx->nDiscardedOvls > 0)
	{
		printf("discarded %d overlaps\n", ctx->nDiscardedOvls);
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
}

static int filter_handler(void* _ctx, Overlap* ovl, int novl)
{
	FilterContext* ctx = (FilterContext*) _ctx;
	int i, j, k;

	// set filter flags
	for (j = 0; j < novl; j++)
	{
		// get rid of all previous filter
		ovl[j].flags = 0;

		ovl[j].flags |= filter(ctx, ovl + j);
	}

	// get rid of contained overlaps
	if (novl > 1)
	{
		for (i = 0; i < novl; i++)
		{
			Overlap *ovl_i = ovl + i;

			if(ovl_i->flags & (OVL_CONT))
				continue;

			for (j = i + 1; j < novl; j++)
			{
				Overlap *ovl_j = ovl + j;

				if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos)
						|| contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
				{
					ovl_j->flags |= (OVL_CONT);
				}
			}
		}
	}

	int k;
	j = k = 0;

	while (j < novl)
	{

		int nAnchorOvls = (ovl[j].flags & OVL_REPEAT) ? 0 : 1;

		while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
		{
			nAnchorOvls += (ovl[k + 1].flags & (OVL_REPEAT | OVL_CONT)) ? 0 : 1;
			k++;
		}

		// ignore all self alignments if those are present
		if (ovl[j].aread == ovl[j].bread)
			nAnchorOvls = 0;

		printf("Aread: %8d, Bread: %8d, novl: %10d, anchorOvl: %5d\n", ovl[j].aread, ovl[j].bread, k - j + 1, nAnchorOvls);
		if (nAnchorOvls)
		{
			chain(ctx, ovl + j, k - j + 1);

			// evaluate user input: nMaxProperChains and get rid (DISCARD all overlaps) of chains that are not wanted
			// if nMaxProperChains == 0, then only THE best chain is kept, all others are discarded
			// if nMaxProperChains == 1, then only THE best chain in forward as well as in reverse complement orientation are kept
			// otherwise keep all chains

			int a, b;
			// get rid of OVL_TEMP marked ovls
			for (a = j; a < k - j + 1; a++)
				;
			if (ovl[a].flags & OVL_TEMP)
				ovl[a].flags |= OVL_DISCARD;

			if (ctx->nMaxProperChains == 0)
			{
				for (a = 1; a < ctx->curChains; a++)
					for (b = 0; b < ctx->ovlChains[a].novl; b++)
						ctx->ovlChains[a].ovls[b]->flags |= OVL_DISCARD;
			}
			else if (ctx->nMaxProperChains == 1)
			{
				int keptOtherOri = 0;
				for (a = 1; a < ctx->curChains; a++)
				{
					if (!keptOtherOri && ((ctx->ovlChains[0].ovls[0]->flags & OVL_COMP) != (ctx->ovlChains[a].ovls[0]->flags & OVL_COMP)))
					{
						keptOtherOri = 1;
						continue;
					}

					for (b = 0; b < ctx->ovlChains[a].novl; b++)
						ctx->ovlChains[a].ovls[b]->flags |= OVL_DISCARD;
				}
			}

			// reset chain and ovl counter
			for (a = 0; a < ctx->curChains; a++)
				ctx->ovlChains[a].novl = 0;
			ctx->curChains = 0;
		}
		else
		{
			// discard all overlaps
			for (i = j; i <= k; i++)
			{
				ovl[i].flags |= OVL_DISCARD;
			}
		}

		j = k + 1;
	}

	// filter by read flags
	int aread = ovl->aread;
	HITS_READ* reads = ctx->db->reads;

	if (reads[aread].flags & READ_DISCARD)
	{
		for (j = 0; j < novl; j++)
		{
			ovl[j].flags |= OVL_DISCARD;
		}
	}
	else
	{
		for (j = 0; j < novl; j++)
		{
			int bread = ovl[j].bread;

			if (reads[bread].flags & READ_DISCARD)
			{
				ovl[j].flags |= OVL_DISCARD;
			}
		}
	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "[-vp] [-nCy <int>] [-r <track>] [-a <file>] <db> <overlaps_in> <overlaps_out>\n");

	fprintf(stderr, "options: -v      	verbose\n");
	fprintf(stderr, "         -n <int>	at least one alignment of a valid chain must have n non-repetitive bases\n");
	fprintf(stderr, "         -p      	purge discarded overlaps\n");
	fprintf(stderr, "         -r <trc>	repeat track name (%s)\n", DEF_ARG_R);
	fprintf(stderr, "         -y <int>	extend repeat mask at contig tips by -y bases)\n");
	fprintf(stderr, "         -a <fle> 	write discarded overlaps that may not symmetrically removed to file\n");
	fprintf(stderr, "         -C <int>  keep valid overlap chains: 0 ... best, 1 ... all\n");
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

	char* pathOutDiscardedOvls = NULL;
	char* pcTrackRepeats = DEF_ARG_R;
	int arg_purge = 0;

	fctx.nMinNonRepeatBases = -1;
	fctx.nVerbose = 0;
	fctx.rp_mergeTips = 0;
	fctx.nMaxProperChains = 1;

	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "vpn:C:r:a:y:")) != -1)
	{
		switch (c)
		{
		case 'a':
			pathOutDiscardedOvls = optarg;
			break;

		case 'y':
			fctx.rp_mergeTips = atoi(optarg);
			break;

		case 'v':
			fctx.nVerbose = 1;
			break;

		case 'p':
			arg_purge = 1;
			break;

		case 'n':
			fctx.nMinNonRepeatBases = atoi(optarg);
			break;

		case 'r':
			pcTrackRepeats = optarg;
			break;

		case 'C':
			fctx.nMaxProperChains = atoi(optarg);
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

	int i;
	for (i = 0; i < DB_NREADS(&db); i++)
	{
		db.reads[i].flags = READ_NONE;
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

	if (pathOutDiscardedOvls)
	{
		FILE* fileOut = fopen(pathOutDiscardedOvls, "w");

		if (fileOut == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathOutDiscardedOvls);
			exit(1);
		}
		fctx.fileOutDiscardedOverlaps = fileOut;
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

	if (fctx.fileOutDiscardedOverlaps)
	{
		fclose(fctx.fileOutDiscardedOverlaps);
	}

	fclose(fileOvlOut);
	fclose(fileOvlIn);

	return 0;
}
