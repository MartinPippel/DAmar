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

#define DEF_ARG_Q TRACK_Q
#define DEF_ARG_T TRACK_TRIM

#define VERBOSE

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
	int nFiltReads;
	int nFiltAlns;
	int nFilteredMinLenReads;
	int nFilteredMaxLenReads;
	int nFilteredMinLenAln;

	// settings
	int nMinChainAlnLength;
	int nMinReadLength, nMaxReadLength;
	int nVerbose;

	HITS_DB* db;
	HITS_TRACK* trackQ;
	HITS_TRACK* trackTrim;

	ovl_header_twidth twidth;

	Chain *ovlChains;
	int curChains;
	int maxChains;
} FilterContext;

extern char* optarg;
extern int optind, opterr, optopt;

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

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, ovls->bread, &trim_bb, &trim_be);
	}
	else
	{
		trim_bb = 0;
		trim_be = DB_READ_LEN(ctx->db, ovls->bread);
	}
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

#ifdef DEBUG_CHAIN
			printf("%d - %d [%d, %d] [%d, %d], aR %d/%d, bR %d/%d\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, aLen, aRep,
					bLen, bRep);
#endif

			int tmpBases = MAX(aLen, bLen);
			if (tmpBases > longestOvlBases)
			{
				longestOvlBases = tmpBases;
				longestOvlIdx = i;
			}
		}

		assert(longestOvlIdx >= 0);

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

		chain->ovls[0] = ovls + longestOvlIdx;
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

		if (nremain && longestOvlIdx + 1 < n)
		{
			ab1 = ovls[longestOvlIdx].path.abpos;
			ae1 = ovls[longestOvlIdx].path.aepos;

			bb1 = ovls[longestOvlIdx].path.bbpos;
			be1 = ovls[longestOvlIdx].path.bepos;

#ifdef DEBUG_CHAIN
			printf("extend longest overlap in right direction\n");
#endif
// 1st right
			int cont = 1;
			int curBestBases = -1;
			int curBestOffset = 1;
			int curBestIntersection = MAX(alen, blen);

			while (cont)
			{
				int stepSize;
				for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
				{
#ifdef DEBUG_CHAIN
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif
					for (i = longestOvlIdx + curBestOffset; i < n; i++)
					{
						Overlap * ovl = ovls + i;

						ab2 = ovl->path.abpos;
						ae2 = ovl->path.aepos;

						bb2 = ovl->path.bbpos;
						be2 = ovl->path.bepos;

						if ((ovl->flags & OVL_COMP) != (ovls[longestOvlIdx].flags & OVL_COMP))
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

						int curBasesInAIvl = (ae2 - ab2);
						int curBasesInBIvl = (be2 - bb2);

						if (curBestIntersection > MAX(intersect(ab1, ae1, ab2, ae2), intersect(bb1, be1, bb2, be2)) && curBestBases < MIN(ae2 - ab2, be2 - bb2))
						{
							curBestBases = MIN(ae2 - ab2, be2 - bb2);
							curBestOffset = i - longestOvlIdx;
							curBestIntersection = MAX(intersect(ab1, ae1, ab2, ae2), intersect(bb1, be1, bb2, be2));
						}

						else if (curBestBases < MIN(curBasesInAIvl, curBasesInBIvl))
						{
#ifdef DEBUG_CHAIN
							printf("found right current best overlap %d %d %d: [%d, %d] [%d, %d] right side\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							curBestOffset = i - longestOvlIdx;
							curBestBases = MIN(curBasesInAIvl, curBasesInBIvl);
						}
						else if (curBestBases == -1 && stepSize + STRIDE_OVL_LOOKAHEAD > MAX_OVL_LOOKAHEAD) // for repetitive genomes
						{
							Overlap *tmpOvl = ovls + (longestOvlIdx + curBestOffset);

							if ((intersect(ab1, ae1, tmpOvl->path.abpos, tmpOvl->path.aepos) < ae1 - tmpOvl->path.abpos)
									&& (intersect(bb1, be1, tmpOvl->path.bbpos, tmpOvl->path.bepos) < be1 - tmpOvl->path.bbpos))
							{
								curBestBases = 1;
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
				if (curBestBases < 0) // i.e. there was no good overlap at right side
				{
#ifdef DEBUG_CHAIN
					printf("could not extend ovlgroup on right side with proper overlap (with stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
					break;
				}

				/// todo further sanity check necessary ???
				ab2 = ovls[longestOvlIdx + curBestOffset].path.abpos;
				ae2 = ovls[longestOvlIdx + curBestOffset].path.aepos;

				bb2 = ovls[longestOvlIdx + curBestOffset].path.bbpos;
				be2 = ovls[longestOvlIdx + curBestOffset].path.bepos;
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
				chain->ovls[chain->novl] = ovls + (longestOvlIdx + curBestOffset);
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

				curBestOffset++;

				curBestBases = -1;

				curBestIntersection = MAX(alen, blen);

				if (longestOvlIdx + curBestOffset >= n)
				{
					cont = 0;
				}
			}
		}

		if (nremain && longestOvlIdx > 0)
		{
			ab1 = ovls[longestOvlIdx].path.abpos;
			ae1 = ovls[longestOvlIdx].path.aepos;

			bb1 = ovls[longestOvlIdx].path.bbpos;
			be1 = ovls[longestOvlIdx].path.bepos;

#ifdef DEBUG_CHAIN
			printf("extend longest overlap in left direction\n");
#endif
// 2nd left side
			int cont = 1;
			int curBestBases = -1;
			int curBestOffset = 1;
			int curBestIntersection = MAX(alen, blen);

			while (cont)
			{
				int stepSize;
				for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
				{
#ifdef DEBUG_CHAIN
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif

					// try to find next best overlap with lookahead of stepSize bases
					for (i = longestOvlIdx - curBestOffset; i >= 0; --i)
					{
						Overlap * ovl = ovls + i;

#ifdef DEBUG_CHAIN
						printf("LEFT: Check ovl: a[%d, %d] b[%d,%d]\n", ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos);
#endif

						ab2 = ovl->path.abpos;
						ae2 = ovl->path.aepos;

						bb2 = ovl->path.bbpos;
						be2 = ovl->path.bepos;

						if ((ovl->flags & OVL_COMP) != (ovls[longestOvlIdx].flags & OVL_COMP))
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

						int curBasesInAIvl = (ae2 - ab2);
						int curBasesInBIvl = (be2 - bb2);

						if (curBestIntersection > MAX(intersect(ab2, ae2, ab1, ae1), intersect(bb2, be2, bb1, be1)) && curBestBases < MIN(ae2 - ab2, be2 - bb2))
						{
							curBestBases = MIN(ae2 - ab2, be2 - bb2);
							curBestOffset = longestOvlIdx - i;
							curBestIntersection = MAX(intersect(ab2, ae2, ab1, ae1), intersect(bb2, be2, bb1, be1));
						}

						else if (curBestBases < MIN(curBasesInAIvl, curBasesInBIvl))
						{
#ifdef DEBUG_CHAIN
							printf("found left current best overlap %d %d %d: [ab2 %d, ae2 %d] [bb2 %d, be2 %d] left side\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
#endif
							curBestOffset = longestOvlIdx - i;
							curBestBases = curBasesInAIvl + curBasesInBIvl;
						}
						else if (curBestBases == -1 && stepSize + STRIDE_OVL_LOOKAHEAD > MAX_OVL_LOOKAHEAD) // for repetitive genomes
						{
							Overlap *tmpOvl = ovls + (longestOvlIdx - curBestOffset);

							if ((intersect(tmpOvl->path.abpos, tmpOvl->path.aepos, ab1, ae1) < ae1 - tmpOvl->path.abpos)
									&& (intersect(tmpOvl->path.bbpos, tmpOvl->path.bepos, bb1, be1) < be1 - tmpOvl->path.bbpos))
							{
								curBestBases = 1;
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
				if (curBestBases < 0) // i.e. there was no good overlap at left side
				{
#ifdef DEBUG_CHAIN
					printf("could not extend ovlgroup on left side with proper overlap (stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
					break;
				}

				/// todo further sanity check necessary ???
				ab2 = ovls[longestOvlIdx - curBestOffset].path.abpos;
				ae2 = ovls[longestOvlIdx - curBestOffset].path.aepos;

				bb2 = ovls[longestOvlIdx - curBestOffset].path.bbpos;
				be2 = ovls[longestOvlIdx - curBestOffset].path.bepos;

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
				chain->ovls[chain->novl] = ovls + (longestOvlIdx - curBestOffset);
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

				curBestOffset++;

				curBestBases = -1;

				curBestIntersection = MAX(alen, blen);

				if (longestOvlIdx - curBestOffset < 0)
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

static void filter_pre(PassContext* pctx, FilterContext* fctx)
{
#ifdef VERBOSE
	printf( ANSI_COLOR_GREEN "PASS filtering\n" ANSI_COLOR_RESET);
#endif

	fctx->twidth = pctx->twidth;

	fctx->curChains = 0;
	fctx->maxChains = 5;
	fctx->ovlChains = (Chain*) malloc(sizeof(Chain) * MAX(fctx->maxChains, 5));
	bzero(fctx->ovlChains, sizeof(Chain) * MAX(fctx->maxChains, 5));

}

static void filter_post(FilterContext* ctx)
{
	if (ctx->nVerbose)
	{
		if (ctx->nFiltAlns > 0)
		{
			printf("total discarded alignments %10d\n", ctx->nFiltAlns);
		}
	}

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
	int j;

	{
		int k;
		j = k = 0;

		int trimABeg, trimAEnd;

		trimABeg = 0;
		trimAEnd = DB_READ_LEN(ctx->db, ovl->aread);

		if (ctx->trackTrim)
			get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trimABeg, &trimAEnd);

		while (j < novl)
		{
			int bReadLen = DB_READ_LEN(ctx->db, ovl[j].bread);
			int skipAlns = 0;

			if (ctx->nMinReadLength != -1 && (bReadLen < ctx->nMinReadLength || bReadLen > ctx->nMaxReadLength))
			{
				skipAlns = 1;
				ovl[j].flags |= OVL_DISCARD | OVL_RLEN;
			}

			while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
			{
				if(skipAlns)
					ovl[k+1].flags |= OVL_DISCARD | OVL_RLEN;
				k++;
			}

			if(!skipAlns)
			{
				chain(ctx, ovl + j, k - j + 1);
				int a,b;

				printf("Chains for %d vs %d \n", ovl[j].aread, ovl[j].bread);
				for (a = 0; a < ctx->curChains; a++)
				{
					printf(" chains %d of %d\n", a+1, ctx->curChains);
					for (b = 0; b < ctx->ovlChains[a].novl; b++)
					{
						printf("  a[%d, %d] b[%d, %d] %c l[%d, %d]\n", ctx->ovlChains[a].ovls[b]->path.abpos, ctx->ovlChains[a].ovls[b]->path.aepos,
								ctx->ovlChains[a].ovls[b]->path.bbpos, ctx->ovlChains[a].ovls[b]->path.bepos, ctx->ovlChains[a].ovls[b]->flags & OVL_COMP ? 'C' : 'N',
										DB_READ_LEN(ctx->db, ovl[j].aread), DB_READ_LEN(ctx->db, ovl[j].bread));
					}

					// keep alignment chains if:
					// and only if:
					// 		1. single front to end alignment (TODO Q1: Allow unaligned front ends < 1K?, Q2: run forcealign + do not allow unaligned tips)
					// 		2. 2 chains with same orientation and circular mapping
					// otherwise
					// 	delete all alignments


				}
				// reset chain and ovl counter
				for (a = 0; a < ctx->curChains; a++)
					ctx->ovlChains[a].novl = 0;
				ctx->curChains = 0;
			}

			j = k + 1;
		}
	}

	// set filter flags
	for (j = 0; j < novl; j++)
	{
		ovl[j].flags |= filter(ctx, ovl + j);
	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "[-vp] [-lLdDo <int>] [-q <track>] <db> <overlaps_in> <overlaps_out>\n");

	fprintf(stderr, "options: -v ... verbose\n");
	fprintf(stderr, "         -d ... max divergence allowed [0,100]\n");
	fprintf(stderr, "         -o ... min overlap length\n");
	fprintf(stderr, "         -l ... min read length\n");
	fprintf(stderr, "         -L ... max read length\n");
	fprintf(stderr, "         -p ... purge discarded overlaps\n");
	fprintf(stderr, "         -q ... q track name (%s)\n", DEF_ARG_Q);
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

	char* pcTrackQ = DEF_ARG_Q;
	char* pcTrackT = DEF_ARG_T;
	int arg_purge = 0;

	fctx.nMinChainAlnLength = -1;
	fctx.nMinReadLength = -1;
	fctx.nMaxReadLength = -1;
	fctx.nVerbose = 0;

	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "vpl:L:d:o:q:t:")) != -1)
	{
		switch (c)
		{
		case 'v':
			fctx.nVerbose += 1;
			break;

		case 'p':
			arg_purge = 1;
			break;

		case 'o':
			fctx.nMinChainAlnLength = atoi(optarg);
			break;

		case 'l':
			fctx.nMinReadLength = atoi(optarg);
			break;

		case 'L':
			fctx.nMaxReadLength = atoi(optarg);
			break;

		case 'q':
			pcTrackQ = optarg;
			break;

		case 't':
			pcTrackT = optarg;
			break;

		default:
			fprintf(stderr, "unknown option %c\n", c);
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

	fctx.trackTrim = track_load(&db, pcTrackT);

	if (!fctx.trackTrim)
	{
		fprintf(stderr, "could not load track %s\n", pcTrackT);
	}

	fctx.trackQ = track_load(&db, pcTrackQ);

	if (!fctx.trackQ)
	{
		fprintf(stderr, "could not load track %s\n", pcTrackQ);
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

	fclose(fileOvlOut);
	fclose(fileOvlIn);

	return 0;
}
