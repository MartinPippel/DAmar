/*******************************************************************************************
 *
 *  Author  :  MARVEL Team
 *
 *  Date    :  February 2016
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <assert.h>
#include <unistd.h>

#include "lib/tracks.h"
#include "lib/pass.h"
#include "lib/oflags.h"
#include "lib/colors.h"
#include "lib/utils.h"
#include "lib/trim.h"

#include "db/DB.h"
#include "dalign/align.h"

// argument defaults

#define DEF_ARG_P 0
#define DEF_ARG_S 0

// thresholds

#define BIN_SIZE    100
#define MIN_LR      450

// switches

#define VERBOSE

#define DEBUG_GAPS
#define DEBUG_CGAPS
#undef DEBUG_CHAIN
#undef DEBUG_CHIMER

// constants

typedef struct
{
	Overlap **ovls;
	int novl;
	int maxOvl;
} Chain;

typedef struct
{
	HITS_DB* db;
	int twidth;
	HITS_TRACK* trackExclude;
	HITS_TRACK* trackRepeat;

	unsigned int stats_contained;
	unsigned int stats_breaks;
	unsigned int stats_breaks_novl;

	unsigned int stats_chimers;
	unsigned int stats_chimers_novl;

	int stitch;

	uint64* rm_bins;
	int rm_maxbins;

	int* rc_left;
	int rc_maxright;

	int* rc_right;
	int rc_maxleft;

	// for trimming
	HITS_TRACK* trackTrim;
	int useRLoader;
	TRIM* trim;
	Read_Loader *rl;

	FILE* chimer;
	Chain *ovlChains;
	int curChains;
	int maxChains;

} GapContext;

// for getopt()

extern char* optarg;
extern int optind, opterr, optopt;

static int loader_handler(void* _ctx, Overlap* ovl, int novl)
{
	GapContext* ctx = (GapContext*) _ctx;
	Read_Loader* rl = ctx->rl;

	int i;
	for (i = 0; i < novl; i++)
	{
		int b = ovl[i].bread;

		int trim_b_left, trim_b_right;
		get_trim(ctx->db, ctx->trackTrim, b, &trim_b_left, &trim_b_right);

		if (ovl[i].flags & OVL_COMP)
		{
			int tmp = trim_b_left;
			int blen = DB_READ_LEN(ctx->db, ovl[i].bread);
			trim_b_left = blen - trim_b_right;
			trim_b_right = blen - tmp;
		}

		if (trim_b_left >= trim_b_right)
		{
			continue;
		}

		int bbt = MAX(trim_b_left, ovl[i].path.bbpos);
		int bet = MIN(trim_b_right, ovl[i].path.bepos);

		if (bbt >= bet)
		{
			continue;
		}

		if (bbt == ovl[i].path.bbpos && bet == ovl[i].path.bepos)
		{
			continue;
		}

		bbt = MAX(trim_b_left, ovl[i].path.bbpos);
		bet = MIN(trim_b_right, ovl[i].path.bepos);

		if (bbt < bet && (bbt != ovl[i].path.bbpos || bet != ovl[i].path.bepos))
		{
			rl_add(rl, ovl[i].aread);
			rl_add(rl, ovl[i].bread);

			continue;
		}

		int bepos = ovl[i].path.bepos;

		if (bepos > bet)
		{
			rl_add(rl, ovl[i].aread);
			rl_add(rl, ovl[i].bread);
		}
	}

	return 1;
}

static int drop_break(Overlap* pOvls, int nOvls, int p1, int p2)
{
	int max = 0;
	int left = 0;
	int i;

	for (i = 0; i < nOvls; i++)
	{
		if (pOvls[i].flags & OVL_DISCARD)
		{
			continue;
		}

		int len = pOvls[i].path.aepos - pOvls[i].path.abpos;
		if (len > max)
		{
			max = len;

			if (pOvls[i].path.abpos < p1)
			{
				left = 1;
			}
			else
			{
				left = 0;
			}
		}
	}

	if (max == 0)
	{
		return 0;
	}

	int dropped = 0;

	if (!left)
	{
#ifdef DEBUG_GAPS
		printf("DROP L @ %5d..%5d\n", p1, p2);
#endif

		for (i = 0; i < nOvls; i++)
		{
			if (pOvls[i].path.abpos < p1)
			{
				if (!(pOvls[i].flags & OVL_DISCARD))
				{
					dropped++;
				}

				pOvls[i].flags |= OVL_DISCARD | OVL_GAP;
			}
		}
	}
	else
	{
#ifdef DEBUG_GAPS
		printf("DROP R @ %5d..%5d\n", p1, p2);
#endif

		for (i = 0; i < nOvls; i++)
		{
			if (pOvls[i].flags & OVL_DISCARD)
				continue;

			if (pOvls[i].path.aepos > p2)
			{
				if (!(pOvls[i].flags & OVL_DISCARD))
				{
					dropped++;
				}

				pOvls[i].flags |= OVL_DISCARD | OVL_GAP;
			}
		}
	}

	return dropped;
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

static int contained(int ab, int ae, int bb, int be)
{
	if (ab >= bb && ae <= be)
	{
		return 1;
	}

	if (bb >= ab && be <= ae)
	{
		return 2;
	}

	return 0;
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

static int stitchable(Overlap* pOvls, int n, int fuzz, int beg, int end)
{
	if (n < 2)
	{
		return 0;
	}

	int t = 0;
	int k, b;
	int ab2, ae1;
	int bb2, be1;

	const int ignore_mask = OVL_TEMP | OVL_CONT | OVL_TRIM | OVL_STITCH;

	int i;
	for (i = 0; i < n; i++)
	{
		if (pOvls[i].flags & ignore_mask)
		{
			continue;
		}

		b = pOvls[i].bread;

		ae1 = pOvls[i].path.aepos;
		be1 = pOvls[i].path.bepos;

		for (k = i + 1; k < n && pOvls[k].bread == b; k++)
		{
			if ((pOvls[k].flags & ignore_mask) || ((pOvls[i].flags & OVL_COMP) != (pOvls[k].flags & OVL_COMP)))
			{
				continue;
			}

			ab2 = pOvls[k].path.abpos;
			bb2 = pOvls[k].path.bbpos;

			int deltaa = abs(ae1 - ab2);
			int deltab = abs(be1 - bb2);

			if (deltaa < fuzz && deltab < fuzz)
			{
				if (pOvls[i].path.abpos < beg - MIN_LR && pOvls[k].path.aepos > end + MIN_LR)
				{
					t++;

					// printf("stitch using b %d\n", b);
				}
			}

		}
	}

	return t;
}

static int ovl_intersect(Overlap* a, Overlap* b)
{
	return intersect(a->path.abpos, a->path.aepos, b->path.abpos, b->path.aepos);
}

static int getRepeatCount(GapContext* ctx, int readID, int beg, int end)
{
	if (ctx->trackRepeat == NULL)
		return 0;

	track_anno* rep_anno = ctx->trackRepeat->anno;
	track_data* rep_data = ctx->trackRepeat->data;

	track_anno rb, re;

	int repBases = 0;
	int rBeg, rEnd;

	if (readID < 0 || readID >= DB_NREADS(ctx->db))
	{
		fprintf(stderr, "[ERROR] - getRepeatCount readID: %d out of bounds [0, %d]\n", readID, DB_NREADS(ctx->db) - 1);
		fflush(stderr);
		exit(1);
	}

	// repeat bases in a-read
	rb = rep_anno[readID] / sizeof(track_data);
	re = rep_anno[readID + 1] / sizeof(track_data);

	while (rb < re)
	{
		rBeg = rep_data[rb];
		rEnd = rep_data[rb + 1];

		repBases += intersect(beg, end, rBeg, rEnd);

		rb += 2;
	}

	return repBases;
}

static int getRepeatBases(GapContext *ctx, Overlap *ovl, int read)
{
	if (ctx->trackRepeat == NULL)
	{
		return 0;
	}

	assert(ovl->aread == read || ovl->bread == read);

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

static int spanningChain(GapContext*ctx, Overlap *ovls, int n, int chimerBeg, int chimerEnd, int trim_ab, int trim_ae)
{
	/// TODO hard coded
	int MIN_OVL_LOOKAHEAD = 2000;
	int MAX_OVL_LOOKAHEAD = 10000;
	int STRIDE_OVL_LOOKAHEAD = 2000;
	int MIN_ANCHOR = 800;
	int MIN_CHAIN_LEN = 3000;

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
	printf("chain(%d,%d,%d) CHAIN: n%d m%d chim [%d, %d] trim [%d, %d]\n", ovls->aread, ovls->bread, n, ctx->curChains, ctx->maxChains, chimerBeg, chimerEnd, trim_ab, trim_ae);
#endif
	if (n < 2)
	{
		if (
				(ovls->path.abpos < chimerBeg && ovls->path.aepos > chimerEnd && (ovls->path.abpos == trim_ab || ovls->path.bbpos <= trim_bb) && (ovls->path.aepos == trim_ae || ovls->path.bepos >= trim_be))
				|| (ovls->path.abpos == trim_ab && ovls->path.aepos == trim_ae) // fully spanning trim interval
			 )
		{
#ifdef DEBUG_CHAIN
			printf("found spanning overlap %d vs %d\n", ovls->aread, ovls->bread);
#endif
			return 1;
		}
		return 0;
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
		return 0;

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

	// check if chain spans putative chimer
	Chain* bestChain = ctx->ovlChains;

	if (bestChain->ovls[0]->path.abpos + MIN(MIN_ANCHOR, chimerBeg - trim_ab) > chimerBeg)
		return 0;

	if (bestChain->ovls[bestChain->novl - 1]->path.aepos - MIN(MIN_ANCHOR, trim_ae - chimerEnd) < chimerEnd)
		return 0;

	int chainLen = 0;
	for (i = 0; i < bestChain->novl; i++)
		chainLen += bestChain->ovls[i]->path.aepos - bestChain->ovls[i]->path.abpos;

	if (chainLen < MIN_CHAIN_LEN)
		return 0;

	// chain must be proper

	// todo
	// 1. allow not fully proper chain? i.e. that start/end near the tips of the trim intervals????
	// 2. set a minimum overlaps length to avoid 'sparse chains' (only induced be repeats)
	// 3. set a minimum of non-repeat bases for a chain
	if (((bestChain->ovls[0]->path.abpos <= trim_ab) || (bestChain->ovls[0]->path.bbpos <= trim_bb))
			&& ((bestChain->ovls[bestChain->novl - 1]->path.aepos >= trim_ae) || (bestChain->ovls[bestChain->novl - 1]->path.bepos >= trim_be)))
	{
#ifdef DEBUG_CHAIN
		printf("found spanning chain %d vs %d\n", ovls->aread, ovls->bread);
#endif
		return 1;
	}

	return 0;
}

static int oChainIntervalIntersection(GapContext *fctx, Overlap* ovls, int novl, int b, int e, int trim_ab, int trim_ae)
{
	int spanningChains = 0;

	int j, k;
	j = k = 0;

	while (j < novl)
	{
		while (k < novl - 1 && ovls[j].bread == ovls[k + 1].bread)
		{
			k++;
		}

		spanningChains += spanningChain(fctx, ovls + j, k - j + 1, b, e, trim_ab, trim_ae);
		{
			// reset chain and ovl counter
			int i;
			for (i = 0; i < fctx->curChains; i++)
				fctx->ovlChains[i].novl = 0;
			fctx->curChains = 0;
		}

		if(spanningChains)
			return spanningChains;

		j = k + 1;
	}
	return spanningChains;
}

static int handle_chimers(GapContext* ctx, Overlap* ovls, int novl)
{
	int *alnBeg = (int *) malloc((DB_READ_MAXLEN(ctx->db) / ctx->twidth + 1) * sizeof(int));
	memset(alnBeg, 0, (DB_READ_MAXLEN(ctx->db) / ctx->twidth + 1) * sizeof(int));

	int *alnEnd = (int *) malloc((DB_READ_MAXLEN(ctx->db) / ctx->twidth + 1) * sizeof(int));
	memset(alnEnd, 0, (DB_READ_MAXLEN(ctx->db) / ctx->twidth + 1) * sizeof(int));

	int trim_ab, trim_ae;
	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, ovls->aread, &trim_ab, &trim_ae);
	}
	else
	{
		trim_ab = 0;
		trim_ae = DB_READ_LEN(ctx->db, ovls->aread);
	}

	if (trim_ab >= trim_ae)
	{
		return 1;
	}

#ifdef DEBUG_CHIMER
	printf("filter_chimers aread %d #ovls %d trim [%d, %d]\n", ovls->aread, novl, trim_ab, trim_ae);
#endif
	int i;

	for (i = 1; i < novl - 1; i++)
	{
		Overlap *opre = ovls + i - 1;
		Overlap *ocur = ovls + i;
		Overlap *opas = ovls + i + 1;

		int prematEnd = (ocur->path.bepos + 1000 < DB_READ_LEN(ctx->db, ocur->bread)) && (ocur->path.aepos + 1000 < trim_ae) && (ocur->path.aepos - 100 > trim_ab);
		int prematBeg = (ocur->path.bbpos > 1000) && (ocur->path.abpos > trim_ab + 1000) && (ocur->path.abpos + 100 < trim_ae);

//		// check weird case first
		if (prematBeg && prematEnd && (opre->bread != ocur->bread) && (ocur->bread != opas->bread)) // most probably a repeat
			continue;

		if ((ocur->bread != opas->bread) && prematEnd)
		{
			alnEnd[(ocur->path.aepos / ctx->twidth)] += 1;
#ifdef DEBUG_CHIMER
			printf("found end %d vs %d a[%d, %d] b[%d, %d] ends[%d]=%d\n", ocur->aread, ocur->bread, ocur->path.abpos, ocur->path.aepos, ocur->path.bbpos,
					ocur->path.bepos, (ocur->path.aepos / ctx->twidth), alnEnd[(ocur->path.aepos / ctx->twidth)]);
#endif
		}

		else if ((opre->bread != ocur->bread) && prematBeg)
		{
			alnBeg[(ocur->path.abpos / ctx->twidth)] += 1;
#ifdef DEBUG_CHIMER
			printf("found beg %d vs %d a[%d, %d] b[%d, %d] beg[%d]=%d\n", ocur->aread, ocur->bread, ocur->path.abpos, ocur->path.aepos, ocur->path.bbpos,
					ocur->path.bepos, (ocur->path.abpos / ctx->twidth), alnBeg[(ocur->path.abpos / ctx->twidth)]);
#endif
		}
	}

	int chimBeg = -1;
	int chimEnd = -1;

	int numChim = 0;

	int j;
	i = j = 0;

	// todo border coverage hard coded - find a dynamic approach
	int BORDER_COV = 5;
	// todo hard coded -
	int MAX_CHIMER_LEN = 12000;
	for (i = 0; i <= trim_ae / ctx->twidth && numChim == 0; i++)
	{
		if (alnEnd[i] >= BORDER_COV)
		{
			chimBeg = i;
		}
		else if (alnEnd[i + 1] >= BORDER_COV)
		{
			chimBeg = ++i;
		}
		else if (alnEnd[i] + alnEnd[i + 1] >= BORDER_COV)
		{
			chimBeg = i;
		}

		if (chimBeg < 0)
			continue;

		for (j = 0; j <= trim_ae / ctx->twidth && numChim == 0; j++)
		{
			if (alnBeg[j] >= BORDER_COV)
			{
				chimEnd = j;
			}

			else if (alnBeg[j + 1] >= BORDER_COV)
			{
				chimEnd = ++j;
			}

			else if (alnBeg[j] + alnBeg[j + 1] >= BORDER_COV)
			{
				chimEnd = j + 1;
			}

			if (chimEnd < 0)
				continue;

#ifdef DEBUG_CHIMER
			printf("check putative chimer indexes [%d, %d]\n", chimBeg, chimEnd);
#endif

			chimBeg = i * ctx->twidth;
			chimEnd *= ctx->twidth;

			if (chimBeg > chimEnd)
			{
				int tmp = chimBeg;
				chimBeg = chimEnd;
				chimEnd = tmp;
			}

#ifdef DEBUG_CHIMER
			printf("check putative chimer interval [%d, %d]\n", chimBeg, chimEnd);
#endif

			if((chimEnd - chimBeg) > MAX_CHIMER_LEN)
			{
#ifdef DEBUG_CHIMER
				printf("chimer interval too large [%d, %d] len %d-- IGNORED\n", chimBeg, chimEnd, (chimEnd - chimBeg));
#endif
				chimEnd = -1;
				continue;

			}

			int tmpRep = getRepeatCount(ctx, ovls->aread, chimBeg, chimEnd);
#ifdef DEBUG_CHIMER
			printf("getRepeatCount %d %d %d: %d\n", ovls->aread, chimBeg, chimEnd, tmpRep);
#endif

			if ((chimEnd - chimBeg) > 1000 && tmpRep * 1.0 / (chimEnd - chimBeg) < 0.7)
			{
#ifdef DEBUG_CHIMER
				printf("chimer interval too large [%d, %d] len %d but not repetitive %.2f%%-- IGNORED\n", chimBeg, chimEnd, (chimEnd - chimBeg),
						tmpRep * 1.0 / (chimEnd - chimBeg));
#endif
				chimEnd = -1;
				continue;
			}

			int intersectBreads = oChainIntervalIntersection(ctx, ovls, novl, chimBeg, chimEnd, trim_ab, trim_ae);
#ifdef DEBUG_CHIMER
			printf("intersectBreads %d\n", intersectBreads);
#endif
			if (intersectBreads < 1)
			{
				++numChim;

#ifdef DEBUG_CHIMER
				printf("FOUND CHIMER %d [%d, %d]\n", ovls->aread, chimBeg, chimEnd);
#endif
				break;
			}
			chimEnd = -1;

			j+=2; // don't check neighboring walls
		}

		chimBeg = -1;
		i+=2; // don't check neighboring walls
	}

	if (numChim)
	{
		// discard all overlaps from current aread
		// write and aread id to file
		for (i = 0; i < novl; i++)
		{
			ovls[i].flags |= (OVL_GAP | OVL_DISCARD);
		}
#ifdef DEBUG_CHIMER
		fprintf(ctx->chimer, "%d\n", ovls->aread);
#endif
#ifdef VERBOSE
		printf("CHIMER %d [%d, %d] - discarded overlaps %d\n", ovls->aread, chimBeg, chimEnd, novl);
#endif
		ctx->stats_chimers++;
		ctx->stats_chimers_novl += novl;
	}

	free(alnBeg);
	free(alnEnd);

	return 1;
}

static int handle_repeat_chimers(GapContext* ctx, Overlap* ovls, int novl)
{
	HITS_TRACK* trackRepeat = ctx->trackRepeat;
	int a = ovls->aread;

	int trim_left, trim_right;
	get_trim(ctx->db, ctx->trackTrim, a, &trim_left, &trim_right);

	track_anno* ta = trackRepeat->anno;
	track_data* td = trackRepeat->data;

	track_anno tab = ta[a] / sizeof(track_data);
	track_anno tae = ta[a + 1] / sizeof(track_data);

	while (tab < tae)
	{
		track_data rb = td[tab];
		track_data re = td[tab + 1];

		while (tab < tae - 2 && re + 100 > td[tab + 2])
		{
			tab += 2;
			re = td[tab + 1];
		}

		tab += 2;

		if (rb - 100 <= trim_left || re + 100 >= trim_right)
		{
			continue;
		}

		int cross_l, cross_r;
		int terminate_l, terminate_r;
		cross_l = cross_r = terminate_l = terminate_r = 0;

		int i;
		for (i = 0; i < novl; i++)
		{
			Overlap* ovl = ovls + i;

			if (ovl->path.abpos < rb - 100)
			{
				if (ovl->path.aepos > re - 100 && ovl->path.aepos < re + 100)
				{
					terminate_r += 1;
				}
				else if (ovl->path.aepos > re + 100)
				{
					cross_r += 1;
				}
			}

			if (ovl->path.aepos > re + 100)
			{
				if (ovl->path.abpos > rb - 100 && ovl->path.abpos < rb + 100)
				{
					terminate_l += 1;
				}
				else if (ovl->path.abpos < rb - 100)
				{
					cross_l += 1;
				}
			}
		}

		if (cross_l < 2 && cross_r < 2 && terminate_l > 2 && terminate_r > 2)
		{
			drop_break(ovls, novl, rb, re);

#ifdef DEBUG_CGAPS
			printf("%7d %5d..%5d CL %3d TL %3d CR %3d TR %3d\n", a, rb, re, cross_l, terminate_l, cross_r, terminate_r);
#endif
		}

	}

	return 1;
}

static int handle_gaps_and_breaks(GapContext* ctx, Overlap* ovls, int novl)
{
	int stitch = ctx->stitch;
	HITS_TRACK* trackExclude = ctx->trackExclude;

	bzero(ctx->rm_bins, sizeof(uint64) * ctx->rm_maxbins);

	int trim_b = INT_MAX;
	int trim_e = 0;

	int i;
	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;

		if ((ovl->flags & OVL_DISCARD) || (ovl->aread == ovl->bread))
		{
			continue;
		}

		int j;
		for (j = i + 1; j < novl && ovls[j].bread == ovl->bread; j++)
		{
			if (ovl_intersect(ovl, ovls + j))
			{
				if (ovl->path.aepos - ovl->path.abpos < ovls[j].path.aepos - ovls[j].path.abpos)
				{
					ovl->flags |= OVL_TEMP;
				}
				else
				{
					ovls[j].flags |= OVL_TEMP;
				}
			}
		}

		if (ovl->flags & OVL_TEMP)
		{
			continue;
		}

		int b = (ovl->path.abpos + MIN_LR) / BIN_SIZE;
		int e = (ovl->path.aepos - MIN_LR) / BIN_SIZE;

		trim_b = MIN(ovl->path.abpos, trim_b);
		trim_e = MAX(ovl->path.aepos, trim_e);

		while (b < e)
		{
			ctx->rm_bins[b]++;
			b++;
		}
	}

	if (trim_b >= trim_e)
	{
		return 1;
	}

	// printf("trim %d..%d\n", trim_b, trim_e);

	int b = (trim_b + MIN_LR) / BIN_SIZE;
	int e = (trim_e - MIN_LR) / BIN_SIZE;

	int beg = -1;

	while (b < e)
	{
		// if (ctx->rm_bins[b] != 0)
		if (ctx->rm_bins[b] > 1)
		{
			if (beg != -1)
			{
				int breakb = (beg - 1) * BIN_SIZE;
				int breake = (b + 1) * BIN_SIZE;

#ifdef DEBUG_GAPS
				printf("READ %7d BREAK %3d..%3d ", ovls->aread, beg, b);
#endif
				// break covered by track interval

				int skip = 0;

				if (trackExclude)
				{
					track_anno* ta = trackExclude->anno;
					track_data* td = trackExclude->data;

					track_anno tab = ta[ovls->aread] / sizeof(track_data);
					track_anno tae = ta[ovls->aread + 1] / sizeof(track_data);

					int masked = 0;

					track_data maskb, maske;

					while (tab < tae)
					{
						maskb = td[tab];
						maske = td[tab + 1];

						if (breakb > maskb && breake < maske)
						{
							masked = 1;
							break;
						}

						tab += 2;
					}

					if (masked)
					{
#ifdef DEBUG_GAPS
						printf(" MASKED %5d..%5d\n", maskb, maske);
#endif
						skip = 1;
					}
				}

				if (!skip && stitch > 0)
				{
					int nstitch = stitchable(ovls, novl, stitch, breakb, breake);

					if (nstitch)
					{
#ifdef DEBUG_GAPS
						printf(" STITCHABLE using %d\n", nstitch);
#endif
						skip = 1;
					}
				}

				if (!skip)
				{
#ifdef DEBUG_GAPS
					printf("\n");
#endif
					ctx->stats_breaks += 1;
					ctx->stats_breaks_novl += drop_break(ovls, novl, breakb, breake);
				}

				beg = -1;
			}
		}
		else
		{
			if (beg == -1)
			{
				beg = b;
			}
		}

		b++;
	}

	if (beg != -1)
	{
		int breakb = (beg - 1) * BIN_SIZE;
		int breake = (b + 1) * BIN_SIZE;

		ctx->stats_breaks += 1;
		ctx->stats_breaks_novl += drop_break(ovls, novl, breakb, breake);
	}

	return 1;
}

static int drop_containments(GapContext* ctx, Overlap* ovl, int novl)
{
	if (novl < 2)
	{
		return 1;
	}

	int i;
	int ab1, ab2, ae1, ae2;
	int bb1, bb2, be1, be2;

	for (i = 0; i < novl; i++)
	{
		int bread = ovl[i].bread;

		if ((ovl[i].flags & OVL_DISCARD) || ovl[i].aread == bread)
		{
			continue;
		}

		int blen = DB_READ_LEN(ctx->db, bread);

		ab1 = ovl[i].path.abpos;
		ae1 = ovl[i].path.aepos;

		if (ovl[i].flags & OVL_COMP)
		{
			bb1 = blen - ovl[i].path.bepos;
			be1 = blen - ovl[i].path.bbpos;
		}
		else
		{
			bb1 = ovl[i].path.bbpos;
			be1 = ovl[i].path.bepos;
		}

		int k;
		for (k = i + 1; k < novl; k++)
		{
			if (ovl[k].flags & OVL_DISCARD)
			{
				continue;
			}

			if (ovl[k].bread != bread)
			{
				break;
			}

			ab2 = ovl[k].path.abpos;
			ae2 = ovl[k].path.aepos;

			if (ovl[k].flags & OVL_COMP)
			{
				bb2 = blen - ovl[k].path.bepos;
				be2 = blen - ovl[k].path.bbpos;
			}
			else
			{
				bb2 = ovl[k].path.bbpos;
				be2 = ovl[k].path.bepos;
			}

			int cont = contained(ab1, ae1, ab2, ae2);
			if (cont && contained(bb1, be1, bb2, be2))
			{
#ifdef VERBOSE_CONTAINMENT
				printf("CONTAINMENT %8d @ %5d..%5d -> %8d @ %5d..%5d\n"
						"                       %5d..%5d -> %8d @ %5d..%5d\n",
						a, ab1, ae1, ovl[i].bread, bb1, be1,
						ab2, ae2, ovl[k].bread, bb2, be2);
#endif

				if (cont == 1)
				{
					ovl[i].flags |= OVL_DISCARD | OVL_CONT;
					ctx->stats_contained++;
					break;
				}
				else if (cont == 2)
				{
					ovl[k].flags |= OVL_DISCARD | OVL_CONT;
					ctx->stats_contained++;
				}
			}
		}
	}

	return 1;
}

static void gaps_pre(PassContext* pctx, GapContext* ctx)
{
#ifdef VERBOSE
	printf(ANSI_COLOR_GREEN "PASS gaps" ANSI_COLOR_RESET "\n");
#endif
	ctx->twidth = pctx->twidth;

	ctx->rm_maxbins = ( DB_READ_MAXLEN(ctx->db) + BIN_SIZE) / BIN_SIZE;
	ctx->rm_bins = malloc(sizeof(uint64) * ctx->rm_maxbins);

	ctx->rc_left = NULL;
	ctx->rc_maxleft = 0;

	ctx->rc_right = NULL;
	ctx->rc_maxright = 0;

	// trim

	ctx->trim = trim_init(ctx->db, pctx->twidth, ctx->trackTrim, ctx->rl);

	// chimer detection

	if (ctx->chimer != NULL)
	{
		ctx->curChains = 0;
		ctx->maxChains = 5;
		ctx->ovlChains = (Chain*) malloc(sizeof(Chain) * ctx->maxChains);
		bzero(ctx->ovlChains, sizeof(Chain) * ctx->maxChains);
	}
}

static void gaps_post(GapContext* ctx)
{
#ifdef VERBOSE
	if (ctx->trackTrim)
	{
		printf("%13lld of %13lld overlaps trimmed\n", ctx->trim->nTrimmedOvls, ctx->trim->nOvls);
		printf("%13lld of %13lld bases trimmed\n", ctx->trim->nTrimmedBases, ctx->trim->nOvlBases);
	}

	printf("dropped %d containments\n", ctx->stats_contained);
	printf("dropped %d overlaps in %d break/gaps\n", ctx->stats_breaks_novl, ctx->stats_breaks);

	if (ctx->chimer)
	{
		printf("dropped %d overlaps in %d chimers\n", ctx->stats_chimers_novl, ctx->stats_chimers);
	}
#endif

	free(ctx->rm_bins);

	free(ctx->rc_left);
	free(ctx->rc_right);

	if (ctx->trackTrim)
	{
		trim_close(ctx->trim);
	}

	if (ctx->chimer)
	{
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
}

static int gaps_handler(void* _ctx, Overlap* ovl, int novl)
{
	GapContext* ctx = (GapContext*) _ctx;

	// trim
	if (ctx->trackTrim)
	{
		int i;
		for (i = 0; i < novl; i++)
			trim_overlap(ctx->trim, ovl + i);
	}

	drop_containments(ctx, ovl, novl);

	handle_gaps_and_breaks(ctx, ovl, novl);

//	if (ctx->trackRepeat)
//	{
//		handle_repeat_chimers(ctx, ovl, novl);
//	}
	if (ctx->chimer)
	{
		handle_chimers(ctx, ovl, novl);
	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "usage:   [-s <int>] [-pL] [-etR <track>] [-C file] <db> <ovl.in> <ovl.out>\n");
	fprintf(stderr, "options: -s ... stitch distance (%d)\n", DEF_ARG_S);
	fprintf(stderr, "         -p ... purge discarded overlaps\n");
	fprintf(stderr, "         -e ... exclude gaps in track intervals\n");
	fprintf(stderr, "         -t ... trim overlaps before gap detection\n");
	fprintf(stderr, "         -L ... two pass processing with read caching\n");

	fprintf(stderr, "EXPERIMENTAL: \n");
	fprintf(stderr, "         -C ... discard chimeric reads and overlaps. Write chimer reads IDs into <file>.\n");
	fprintf(stderr, "         -R ... repeat track - used to detect chimers in large repeats\n");

}

int main(int argc, char* argv[])
{
	HITS_DB db;
	PassContext* pctx;
	FILE* fileOvlIn;
	FILE* fileOvlOut;
	GapContext gctx;

	// process arguments
	int arg_purge = DEF_ARG_P;
	char* arg_trimTrack = NULL;
	char* excludeTrack = NULL;
	char* repeatTrack = NULL;
	char* chimerFile = NULL;

	bzero(&gctx, sizeof(GapContext));

	gctx.stitch = DEF_ARG_S;

	opterr = 0;

	int c;
	while ((c = getopt(argc, argv, "Ls:pe:R:t:C:")) != -1)
	{
		switch (c)
		{
		case 'R':
			repeatTrack = optarg;
			break;

		case 'e':
			excludeTrack = optarg;
			break;

		case 'C':
			chimerFile = optarg;
			break;

		case 't':
			arg_trimTrack = optarg;
			break;

		case 's':
			gctx.stitch = atoi(optarg);
			break;

		case 'p':
			arg_purge = 1;
			break;

		case 'L':
			gctx.useRLoader = 1;
			break;

		default:
			usage();
			exit(1);
		}
	}

	if (argc - optind < 3)
	{
		usage();
		exit(1);
	}

	char* pcPathReadsIn = argv[optind++];
	char* pcPathOverlapsIn = argv[optind++];
	char* pcPathOverlapsOut = argv[optind++];

	if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL)
	{
		fprintf(stderr, "could not open input track '%s'\n", pcPathOverlapsIn);
		exit(1);
	}

	if ((fileOvlOut = fopen(pcPathOverlapsOut, "w")) == NULL)
	{
		fprintf(stderr, "could not open output track '%s'\n", pcPathOverlapsOut);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "could not open database '%s'\n", pcPathReadsIn);
		exit(1);
	}

	if (excludeTrack)
	{
		gctx.trackExclude = track_load(&db, excludeTrack);

		if (!gctx.trackExclude)
		{
			fprintf(stderr, "could not open track '%s'\n", excludeTrack);
			exit(1);
		}
	}

	if (chimerFile)
	{
		gctx.chimer = fopen(chimerFile, "w");

		if (gctx.chimer == NULL)
		{
			fprintf(stderr, "could not open chimer file '%s'\n", chimerFile);
			exit(1);
		}

		if (repeatTrack == NULL)
		{
			fprintf(stderr, "to detect chimers in large repeat regions a repeat track is necessary! \n");
			exit(1);
		}

		if (arg_trimTrack == NULL)
		{
			fprintf(stderr, "to detect chimers in repeat regions a trim track is necessary! \n");
			exit(1);
		}
	}

	if (repeatTrack)
	{
		gctx.trackRepeat = track_load(&db, repeatTrack);

		if (!gctx.trackRepeat)
		{
			fprintf(stderr, "could not open track '%s'\n", repeatTrack);
			exit(1);
		}
	}

	gctx.db = &db;

	if (arg_trimTrack != NULL)
	{
		gctx.trackTrim = track_load(gctx.db, arg_trimTrack);
		if (!gctx.trackTrim)
		{
			fprintf(stderr, "could not open track '%s'\n", arg_trimTrack);
			exit(1);
		}

		if (gctx.useRLoader)
		{
			gctx.rl = rl_init(&db, 1);

			pctx = pass_init(fileOvlIn, NULL);

			pctx->data = &gctx;
			pctx->split_b = 1;
			pctx->load_trace = 0;

			pass(pctx, loader_handler);
			rl_load_added(gctx.rl);
			pass_free(pctx);
		}
	}

	pctx = pass_init(fileOvlIn, fileOvlOut);

	pctx->split_b = 0;
	pctx->load_trace = 1;
	pctx->unpack_trace = (gctx.trackTrim == NULL ? 0 : 1);
	pctx->data = &gctx;
	pctx->write_overlaps = 1;
	pctx->purge_discarded = arg_purge;

	gaps_pre(pctx, &gctx);

	pass(pctx, gaps_handler);

	gaps_post(&gctx);

	// cleanup

	Close_DB(&db);

	if (gctx.useRLoader)
	{
		rl_free(gctx.rl);
	}

	pass_free(pctx);

	fclose(fileOvlIn);
	fclose(fileOvlOut);
	if (gctx.chimer)
		fclose(gctx.chimer);

	return 0;
}

