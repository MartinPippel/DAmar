/*
 * LAseparate.c
 *
 *  Created on: 10 Jul 2018
 *      Author: pippel
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/param.h>
#include <unistd.h>

#include "../dalign/align.h"
#include "../db/DB.h"
#include "../lib/colors.h"
#include "../lib/oflags.h"
#include "../lib/pass.h"
#include "../lib/read_loader.h"
#include "../lib/tracks.h"
#include "../lib/utils.h"

#define DEF_ARG_L 0
#define DEF_ARG_O 0
#define DEF_ARG_R "rep"

#undef DEBUG_CHAIN

static void usage()
{
	fprintf(stderr, "[-vLT] [-lo <int>] [-r <track>] <db> <overlaps_in> <overlaps_out> <overlaps_out_discard>\n");

	fprintf(stderr, "options: -v ... verbose\n");
	fprintf(stderr, "         -L ... two pass processing with read caching\n");
	fprintf(stderr, "         -o ... min overlap length (default: %d)\n", DEF_ARG_O);
	fprintf(stderr, "         -l ... min read length (default: %d)\n", DEF_ARG_L);
	fprintf(stderr, "         -r ... repeat track name (default %s)\n", DEF_ARG_R);
	fprintf(stderr, "         -T ... separate Type: 0 - for RepComp, 1 - for ForceAlign\n");
}

typedef struct
{
	Overlap **ovls;
	int novl;
	int maxOvl;
} Chain;

typedef struct
{
	// settings
	int nMinAlnLength;
	int nMinReadLength;
	int nVerbose;
	int type;

	// status
	int nKeptOvls;
	int nDiscardOvls;

	HITS_DB* db;
	HITS_TRACK* trackRepeat;

	int useRLoader;
	Read_Loader* rl;

	ovl_header_twidth twidth;

	Chain *ovlChains;
	int curChains;
	int maxChains;

} SeparateContext;

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

static int contained(int ab, int ae, int bb, int be)
{
	if (ab >= bb && ae <= be)
	{
		return 1;
	}

	return 0;
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

static int getRepeatBases(SeparateContext *ctx, Overlap *ovl, int read)
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

static int loader_handler(void* _ctx, Overlap* ovl, int novl)
{
	SeparateContext* ctx = (SeparateContext*) _ctx;
	Read_Loader* rl = ctx->rl;

	int i;
	for (i = 0; i < novl; i++)
	{
		int b = ovl[i].bread;

		int trim_b_left, trim_b_right;

		trim_b_left = 0;
		trim_b_right = DB_READ_LEN(ctx->db, b);

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

static void separate_pre(PassContext* pctx, SeparateContext* fctx)
{
	printf( ANSI_COLOR_GREEN "PASS separate\n" ANSI_COLOR_RESET);

	fctx->twidth = pctx->twidth;
	fctx->nKeptOvls = 0;
	fctx->nDiscardOvls = 0;

	fctx->curChains = 0;
	fctx->maxChains = 5;
	fctx->ovlChains = (Chain*) malloc(sizeof(Chain) * fctx->maxChains);
	bzero(fctx->ovlChains, sizeof(Chain) * fctx->maxChains);
}

static void separate_post(SeparateContext* ctx)
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

	printf("Kept    ovls: %10d\n", ctx->nKeptOvls);
	printf("Discard ovls: %10d\n", ctx->nDiscardOvls);
}

static void chain(SeparateContext*ctx, Overlap *ovls, int n)
{
	/// TODO hard coded - was taken from LAanalyze which is adapted for contig vs contig alignments
	int MIN_OVL_LOOKAHEAD = 5000;
	int MAX_OVL_LOOKAHEAD = 30000;
	int STRIDE_OVL_LOOKAHEAD = 5000;

	if (n < 2)
		return;
#ifdef DEBUG_CHAIN
	printf("chain(%d,%d,%d) CHAIN: n%d m%d\n", ovls->aread, ovls->bread, n, ctx->curChains, ctx->maxChains);
#endif
	int aread, bread;
	int alen, blen;
	int i;

	aread = ovls->aread;
	bread = ovls->bread;

	alen = DB_READ_LEN(ctx->db, aread);
	blen = DB_READ_LEN(ctx->db, bread);

	int nremain = n;

// mark contained overlaps
#ifdef DEBUG_CHAIN
	printf("mark contained overlaps\n");
#endif

	{
		int j;
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = ovls + i;

			if (ovl_i->flags & (OVL_CONT))
				continue;

			for (j = i + 1; j < n; j++)
			{
				Overlap *ovl_j = ovls + j;

				if (ovl_j->flags & (OVL_CONT))
					continue;

				if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos)
						&& contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
				{
					ovl_j->flags |= (OVL_CONT);
				}
			}

		}
	}

	{
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = ovls + i;

			if (ovl_i->flags & (OVL_CONT | OVL_DISCARD))
			{
				ovl_i->flags |= OVL_DISCARD;
				nremain--;
			}
		}
	}

#ifdef DEBUG_CHAIN
	printf("nremain %d\n", nremain);
#endif
	if (nremain < 1)
		return;

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

			if (ovl->flags & (OVL_CONT | OVL_DISCARD | OVL_TEMP))
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
				if ((ovl->flags & (OVL_TEMP | OVL_CONT | OVL_DISCARD)) || ((ovl->flags & OVL_COMP) != (chain->ovls[chainIdx]->flags & OVL_COMP)))
					continue;

				if (ovl->path.abpos < chain->ovls[chainIdx]->path.abpos)
					continue;

				if (ovl->path.abpos > chain->ovls[chainLastIdx]->path.abpos)
					break;

				int ovhBases = 100;
				for (j = chainIdx; j < chainLastIdx; j++)
				{
					if (chain->ovls[j]->path.aepos - ovhBases < ovl->path.abpos && chain->ovls[j + 1]->path.abpos + ovhBases > ovl->path.aepos
							&& chain->ovls[j]->path.bepos - ovhBases < ovl->path.bbpos && chain->ovls[j + 1]->path.bbpos + ovhBases > ovl->path.bepos)
					{
						Overlap *lastAddedOvl = chain->ovls[chain->novl - 1];

						if (intersect(ovl->path.abpos, ovl->path.aepos, lastAddedOvl->path.abpos, lastAddedOvl->path.aepos) > ovhBases
								|| intersect(ovl->path.bbpos, ovl->path.bepos, lastAddedOvl->path.bbpos, lastAddedOvl->path.bepos) > ovhBases)
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
			// mark remaining ovls as DISC if the overlap with a chain overlap !!
#ifdef DEBUG_CHAIN
			printf("// mark remaining ovls as DISC if they overlap with a chain overlap !!\n");
#endif
			int chainIdx = 0;
			int chainLastIdx = chain->novl - 1;
			int j;
			for (i = 0; i < n; i++)
			{
				Overlap *ovl = ovls + i;
				if ((ovl->flags & (OVL_TEMP | OVL_CONT | OVL_DISCARD)) || ((ovl->flags & OVL_COMP) != (chain->ovls[chainIdx]->flags & OVL_COMP)))
					continue;

				for (j = chainIdx; j <= chainLastIdx; j++)
				{
					if (intersect(chain->ovls[j]->path.abpos, chain->ovls[j]->path.aepos, ovl->path.abpos, ovl->path.aepos)
							|| intersect(chain->ovls[j]->path.bbpos, chain->ovls[j]->path.bepos, ovl->path.bbpos, ovl->path.bepos))
					{
						ovl->flags |= OVL_DISCARD;
						nremain--;
#ifdef DEBUG_CHAIN
						printf("DISCARD [%d, %d] [%d, %d] nremain %d\n", ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, nremain);
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
			int j;
			for (j = 0; j < chain->novl; j++)
			{
				chain->ovls[j]->flags |= OVL_DISCARD;
			}
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

static int separate_handler(void* _ctx, Overlap* ovl, int novl)
{
	SeparateContext* ctx = (SeparateContext*) _ctx;
	int i;

	chain(ctx, ovl, novl);

	if (ctx->type == 0) // separate for repcomp
	{
		if (ctx->curChains)
		{
			Chain *bestChain = ctx->ovlChains;

			int properBeg = 0;
			int properEnd = 0;
			int gapBasesInA = 0;
			int gapBasesInB = 0;
			int itsBasesInA = 0;
			int itsBasesInB = 0;
			int overlapBases = 0;

			// check for proper begin
			if (MIN(bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bbpos) < 1000)
				properBeg = 1;
			// check for proper end
			if (bestChain->ovls[bestChain->novl - 1]->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread)
					|| bestChain->ovls[bestChain->novl - 1]->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread))
				properEnd = 1;
#ifdef DEBUG
			printf("properBeg: %d, properEnd %d\n", properBeg, properEnd);
#endif
			if (properBeg && properEnd)
			{
				int i, its_a, its_b;
				overlapBases = MAX(bestChain->ovls[0]->path.aepos - bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bepos - bestChain->ovls[0]->path.bbpos);

				for (i = 1; i < bestChain->novl; i++)
				{
					its_a = its_b = 0;
					overlapBases += MAX(bestChain->ovls[i]->path.aepos - bestChain->ovls[i]->path.abpos, bestChain->ovls[i]->path.bepos - bestChain->ovls[i]->path.bbpos);
					// check for intersection in A
					if (bestChain->ovls[i]->path.abpos < bestChain->ovls[i - 1]->path.aepos)
					{
						its_a = bestChain->ovls[i - 1]->path.aepos - bestChain->ovls[i]->path.abpos;
						if (its_a > 1000)
						{
							itsBasesInA = -1;
							break;
						}
						itsBasesInA += its_a;
					}
					// check for gap in A
					else
					{
						int gap = bestChain->ovls[i]->path.abpos - bestChain->ovls[i - 1]->path.aepos;
//						if (gap > 1000)
//						{
//							gapBasesInA = -1;
//							break;
//						}
						gapBasesInA += gap;
					}
					// check for intersection in B
					if (bestChain->ovls[i]->path.bbpos < bestChain->ovls[i - 1]->path.bepos)
					{
						its_b = bestChain->ovls[i - 1]->path.bepos - bestChain->ovls[i]->path.bbpos;
						if (its_b > 1000)
						{
							itsBasesInB = -1;
							break;
						}
						itsBasesInB += its_b;
					}
					// check for gap in B
					else
					{
						int gap = bestChain->ovls[i]->path.bbpos - bestChain->ovls[i - 1]->path.bepos;
//						if (gap > 1000)
//						{
//							gapBasesInB = -1;
//							break;
//						}
						gapBasesInB += gap;
					}
					overlapBases -= MAX(its_a, its_b);
				}
			}

#ifdef DEBUG
			printf("if(%d && %d && %d >=0 && %d >=0 && %d >=0 && %d >=0 && %d * 0.3 > MAX(%d, %d))\n", properBeg, properEnd, itsBasesInA, itsBasesInB, gapBasesInA,
					gapBasesInB, overlapBases, gapBasesInA, gapBasesInB);
#endif
			// if there is a proper chain between A and B reads, then discard all overlaps between A and B for the repcomp step, (otherwise do repcomp)
			if (properBeg && properEnd && itsBasesInA >= 0 && itsBasesInB >= 0 && gapBasesInA >= 0 && gapBasesInB >= 0
					&& overlapBases * 0.3 > MAX(gapBasesInA, gapBasesInB))
			{
#ifdef DEBUG
				printf("FOUND PROPER CHAIN - EXCLUDE ALL OVLS FROM REPCOMP INPUT\n");
#endif
				for (i = 0; i < novl; i++)
				{
					ovl[i].flags |= OVL_DISCARD;
					ctx->nDiscardOvls++;
				}
			}
			else // reset all discard flags that were set during chain detection
			{
				for (i = 0; i < novl; i++)
				{
					if (ovl[i].flags & OVL_DISCARD)
					{
						ovl[i].flags &= ~OVL_DISCARD;
					}
					ctx->nKeptOvls++;
				}
			}
		}
		else // there is only one overlap !!! check this
		{
			int properBeg = (MIN(ovl->path.abpos, ovl->path.bbpos) < 1000);
			int properEnd = (ovl->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread) || ovl->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread));

			if (properBeg && properEnd)
			{
#ifdef DEBUG
				printf("FOUND PROPER SINGLE OVL - EXCLUDE OVL FROM REPCOMP INPUT\n");
#endif
				ovl->flags |= OVL_DISCARD;
				ctx->nDiscardOvls++;
			}
			else
			{
				ctx->nKeptOvls++;
			}
		}
	}
	else if (ctx->type == 1)
	{
		if (ctx->curChains == 1) // ignore everything, where you have more then 1 chain
		{
			Chain *bestChain = ctx->ovlChains;

			int properBeg = 0;
			int properEnd = 0;
			int gapBasesInA = 0;
			int gapBasesInB = 0;
			int itsBasesInA = 0;
			int itsBasesInB = 0;
			int overlapBases = 0;

			// check for proper begin
			if (MIN(bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bbpos) < 1000)
				properBeg = 1;
			// check for proper end
			if (bestChain->ovls[bestChain->novl - 1]->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread)
					|| bestChain->ovls[bestChain->novl - 1]->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread))
				properEnd = 1;
#ifdef DEBUG
			printf("properBeg: %d, properEnd %d\n", properBeg, properEnd);
#endif
			if (properBeg && properEnd)
			{
				int i, its_a, its_b;
				overlapBases = MAX(bestChain->ovls[0]->path.aepos - bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bepos - bestChain->ovls[0]->path.bbpos);

				for (i = 1; i < bestChain->novl; i++)
				{
					its_a = its_b = 0;
					overlapBases += MAX(bestChain->ovls[i]->path.aepos - bestChain->ovls[i]->path.abpos, bestChain->ovls[i]->path.bepos - bestChain->ovls[i]->path.bbpos);
					// check for intersection in A
					if (bestChain->ovls[i]->path.abpos < bestChain->ovls[i - 1]->path.aepos)
					{
						its_a = bestChain->ovls[i - 1]->path.aepos - bestChain->ovls[i]->path.abpos;
						if (its_a > 1000)
						{
							itsBasesInA = -1;
							break;
						}
						itsBasesInA += its_a;
					}
					// check for gap in A
					else
					{
						int gap = bestChain->ovls[i]->path.abpos - bestChain->ovls[i - 1]->path.aepos;
//						if (gap > 1000)
//						{
//							gapBasesInA = -1;
//							break;
//						}
						gapBasesInA += gap;
					}
					// check for intersection in B
					if (bestChain->ovls[i]->path.bbpos < bestChain->ovls[i - 1]->path.bepos)
					{
						its_b = bestChain->ovls[i - 1]->path.bepos - bestChain->ovls[i]->path.bbpos;
						if (its_b > 1000)
						{
							itsBasesInB = -1;
							break;
						}
						itsBasesInB += its_b;
					}
					// check for gap in B
					else
					{
						int gap = bestChain->ovls[i]->path.bbpos - bestChain->ovls[i - 1]->path.bepos;
//						if (gap > 1000)
//						{
//							gapBasesInB = -1;
//							break;
//						}
						gapBasesInB += gap;
					}
					overlapBases -= MAX(its_a, its_b);
				}
			}

#ifdef DEBUG
			printf("if(%d && %d && %d >=0 && %d >=0 && %d >=0 && %d >=0 && %d * 0.3 > MAX(%d, %d))\n", properBeg, properEnd, itsBasesInA, itsBasesInB, gapBasesInA,
					gapBasesInB, overlapBases, gapBasesInA, gapBasesInB);
#endif
			if (properBeg && properEnd && itsBasesInA >= 0 && itsBasesInB >= 0 && gapBasesInA >= 0 && gapBasesInB >= 0
					&& overlapBases * 0.3 > MAX(gapBasesInA, gapBasesInB))
			{
#ifdef DEBUG
				printf("FOUND PROPER CHAIN - INCLUDE ALL CHAIN OVLS FOR FORCEALIGN INPUT\n");
#endif
				for (i = 0; i < novl; i++)
				{
					if (ovl[i].flags & OVL_DISCARD)
						ctx->nDiscardOvls++;
					else
						ctx->nKeptOvls++;
				}
			}
			else // multiple chains: discard all overlaps, i.e. no forcealign
			{
				for (i = 0; i < novl; i++)
				{
					ovl[i].flags |= OVL_DISCARD;
					ctx->nDiscardOvls++;
				}
			}
		}
		else if (ctx->curChains > 1) // multiple chains: discard all overlaps, i.e. no forcealign
		{
			for (i = 0; i < novl; i++)
			{
				ovl[i].flags |= OVL_DISCARD;
				ctx->nDiscardOvls++;
			}
		}
		else // there is only one overlap !!! check this
		{
			int properBeg = (MIN(ovl->path.abpos, ovl->path.bbpos) < 1000);
			int properEnd = (ovl->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread) || ovl->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread));

			int realBeg = (MIN(ovl->path.abpos, ovl->path.bbpos) == 0);
			int realEnd = (ovl->path.aepos == DB_READ_LEN(ctx->db, ovl->aread) || ovl->path.bepos == DB_READ_LEN(ctx->db, ovl->bread));

			if ((properBeg && properEnd) && !(realBeg && realEnd))
			{
#ifdef DEBUG
				printf("FOUND PROPER SINGLE OVL - INCLUDE OVL FOR FORCEALIGN INPUT\n");
#endif
				ctx->nKeptOvls++;
			}
			else
			{
				ovl->flags |= OVL_DISCARD;
				ctx->nDiscardOvls++;
			}
		}
	}

	// reset chain and ovl counter
	for (i = 0; i < ctx->curChains; i++)
		ctx->ovlChains[i].novl = 0;
	ctx->curChains = 0;

	return 1;
}

static int separateDiscard_handler(void* _ctx, Overlap* ovl, int novl)
{
	SeparateContext* ctx = (SeparateContext*) _ctx;
	int i;

	chain(ctx, ovl, novl);

	if (ctx->type == 0) // separate for repcomp
	{
		if (ctx->curChains)
		{
			Chain *bestChain = ctx->ovlChains;

			int properBeg = 0;
			int properEnd = 0;
			int gapBasesInA = 0;
			int gapBasesInB = 0;
			int itsBasesInA = 0;
			int itsBasesInB = 0;
			int overlapBases = 0;

			// check for proper begin
			if (MIN(bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bbpos) < 1000)
				properBeg = 1;
			// check for proper end
			if (bestChain->ovls[bestChain->novl - 1]->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread)
					|| bestChain->ovls[bestChain->novl - 1]->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread))
				properEnd = 1;
#ifdef DEBUG
			printf("properBeg: %d, properEnd %d\n", properBeg, properEnd);
#endif
			if (properBeg && properEnd)
			{
				int i, its_a, its_b;
				overlapBases = MAX(bestChain->ovls[0]->path.aepos - bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bepos - bestChain->ovls[0]->path.bbpos);

				for (i = 1; i < bestChain->novl; i++)
				{
					its_a = its_b = 0;
					overlapBases += MAX(bestChain->ovls[i]->path.aepos - bestChain->ovls[i]->path.abpos, bestChain->ovls[i]->path.bepos - bestChain->ovls[i]->path.bbpos);
					// check for intersection in A
					if (bestChain->ovls[i]->path.abpos < bestChain->ovls[i - 1]->path.aepos)
					{
						its_a = bestChain->ovls[i - 1]->path.aepos - bestChain->ovls[i]->path.abpos;
						if (its_a > 1000)
						{
							itsBasesInA = -1;
							break;
						}
						itsBasesInA += its_a;
					}
					// check for gap in A
					else
					{
						int gap = bestChain->ovls[i]->path.abpos - bestChain->ovls[i - 1]->path.aepos;
//						if (gap > 1000)
//						{
//							gapBasesInA = -1;
//							break;
//						}
						gapBasesInA += gap;
					}
					// check for intersection in B
					if (bestChain->ovls[i]->path.bbpos < bestChain->ovls[i - 1]->path.bepos)
					{
						its_b = bestChain->ovls[i - 1]->path.bepos - bestChain->ovls[i]->path.bbpos;
						if (its_b > 1000)
						{
							itsBasesInB = -1;
							break;
						}
						itsBasesInB += its_b;
					}
					// check for gap in B
					else
					{
						int gap = bestChain->ovls[i]->path.bbpos - bestChain->ovls[i - 1]->path.bepos;
//						if (gap > 1000)
//						{
//							gapBasesInB = -1;
//							break;
//						}
						gapBasesInB += gap;
					}
					overlapBases -= MAX(its_a, its_b);
				}
			}
#ifdef DEBUG
			printf("if(%d && %d && %d >=0 && %d >=0 && %d >=0 && %d >=0 && %d * 0.3 > MAX(%d, %d))\n", properBeg, properEnd, itsBasesInA, itsBasesInB, gapBasesInA,
					gapBasesInB, overlapBases, gapBasesInA, gapBasesInB);
#endif
			// if there is a proper chain between A and B reads, then discard all overlaps between A and B for the repcomp step, (otherwise do repcomp)
			if (properBeg && properEnd && itsBasesInA >= 0 && itsBasesInB >= 0 && gapBasesInA >= 0 && gapBasesInB >= 0
					&& overlapBases * 0.3 > MAX(gapBasesInA, gapBasesInB))
			{
#ifdef DEBUG
				printf("FOUND PROPER CHAIN - EXCLUDE ALL OVLS FROM REPCOMP INPUT\n");
#endif
				for (i = 0; i < novl; i++)
				{
					ovl[i].flags |= OVL_DISCARD;
				}
			}
			else // reset all discard flags that were set during chain detection
			{
				for (i = 0; i < novl; i++)
				{
					if (ovl[i].flags & OVL_DISCARD)
					{
						ovl[i].flags &= ~OVL_DISCARD;
					}
				}
			}
		}
		else // there is only one overlap !!! check this
		{
			int properBeg = (MIN(ovl->path.abpos, ovl->path.bbpos) < 1000);
			int properEnd = (ovl->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread) || ovl->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread));

			if (properBeg && properEnd)
			{
#ifdef DEBUG
				printf("FOUND PROPER SINGLE OVL - EXCLUDE OVL FROM REPCOMP INPUT\n");
#endif
				ovl->flags |= OVL_DISCARD;
			}
		}

		// revert flags
		for (i = 0; i < novl; i++)
		{
			if (ovl[i].flags & OVL_DISCARD)
			{
				ovl[i].flags &= ~OVL_DISCARD;
				ctx->nKeptOvls++;
			}
			else
			{
				ovl[i].flags |= (OVL_DISCARD);
				ctx->nDiscardOvls++;
			}
		}
	}
	else if (ctx->type == 1) // separate for forcealign
	{
		if (ctx->curChains == 1) // ignore everything, where you have more then 1 chain
		{
			Chain *bestChain = ctx->ovlChains;

			int properBeg = 0;
			int properEnd = 0;
			int gapBasesInA = 0;
			int gapBasesInB = 0;
			int itsBasesInA = 0;
			int itsBasesInB = 0;
			int overlapBases = 0;

			// check for proper begin
			if (MIN(bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bbpos) < 1000)
				properBeg = 1;
			// check for proper end
			if (bestChain->ovls[bestChain->novl - 1]->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread)
					|| bestChain->ovls[bestChain->novl - 1]->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread))
				properEnd = 1;

#ifdef DEBUG
			printf("properBeg: %d, properEnd %d\n", properBeg, properEnd);
#endif
			if (properBeg && properEnd)
			{
				int i, its_a, its_b;
				overlapBases = MAX(bestChain->ovls[0]->path.aepos - bestChain->ovls[0]->path.abpos, bestChain->ovls[0]->path.bepos - bestChain->ovls[0]->path.bbpos);

				for (i = 1; i < bestChain->novl; i++)
				{
					its_a = its_b = 0;
					overlapBases += MAX(bestChain->ovls[i]->path.aepos - bestChain->ovls[i]->path.abpos, bestChain->ovls[i]->path.bepos - bestChain->ovls[i]->path.bbpos);
					// check for intersection in A
					if (bestChain->ovls[i]->path.abpos < bestChain->ovls[i - 1]->path.aepos)
					{
						its_a = bestChain->ovls[i - 1]->path.aepos - bestChain->ovls[i]->path.abpos;
						if (its_a > 1000)
						{
							itsBasesInA = -1;
							break;
						}
						itsBasesInA += its_a;
					}
					// check for gap in A
					else
					{
						int gap = bestChain->ovls[i]->path.abpos - bestChain->ovls[i - 1]->path.aepos;
//						if (gap > 1000)
//						{
//							gapBasesInA = -1;
//							break;
//						}
						gapBasesInA += gap;
					}
					// check for intersection in B
					if (bestChain->ovls[i]->path.bbpos < bestChain->ovls[i - 1]->path.bepos)
					{
						its_b = bestChain->ovls[i - 1]->path.bepos - bestChain->ovls[i]->path.bbpos;
						if (its_b > 1000)
						{
							itsBasesInB = -1;
							break;
						}
						itsBasesInB += its_b;
					}
					// check for gap in B
					else
					{
						int gap = bestChain->ovls[i]->path.bbpos - bestChain->ovls[i - 1]->path.bepos;
//						if (gap > 1000)
//						{
//							gapBasesInB = -1;
//							break;
//						}
						gapBasesInB += gap;
					}
					overlapBases -= MAX(its_a, its_b);
				}
			}
#ifdef DEBUG
			printf("if(%d && %d && %d >=0 && %d >=0 && %d >=0 && %d >=0 && %d * 0.3 > MAX(%d, %d))\n", properBeg, properEnd, itsBasesInA, itsBasesInB, gapBasesInA,
					gapBasesInB, overlapBases, gapBasesInA, gapBasesInB);
#endif
			if (properBeg && properEnd && itsBasesInA >= 0 && itsBasesInB >= 0 && gapBasesInA >= 0 && gapBasesInB >= 0
					&& overlapBases * 0.3 > MAX(gapBasesInA, gapBasesInB))
			{
				;
#ifdef DEBUG
				printf("FOUND PROPER CHAIN - INCLUDE ALL CHAIN OVLS FOR FORCEALIGN INPUT\n");
#endif
			}
			else // multiple chains: discard all overlaps, i.e. no forcealign
			{
				for (i = 0; i < novl; i++)
				{
					ovl[i].flags |= OVL_DISCARD;
				}
			}
		}
		else if (ctx->curChains > 1) // multiple chains: discard all overlaps, i.e. no forcealign
		{
			for (i = 0; i < novl; i++)
			{
				ovl[i].flags |= OVL_DISCARD;
			}
		}
		else // there is only one overlap !!! check this
		{
			int properBeg = (MIN(ovl->path.abpos, ovl->path.bbpos) < 1000);
			int properEnd = (ovl->path.aepos + 1000 > DB_READ_LEN(ctx->db, ovl->aread) || ovl->path.bepos + 1000 > DB_READ_LEN(ctx->db, ovl->bread));

			int realBeg = (MIN(ovl->path.abpos, ovl->path.bbpos) == 0);
			int realEnd = (ovl->path.aepos == DB_READ_LEN(ctx->db, ovl->aread) || ovl->path.bepos == DB_READ_LEN(ctx->db, ovl->bread));

			if ((properBeg && properEnd) && !(realBeg && realEnd))
			{
				;
#ifdef DEBUG
				printf("FOUND PROPER SINGLE OVL - INCLUDE OVL FOR FORCEALIGN INPUT\n");
#endif
			}
			else
			{
				ovl->flags |= OVL_DISCARD;
			}
		}

		// revert flags
		for (i = 0; i < novl; i++)
		{
			if (ovl[i].flags & OVL_DISCARD)
			{
				ovl[i].flags &= ~OVL_DISCARD;
				ctx->nKeptOvls++;
			}
			else
			{
				ovl[i].flags |= (OVL_DISCARD);
				ctx->nDiscardOvls++;
			}
		}
	}

	// reset chain and ovl counter
	for (i = 0; i < ctx->curChains; i++)
		ctx->ovlChains[i].novl = 0;
	ctx->curChains = 0;

	return 1;
}

int main(int argc, char* argv[])
{
	HITS_DB db;
	SeparateContext sctx;
	PassContext* pctx;
	FILE* fileOvlIn;
	FILE* fileOvlOut;
	FILE* fileOvlOutDiscard;

	bzero(&sctx, sizeof(SeparateContext));

	sctx.db = &db;

	// args

	char* pcTrackRepeats = DEF_ARG_R;
	int arg_purge = 1;

	sctx.nMinAlnLength = -1;
	sctx.nMinReadLength = -1;
	sctx.nVerbose = 0;
	sctx.useRLoader = 0;

	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "vLo:l:r:T:")) != -1)
	{
		switch (c)
		{

		case 'L':
			sctx.useRLoader = 1;
			break;

		case 'v':
			sctx.nVerbose++;
			break;

		case 'o':
			sctx.nMinAlnLength = atoi(optarg);
			break;

		case 'l':
			sctx.nMinReadLength = atoi(optarg);
			break;

		case 'T':
			sctx.type = atoi(optarg);
			if (sctx.type < 0 || sctx.type > 1)
			{
				fprintf(stderr, "Unsupported type: %d\n", sctx.type);
				exit(1);
			}
			break;

		case 'r':
			pcTrackRepeats = optarg;
			break;

		default:
			fprintf(stderr, "unknown option %c\n", c);
			usage();
			exit(1);
		}
	}

	if (argc - optind != 4)
	{
		usage();
		exit(1);
	}

	char* pcPathReadsIn = argv[optind++];
	char* pcPathOverlapsIn = argv[optind++];
	char* pcPathOverlapsOut = argv[optind++];
	char* pcPathOverlapsOutDiscard = argv[optind++];

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

	if ((fileOvlOutDiscard = fopen(pcPathOverlapsOutDiscard, "w")) == NULL)
	{
		fprintf(stderr, "could not open %s\n", pcPathOverlapsOutDiscard);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "could not open %s\n", pcPathReadsIn);
		exit(1);
	}

	if (pcTrackRepeats)
	{
		sctx.trackRepeat = track_load(&db, pcTrackRepeats);

		if (!sctx.trackRepeat)
		{
			fprintf(stderr, "could not load track %s\n", pcTrackRepeats);
			exit(1);
		}
	}

	// passes

	if (sctx.useRLoader)
	{
		sctx.rl = rl_init(&db, 1);

		pctx = pass_init(fileOvlIn, NULL);

		pctx->data = &sctx;
		pctx->split_b = 1;
		pctx->load_trace = 0;

		pass(pctx, loader_handler);
		rl_load_added(sctx.rl);
		pass_free(pctx);
	}

	pctx = pass_init(fileOvlIn, fileOvlOut);

	pctx->split_b = 1;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &sctx;
	pctx->write_overlaps = 1;
	pctx->purge_discarded = arg_purge;

	separate_pre(pctx, &sctx);
	pass(pctx, separate_handler);
	separate_post(&sctx);

	int ninclude = sctx.nKeptOvls;
	int nexclude = sctx.nDiscardOvls;

	assert(sctx.nKeptOvls + sctx.nDiscardOvls == pctx->novl);

	pass_free(pctx);

	pctx = pass_init(fileOvlIn, fileOvlOutDiscard);

	pctx->split_b = 1;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &sctx;
	pctx->write_overlaps = 1;
	pctx->purge_discarded = arg_purge;

	separate_pre(pctx, &sctx);
	pass(pctx, separateDiscard_handler); // invert discard flag
	separate_post(&sctx);

	assert((sctx.nKeptOvls == nexclude) && (sctx.nDiscardOvls == ninclude));

	pass_free(pctx);

	// cleanup

	if (sctx.useRLoader)
	{
		rl_free(sctx.rl);
	}

	Close_DB(&db);

	fclose(fileOvlOut);
	fclose(fileOvlOutDiscard);
	fclose(fileOvlIn);

	return 0;
}
