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

#define DEF_ARG_S -1
#define DEF_ARG_T TRACK_TRIM
#define DEF_ARG_TT 0
#define DEF_ARG_R TRACK_REPEATS

#define BIN_SIZE 100
#define MIN_LR 500

#define READ_NONE 0x0
#define READ_DISCARD ( 0x1 << 0 )
#define READ_KEEP ( 0x1 << 1 )

#define VERBOSE
#undef VERBOSE_STITCH

#ifdef VERBOSE_STITCH
#define OVL_STRAND( ovl ) ( ( ( ovl )->flags & OVL_COMP ) ? 'c' : 'n' )
#endif

#undef DEBUG_REPEAT_EXTENSION
#undef DEBUG_CHAIN

#define REMOVE_STITCH_OVL ( 1 << 0 )
#define REMOVE_MOD_OVL ( 1 << 1 )
#define REMOVE_TP ( 1 << 2 )
#define REMOVE_NONID_OVL ( 1 << 3 )
#define REMOVE_SPECREP_OVL ( 1 << 4 )
#define REMOVE_ID_OVL ( 1 << 5 )

typedef struct
{
	Overlap **ovls;
	int novl;
	int maxOvl;

} Chain;

typedef struct
{
	// stats counters
	int nFilteredDiffs;
	int nFilteredDiffsSegments;
	int nFilteredUnalignedBases;
	int nFilteredLength;
	int nFilteredRepeat;
	int nFilteredReadLength;
	int nRepeatOvlsKept;
	int nFilteredLocalEnd;
	int nLowCovALn;
	int nSymDiscard;
	int nMultiMapper;
	int nMultiMapperBases;
	int nCovFiltReads;
	int64 nCovFiltOverlaps, nCovFiltBases;

	// settings
	int nStitched;
	float fMaxDiffs;
	int nMaxUnalignedBases, nMinAlnLength;
	int nMinNonRepeatBases, nMinNonRepeatBasesChain, nMinReadLength;
	int nVerbose;
	int stitch;
	int stitch_aggressively;
	int rm_cov;        // repeat modules, coverage
	int rm_aggressive; // -M
	int do_trim;
	int downsample;
	int lowCoverageFilter;
	int hghCoverageFilter;
	char* cov_read_active;
	int remUpToXPercAln;

	int removeFlags; 	// 1 << 0 ... stitched overlaps, 1 << 1 ... module overlaps, 1 << 2 .. trace points, 1 << 3 .. non-identity overlaps,
	// 1 << 4 .. remove B-read repeat overlaps, if a proper overlap between A and B exist i.e A ------------ or A ------------ or A ------------ or A   ------------
	//  																					   B -----E|B				B     E|B-----    B   E|B---E|B		  B ------------------
	// 1 << 5 .. identity overlaps,
	int includeReadFlag;

	int removeLowCoverageOverlaps;
	int removeMultiMappers;
	char* cover_multi_mapper;

	int rp_mergeTips;  // increase repeat interval if it ends with rp_mergeTips Bases according to trim annotation
	// repeat modules - merged repeats
	int rm_merge;
	int rm_mode;

	track_data* rm_repeat;
	unsigned int rm_maxrepeat;

	uint64_t* rm_bins;
	int rm_maxbins;

	// repeat modules - result track

	track_anno* rm_anno;
	track_data* rm_data;
	track_anno rm_ndata;
	track_anno rm_maxdata;

	// local ends
	int* le_lbins;
	int* le_rbins;
	int le_maxbins;

	FILE* fileSpanningReads;

	HITS_DB* db;
	HITS_TRACK* trackRepeat;
	HITS_TRACK* trackTrim;
	HITS_TRACK* trackQ;

	int* r2bin;
	int max_r2bin;

	int useRLoader;
	TRIM* trim;
	Read_Loader* rl;

	ovl_header_twidth twidth;

	FILE* fileOutDiscardedOverlaps;
	int ** discardedAreadList;
	int * discardedBreads;

	int nMaxProperChains;
	Chain *ovlChains;
	int curChains;
	int maxChains;
} FilterContext;

extern char* optarg;
extern int optind, opterr, optopt;

static int loader_handler(void* _ctx, Overlap* ovl, int novl)
{
	FilterContext* ctx = (FilterContext*) _ctx;
	Read_Loader* rl = ctx->rl;

	static int firstCall = 1;

	int i;
	for (i = 0; i < novl; i++)
	{
		int b = ovl[i].bread;

		int trim_b_left, trim_b_right;

		if (ctx->trackTrim)
			get_trim(ctx->db, ctx->trackTrim, b, &trim_b_left, &trim_b_right);
		else
		{
			if (ctx->trackTrim == NULL && firstCall)
			{
				printf("[WARNING] - Read loader is used without trim track. This can cause issues!\n");
				firstCall = 0;
			}

			trim_b_left = 0;
			trim_b_right = DB_READ_LEN(ctx->db, b);
		}

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

static void removeOvls(FilterContext *fctx, Overlap* ovls, int novls, int rmFlag)
{
	int i;

	for (i = 0; i < novls; i++)
	{
		if (rmFlag & REMOVE_TP)
		{
			ovls[i].path.tlen = 0;
			ovls[i].path.trace = NULL;
		}

		if ((rmFlag & REMOVE_MOD_OVL) && (ovls[i].flags & OVL_MODULE))
		{
			ovls[i].flags |= OVL_DISCARD;
		}

		else if ((rmFlag & REMOVE_STITCH_OVL) && (ovls[i].path.tlen == 0))
		{
			ovls[i].flags |= OVL_DISCARD;
		}

		else if ((rmFlag & REMOVE_NONID_OVL) && (ovls[i].aread != ovls[i].bread))
		{
			ovls[i].flags |= OVL_DISCARD;
		}

		else if ((rmFlag & REMOVE_ID_OVL) && (ovls[i].aread == ovls[i].bread))
		{
			ovls[i].flags |= OVL_DISCARD;
		}
	}

	// special case: REMOVE_SPECREP_OVL
	if ((rmFlag & REMOVE_SPECREP_OVL))
	{
		printf("REMOVE_SPECREP_OVL\n");
		int j, k;
		j = k = 0;
		int proper;
		int alen = DB_READ_LEN(fctx->db, ovls->aread);

		while (j < novls)
		{
			int blen = DB_READ_LEN(fctx->db, ovls[j].bread);
			proper = 0;

			if ((ovls[j].path.abpos == 0 || ovls[j].path.bbpos == 0) && (ovls[j].path.aepos == alen || ovls[j].path.bepos == blen))
				proper++;

			while (k < novls - 1 && ovls[j].bread == ovls[k + 1].bread)
			{
				if ((ovls[k + 1].path.abpos == 0 || ovls[k + 1].path.bbpos == 0) && (ovls[k + 1].path.aepos == alen || ovls[k + 1].path.bepos == blen))
					proper++;

				k++;
			}

			printf("aread: %d bread: %d novl: %d proper: %d\n", ovls[j].aread, ovls[j].bread, k - j + 1, proper);

			if (proper == 1 && k - j + 1 > 1)
			{
				int l;
				for (l = 0; l < k - j + 1; l++)
				{
					Overlap *ovl = ovls + j + l;
					if (!((ovl->path.abpos == 0 || ovl->path.bbpos == 0) && (ovl->path.aepos == alen || ovl->path.bepos == blen)))
						ovl->flags |= OVL_DISCARD;
				}
			}
			j = k + 1;
		}
	}
}

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

static void chain(FilterContext *ctx, Overlap *ovls, int n, int trim_ab, int trim_ae)
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

static int stitch(Overlap* ovls, int n, int fuzz, int aggressive)
{
	int stitched = 0;

	if (n < 2)
	{
		return stitched;
	}

	int i, k, b;
	int ab2, ae1, ae2;
	int bb2, be1, be2;

	const int ignore_mask = OVL_CONT | OVL_STITCH | OVL_GAP | OVL_TRIM;

	for (i = 0; i < n; i++)
	{
		Overlap* ovli = ovls + i;

		if (ovli->flags & ignore_mask)
		{
			continue;
		}

		b = ovli->bread;

		ae1 = ovli->path.aepos;
		be1 = ovli->path.bepos;

		int found = 1;

		while (found)
		{
			found = 0;
			int maxk = 0;
			int maxlen = 0;

			for (k = i + 1; k < n && ovls[k].bread <= b; k++)
			{
				Overlap* ovlk = ovls + k;

				if (ovlk->flags & ignore_mask || (ovli->flags & OVL_COMP) != (ovlk->flags & OVL_COMP))
				{
					continue;
				}

				ab2 = ovlk->path.abpos;
				ae2 = ovlk->path.aepos;

				bb2 = ovlk->path.bbpos;
				be2 = ovlk->path.bepos;

				int deltaa = abs(ae1 - ab2);
				int deltab = abs(be1 - bb2);

				if (deltaa < fuzz && deltab < fuzz && (aggressive || abs(deltaa - deltab) < 40))
				{
					if (ae2 - ab2 > maxlen)
					{
						found = 1;

						maxk = k;
						maxlen = ae2 - ab2;
					}
				}
			}

			if (found)
			{
				Overlap* ovlk = ovls + maxk;

				ab2 = ovlk->path.abpos;
				ae2 = ovlk->path.aepos;

				bb2 = ovlk->path.bbpos;
				be2 = ovlk->path.bepos;

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
		}
	}

	return stitched;
}

static int find_repeat_modules(FilterContext* ctx, Overlap* ovls, int novl)
{
	int a = ovls->aread;
	int trim_ab, trim_ae;
	const int exclude_mask = OVL_LOCAL | OVL_TRIM | OVL_CONT | OVL_STITCH | OVL_GAP | OVL_DIFF;
	const int allowed_mask = OVL_DISCARD | OVL_REPEAT | OVL_COMP; // | OVL_OLEN | OVL_RLEN;

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, a, &trim_ab, &trim_ae);
	}
	else
	{
		trim_ab = 0;
		trim_ae = DB_READ_LEN(ctx->db, a);
	}

	if (trim_ab >= trim_ae)
	{
		return 1;
	}

	int i;
	uint32_t left = 0;
	uint32_t right = 0;

	int left_potential = 0;
	int right_potential = 0;

	// if there are non-repeat overlaps, entering left and leaving right, then no module overlaps are needed

	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;
		int abpos = ovl->path.abpos;
		int aepos = ovl->path.aepos;
		int flags = ovl->flags;

		// contained

		if (abpos == trim_ab && aepos == trim_ae)
		{
			left = 1;
			right = 1;
			break;
		}

		// potential exits, only allowed to have discard, repeat or comp flag

		if ((flags & (~allowed_mask)) == 0)
		{
			if (abpos == trim_ab)
			{
				left_potential += 1;
			}

			if (aepos == trim_ae)
			{
				right_potential += 1;
			}
		}

		// exits left / right

		if (flags & OVL_DISCARD)
		{
			continue;
		}

		if (abpos == trim_ab)
		{
			left += 1;
		}

		if (aepos == trim_ae)
		{
			right += 1;
		}
	}

	if ((left > 0 && right > 0) || (left == 0 && right == 0))
	{
		return 1;
	}

	track_anno* ranno = (track_anno*) (ctx->trackRepeat->anno);
	track_data* rdata = (track_data*) (ctx->trackRepeat->data);

	track_anno ob = ranno[a];
	track_anno oe = ranno[a + 1];

	// no repeats, nothing to do

	if (ob >= oe)
	{
		return 1;
	}

	if (oe - ob > ctx->rm_maxrepeat) // bytes
	{
		ctx->rm_maxrepeat = (oe - ob) + 128;
		ctx->rm_repeat = malloc(ctx->rm_maxrepeat);
	}

	// merge repeats close to each other

	int nrepeat = 0;

	ob /= sizeof(track_data);
	oe /= sizeof(track_data);

	int rm_merge = ctx->rm_merge;
	if (rm_merge < ctx->nMinNonRepeatBases)
		rm_merge = ctx->nMinNonRepeatBases;

	while (ob < oe)
	{
		int b = MAX(trim_ab, rdata[ob]);
		int e = MIN(trim_ae, rdata[ob + 1]);
		ob += 2;

		if (b >= e)
		{
			continue;
		}
		/// try to ignore to add up low complexity intervals from DBdust
		if ((nrepeat > 0) && (b - ctx->rm_repeat[nrepeat - 1] < rm_merge) && (rm_merge < 100 || (e - b > 100)))
		{
			// ctx->stats_merged++;
			ctx->rm_repeat[nrepeat - 1] = e;
		}
		else
		{
			ctx->rm_repeat[nrepeat++] = b;
			ctx->rm_repeat[nrepeat++] = e;
		}
	}

	// for each segment count number of reads anchored with at least MIN_LR

	bzero(ctx->rm_bins, sizeof(uint64_t) * ctx->rm_maxbins);

	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;

		if (ovl->flags & ( OVL_STITCH | OVL_TRIM))
		{
			continue;
		}

		int b = (ovl->path.abpos + MIN_LR) / BIN_SIZE;
		int e = (ovl->path.aepos - MIN_LR) / BIN_SIZE;

		// spanning local alignments indicate non-reliable points inside the repeat

		int incr = 1;

		if (ctx->rm_aggressive == 0 && (ovl->flags & OVL_LOCAL))
		{
			incr = ctx->rm_cov;
		}

		while (b < e)
		{
			ctx->rm_bins[b] += incr;
			b++;
		}
	}

	track_anno prev_offset = ctx->rm_ndata;

	for (i = 0; i < nrepeat; i += 2)
	{
		int rb = ctx->rm_repeat[i];
		int re = ctx->rm_repeat[i + 1];

		// skip repeats if there are valid spanning overlaps

		int skip = 0;
		int j;

		for (j = 0; j < novl; j++)
		{
			Overlap* ovl = ovls + j;

			if ((ovl->flags & OVL_DISCARD))
			{
				continue;
			}

			if (ovl->path.abpos + MIN_LR <= rb && ovl->path.aepos - MIN_LR >= re)
			{
				skip = 1;
				break;
			}
		}

		if (skip)
		{
			continue;
		}

		// gaps inside the repeat with < expected coverage are potential repeat modules

		int b = MAX(trim_ab + 1000, rb + MIN_LR) / BIN_SIZE;
		int e = MIN(trim_ae - 1000, re - MIN_LR) / BIN_SIZE;

		int beg = -1;

		while (b < e)
		{
			if (ctx->rm_bins[b] > 0 && ctx->rm_bins[b] < (uint64_t) ctx->rm_cov)
			{
				//printf("READ %7d POINT @ %5d..%5d %3d %2llu\n", a, b * BIN_SIZE - 50, b * BIN_SIZE + 50, b, ctx->rm_bins[b]);

				if (beg == -1)
				{
					beg = b;
				}
			}
			else
			{
				if (beg != -1)
				{
					if (b - beg > 7)
					{
						//printf("MOD  %7d POINT @ %5d %5d %5d..%5d\n", a, (b + beg) / 2, (b + beg) * BIN_SIZE / 2, beg, b);

						if (ctx->rm_ndata + 2 > ctx->rm_maxdata)
						{
							ctx->rm_maxdata = 1.2 * ctx->rm_ndata + 100;
							ctx->rm_data = realloc(ctx->rm_data, ctx->rm_maxdata * sizeof(track_data));
						}

						int intb = (b + beg) / 2 * BIN_SIZE - 50;
						int inte = intb + 100;

						ctx->rm_anno[a] += 2 * sizeof(track_data);

						ctx->rm_data[ctx->rm_ndata++] = intb;
						ctx->rm_data[ctx->rm_ndata++] = inte;
					}

					beg = -1;
				}
			}

			b++;
		}
	}

	// restore discarded overlaps spanning the repeat module junction

	uint32_t enabled = 0;

	while (prev_offset < ctx->rm_ndata)
	{
		int b = ctx->rm_data[prev_offset++];
		int e = ctx->rm_data[prev_offset++];

		// ignore module that would result in excessive coverage at the junction

		//printf("MODULE @ %d..%d\n", b, e);
		int cov = 0;

		for (i = 0; i < novl; i++)
		{
			Overlap* ovl = ovls + i;

			if (ovl->path.abpos + 100 < b && ovl->path.aepos - 100 > e)
			{
				if (!(ovl->flags & OVL_DISCARD) || ((ovl->flags & OVL_REPEAT) && !(ovl->flags & exclude_mask)))
				{
					cov++;
				}
			}
		}

		if (cov > ctx->rm_cov || cov < ctx->rm_cov / 2)
		{
			continue;
		}

		for (i = 0; i < novl; i++)
		{
			Overlap* ovl = ovls + i;

			if (!(ovl->flags & OVL_REPEAT) || (ovl->flags & exclude_mask))
			{
				continue;
			}

			// MARTIN ignore overlaps that do not have an overhang
			if (ovl->path.abpos > trim_ab && ovl->path.aepos < trim_ae)
			{
				continue;
			}

			if (ovl->path.abpos + MIN_LR < b && ovl->path.aepos - MIN_LR > e)
			{
				ovl->flags &= ~OVL_DISCARD;
				ovl->flags |= OVL_MODULE;

				ctx->nRepeatOvlsKept++;
				enabled += 1;
			}
		}
	}

	if (enabled)
	{
		return 1;
	}

	printf("%d | left = %d right %d lp %d rp %d\n", a, left, right, left_potential, right_potential);

	/*

	 // try to relax the min non-repeat bases condition

	 ob = ranno[a] / sizeof(track_data);
	 oe = ranno[a + 1] / sizeof(track_data);

	 int dist_left = INT_MAX;
	 int dist_right = INT_MAX;

	 while ( ob < oe )
	 {
	 int rb = MAX(trim_ab, rdata[ob]);
	 int re = MIN(trim_ae, rdata[ob + 1]);
	 ob += 2;

	 if (rb >= re)
	 {
	 continue;
	 }

	 dist_left = MIN( rb - trim_ab, dist_left );
	 dist_right = MIN( trim_ae - re, dist_right );
	 }

	 printf("%d | min dist left %d right %d\n", a, dist_left, dist_right);

	 if ( (left == 0 && left_potential <= ctx->rm_cov && dist_left > 100 && dist_left != INT_MAX) ||
	 (right == 0 && right_potential <= ctx->rm_cov && dist_right > 100 && dist_right != INT_MAX) )
	 {
	 for ( i = 0 ; i < novl ; i++ )
	 {
	 Overlap* ovl = ovls + i;

	 if ( (left == 0 && ovl->path.abpos == trim_ab && !(ovl->flags & (~allowed_mask)) ) ||
	 (right == 0 && ovl->path.aepos == trim_ae && !(ovl->flags & (~allowed_mask)) ) )
	 {
	 ovl->flags &= ~OVL_DISCARD;
	 ovl->flags |= OVL_MODULE;
	 ctx->nRepeatOvlsKept++;
	 enabled += 1;
	 }
	 }
	 }

	 printf("%d | enabled %d\n", a, enabled);

	 if (enabled)
	 {
	 return 1;
	 }
	 */

	//
	// TODO --- coverage ... don't count multiple overlaps with the same b ????
	//
	if (ctx->rm_mode < 2)
	{
		return 1;
	}

	// look at b reads and see if they would lead is to a unique region

	if (novl > ctx->max_r2bin)
	{
		ctx->max_r2bin = novl * 1.2 + 128;
		ctx->r2bin = realloc(ctx->r2bin, sizeof(int) * ctx->max_r2bin);
	}

	int potential = 0;

	bzero(ctx->r2bin, sizeof(int) * ctx->max_r2bin);
	bzero(ctx->rm_bins, sizeof(uint64_t) * ctx->rm_maxbins);

	int binsize = MAX(BIN_SIZE, 1000);

	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;
		int trim_bb, trim_be;
		int b = ovl->bread;

		if ((ovl->flags & exclude_mask) || (left == 0 && ovl->path.abpos != trim_ab) || (right == 0 && ovl->path.aepos != trim_ae))
		{
			continue;
		}

		get_trim(ctx->db, ctx->trackTrim, b, &trim_bb, &trim_be);

		printf("%7d | btrim  %5d..%5d\n", b, trim_bb, trim_be);

		int bb, be;

		if (ovl->flags & OVL_COMP)
		{
			int blen = DB_READ_LEN(ctx->db, b);
			bb = blen - ovl->path.bepos;
			be = blen - ovl->path.bbpos;
		}
		else
		{
			bb = ovl->path.bbpos;
			be = ovl->path.bepos;
		}

		ob = ranno[b] / sizeof(track_data);
		oe = ranno[b + 1] / sizeof(track_data);
		int rb, re;

		int nrb = 0;
		int prev_end = trim_bb;
		int rlen_in_b = 0;

		while (ob < oe)
		{
			rb = MAX(trim_bb, rdata[ob]);
			re = MIN(trim_be, rdata[ob + 1]);
			ob += 2;

			printf("%d | %7d | b repeat %5d..%5d\n", a, b, rb, re);

			if (rb >= re)
			{
				continue;
			}

			if (left == 0)
			{
				if (ovl->flags & OVL_COMP)
				{
					if (rb < be && be < re)
					{
						rlen_in_b = re - be;
					}
				}
				else
				{
					if (rb < bb && bb < re)
					{
						rlen_in_b = bb - rb;
					}
				}
			}
			else
			{
				if (ovl->flags & OVL_COMP)
				{
					if (rb < bb && bb < re)
					{
						rlen_in_b = bb - rb;
					}
				}
				else
				{
					if (rb < be && be < re)
					{
						rlen_in_b = re - be;
					}
				}
			}

			nrb += rb - prev_end;
			prev_end = re;
		}

		nrb += trim_be - prev_end;

		if (nrb < ctx->nMinNonRepeatBases || rlen_in_b == 0)
		{
			continue;
		}

		printf("%d -> %7d | leftover rlen %5d nrb %5d", a, b, rlen_in_b, nrb);

		if (!(ovl->flags & (~allowed_mask)))
		{
			ovl->flags |= OVL_TEMP;
			potential += 1;

			printf("  YES");

			int bin = rlen_in_b / binsize;

			ctx->rm_bins[bin] += 1;
			ctx->r2bin[i] = bin;
		}

		printf("\n");
	}

	if (potential > 0)
	{
		for (i = 0; i < novl; i++)
		{
			Overlap* ovl = ovls + i;

			if (ovl->flags & OVL_TEMP)
			{
				int bin = ctx->r2bin[i];
				ovl->flags &= ~OVL_TEMP;

				if (ctx->rm_bins[bin] > 2 && ctx->rm_bins[bin] < (uint32_t) ctx->rm_cov)
				{
					ovl->flags &= ~OVL_DISCARD;
					ovl->flags |= OVL_OPTIONAL;
					ctx->nRepeatOvlsKept++;
					enabled += 1;
				}
			}
		}

		if (enabled)
		{
			return 1;
		}
	}

	if (ctx->rm_mode < 3)
	{
		return 1;
	}

	int prevb = -1;
	int distinctb = 1;

	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;
		int b = ovl->bread;

		if ((ovl->flags & exclude_mask) || (left == 0 && ovl->path.abpos != trim_ab) || (right == 0 && ovl->path.aepos != trim_ae))
		{
			continue;
		}

		if (prevb != b)
		{
			distinctb += 1;
			prevb = b;
		}
	}

	if (distinctb < ctx->rm_cov)
	{
		prevb = -1;
		int maxbidx = -1;
		int maxblen = -1;

		for (i = 0; i < novl; i++)
		{
			Overlap* ovl = ovls + i;
			int b = ovl->bread;
			int len = ovl->path.aepos - ovl->path.abpos;

			if ((ovl->flags & exclude_mask) || (left == 0 && ovl->path.abpos != trim_ab) || (right == 0 && ovl->path.aepos != trim_ae))
			{
				continue;
			}

			if (prevb != b)
			{
				if (maxbidx != -1)
				{
					ovls[maxbidx].flags &= ~OVL_DISCARD;
					ovls[maxbidx].flags |= OVL_OPTIONAL;
					ctx->nRepeatOvlsKept++;
					enabled += 1;
				}

				prevb = b;
				maxbidx = i;
				maxblen = len;
			}
			else
			{
				if (len > maxblen)
				{
					maxbidx = i;
					maxblen = len;
				}
			}
		}

		if (maxbidx != -1)
		{
			ovls[maxbidx].flags &= ~OVL_DISCARD;
			ovls[maxbidx].flags |= OVL_OPTIONAL;
			ctx->nRepeatOvlsKept++;
			enabled += 1;
		}
	}

	return 1;
}

static int filter(FilterContext* ctx, Overlap* ovl)
{
	int nLen = ovl->path.aepos - ovl->path.abpos;
	int nLenB = ovl->path.bepos - ovl->path.bbpos;
	int ret = 0;

	if (nLenB < nLen)
	{
		nLen = nLenB;
	}

	int trim_ab, trim_ae, trim_bb, trim_be;
	int trim_alen, trim_blen;

	int ovlALen = DB_READ_LEN(ctx->db, ovl->aread);
	int ovlBLen = DB_READ_LEN(ctx->db, ovl->bread);

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trim_ab, &trim_ae);
		trim_alen = trim_ae - trim_ab;

		get_trim(ctx->db, ctx->trackTrim, ovl->bread, &trim_bb, &trim_be);
		trim_blen = trim_be - trim_bb;

		if (ovl->flags & OVL_COMP)
		{
			int t = trim_bb;
			trim_bb = ovlBLen - trim_be;
			trim_be = ovlBLen - t;
		}
	}
	else
	{
		trim_ab = 0;
		trim_ae = ovlALen;

		trim_bb = 0;
		trim_be = ovlBLen;

		trim_alen = ovlALen;
		trim_blen = ovlBLen;
	}

	// check list of discarded overlaps
	if (ctx->discardedAreadList && ctx->discardedBreads)
	{
		int a, b;
		if (ovl->aread < ovl->bread)
		{
			a = ovl->aread;
			b = ovl->bread;
		}
		else
		{
			a = ovl->bread;
			b = ovl->aread;
		}

		if (ctx->discardedAreadList[a])
		{
			int discard = 0;
			int j;
			for (j = 0; ctx->discardedAreadList[a][j] > -1; j++)
			{
				if (ctx->discardedAreadList[a][j] == b)
				{
					discard = 1;
					break;
				}
			}
			if (discard)
			{
				ret |= OVL_DISCARD | OVL_SYMDISCARD;
				ctx->nSymDiscard++;

				if (ctx->nVerbose)
				{
					printf("overlap (%d x %d): in list of discarded overlaps\n", ovl->aread, ovl->bread);
				}

			}
		}
	}

	if (ctx->nMinReadLength != -1 && (trim_alen < ctx->nMinReadLength || trim_blen < ctx->nMinReadLength))
	{
		ctx->nFilteredReadLength++;
		ret |= OVL_DISCARD | OVL_RLEN;
	}

	if (ctx->nMinAlnLength != -1 && nLen < ctx->nMinAlnLength)
	{
		if (ctx->nVerbose)
		{
			printf("overlap %d -> %d: drop due to length %d\n", ovl->aread, ovl->bread, nLen);
		}

		ctx->nFilteredLength++;
		ret |= OVL_DISCARD | OVL_OLEN;
	}

	if (ctx->nMaxUnalignedBases != -1)
	{
		if (((ovl->path.abpos - trim_ab) > ctx->nMaxUnalignedBases && (ovl->path.bbpos - trim_bb) > ctx->nMaxUnalignedBases)
				|| ((trim_ae - ovl->path.aepos) > ctx->nMaxUnalignedBases && (trim_be - ovl->path.bepos) > ctx->nMaxUnalignedBases))
		{
			if (ctx->nVerbose)
			{
				printf("overlap %d -> %d: drop due to unaligned overhang [%d, %d -> trim %d, %d], [%d, %d -> trim %d, %d]\n", ovl->aread, ovl->bread, ovl->path.abpos,
						ovl->path.aepos, trim_ab, trim_ae, ovl->path.bbpos, ovl->path.bepos, trim_bb, trim_be);
			}

			ctx->nFilteredUnalignedBases++;

			ret |= OVL_DISCARD | OVL_LOCAL;
		}
	}

	if (ctx->fMaxDiffs > 0 && ctx->remUpToXPercAln == 0)
	{
		if (1.0 * ovl->path.diffs / nLen > ctx->fMaxDiffs)
		{
			if (ctx->nVerbose)
			{
				printf("overlap %d -> %d: drop due to diffs %d length %d\n", ovl->aread, ovl->bread, ovl->path.diffs, nLen);
			}

			ctx->nFilteredDiffs++;

			ret |= OVL_DISCARD | OVL_DIFF;

			if (ctx->fileOutDiscardedOverlaps)
			{
				fprintf(ctx->fileOutDiscardedOverlaps, "%d %d\n", ovl->aread, ovl->bread);
			}

		}
	}

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

				// ignore repeat in front of trim intervaL
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

				// ignore repeat behind end of trim intervaL
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
			ret |= OVL_DISCARD | OVL_REPEAT;
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
				ret |= OVL_DISCARD | OVL_REPEAT;
			}
		}

		if (!(ret & OVL_REPEAT))
		{
			while (b < e)
			{
				rb = repeats_data[b];
				re = repeats_data[b + 1];

				if (re < rp_mergeTip_ab || rb > rp_mergeTip_ae)
				{
					b += 2;
					continue;
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
				ret |= OVL_DISCARD | OVL_REPEAT;
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
				ret |= OVL_DISCARD | OVL_REPEAT;
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

				ret |= OVL_DISCARD | OVL_REPEAT;
			}

		}
	}

	// check tracepoints
	if (!(ovl->flags & OVL_DISCARD) && (ctx->removeFlags & REMOVE_TP))
	{
		ovl_trace* trace = (ovl_trace*) ovl->path.trace;

		int bpos = ovl->path.bbpos;

		int j;

		for (j = 0; j < ovl->path.tlen; j += 2)
		{
			bpos += trace[j + 1];
		}

		if (bpos != ovl->path.bepos)
		{
			ret |= OVL_DISCARD;

			if (ctx->nVerbose)
			{
				printf("overlap (%d x %d): pass-through points inconsistent be = %d (expected %d)\n", ovl->aread, ovl->bread, bpos, ovl->path.bepos);
			}
		}
	}

	return ret;
}

static void write_spanning_reads(FilterContext* ctx, int aread)
{
	int b, e, rb, re;
	int trim_ab, trim_ae;

	track_anno* repeats_anno = ctx->trackRepeat->anno;
	track_data* repeats_data = ctx->trackRepeat->data;

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, aread, &trim_ab, &trim_ae);
	}
	else
	{
		trim_ab = 0;
		trim_ae = DB_READ_LEN(ctx->db, aread);
	}

	b = repeats_anno[aread] / sizeof(track_data);
	e = repeats_anno[aread + 1] / sizeof(track_data);

	int minbases = MAX(ctx->nMinNonRepeatBases, 300);

	while (b < e)
	{
		rb = repeats_data[b];
		re = repeats_data[b + 1];

		if (rb - minbases > trim_ab && re + minbases < trim_ae)
		{
			fprintf(ctx->fileSpanningReads, "%d %d\n", aread, re - rb);
		}

		b += 2;
	}
}

static void filter_pre(PassContext* pctx, FilterContext* fctx)
{
#ifdef VERBOSE
	printf( ANSI_COLOR_GREEN "PASS filtering\n" ANSI_COLOR_RESET);
#endif

	fctx->twidth = pctx->twidth;

	// trim

	if (fctx->do_trim)
	{
		fctx->trim = trim_init(fctx->db, pctx->twidth, fctx->trackTrim, fctx->rl);
	}

	// repeat modules

	fctx->rm_anno = (track_anno*) malloc(sizeof(track_anno) * ( DB_NREADS( fctx->db ) + 1));
	bzero(fctx->rm_anno, sizeof(track_anno) * ( DB_NREADS( fctx->db ) + 1));

	fctx->rm_ndata = 0;
	fctx->rm_maxdata = 100;
	fctx->rm_data = (track_data*) malloc(sizeof(track_data) * fctx->rm_maxdata);

	fctx->rm_maxbins = ( DB_READ_MAXLEN( fctx->db ) + BIN_SIZE) / BIN_SIZE;
	fctx->rm_bins = malloc(sizeof(uint64_t) * fctx->rm_maxbins);

	fctx->le_maxbins = ( DB_READ_MAXLEN( fctx->db ) + BIN_SIZE) / BIN_SIZE;
	fctx->le_lbins = malloc(sizeof(int) * fctx->le_maxbins);
	fctx->le_rbins = malloc(sizeof(int) * fctx->le_maxbins);

	if (fctx->removeMultiMappers > 1)
	{
		fctx->cover_multi_mapper = (char*) malloc(DB_READ_MAXLEN(fctx->db));
	}

	if (fctx->nMaxProperChains >= 0)
	{
		fctx->curChains = 0;
		fctx->maxChains = 5;
		fctx->ovlChains = (Chain*) malloc(sizeof(Chain) * MAX(fctx->maxChains, fctx->nMaxProperChains));
		bzero(fctx->ovlChains, sizeof(Chain) * MAX(fctx->maxChains, fctx->nMaxProperChains));
	}
}

static void filter_post(FilterContext* ctx)
{
#ifdef VERBOSE
	if (ctx->trim)
	{
		printf("trimmed %'lld of %'lld overlaps\n", ctx->trim->nTrimmedOvls, ctx->trim->nOvls);
		printf("trimmed %'lld of %'lld bases\n", ctx->trim->nTrimmedBases, ctx->trim->nOvlBases);
	}

	if (ctx->nFilteredReadLength > 0)
	{
		printf("min read length of %5d discarded         %10d\n", ctx->nMinReadLength, ctx->nFilteredReadLength);
	}

	if (ctx->nFilteredRepeat > 0)
	{
		printf("min non-repeat bases of %4d discarded     %10d\n", ctx->nMinNonRepeatBases, ctx->nFilteredRepeat);
	}

	if (ctx->nRepeatOvlsKept > 0)
	{
		printf("  kept %d repeat overlaps (mode %d)\n", ctx->nRepeatOvlsKept, ctx->rm_mode);
	}

	if (ctx->nFilteredDiffs > 0)
	{
		printf("diff threshold of %.1f discarded             %d\n", ctx->fMaxDiffs * 100., ctx->nFilteredDiffs);
	}

	if (ctx->nFilteredDiffsSegments > 0)
	{
		printf("diff threshold of %3.1f on segments discarded  %10d\n", ctx->fMaxDiffs * 100., ctx->nFilteredDiffsSegments);
	}

	if (ctx->nFilteredLength > 0)
	{
		printf("min overlap length of %d discarded                %10d\n", ctx->nMinAlnLength, ctx->nFilteredLength);
	}

	if (ctx->nFilteredUnalignedBases > 0)
	{
		printf("unaligned bases threshold of %4d discarded    %10d\n", ctx->nMaxUnalignedBases, ctx->nFilteredUnalignedBases);
	}

	if (ctx->nMultiMapper > 0)
	{
		printf("multiMapper of %4d discarded    %10d\n", ctx->nMultiMapper, ctx->nMultiMapperBases);
	}

	if (ctx->nStitched > 0)
	{
		printf("stitched %4d at fuzzing of                   %10d\n", ctx->nStitched, ctx->stitch);
	}

	if (ctx->nFilteredLocalEnd)
	{
		printf("local ends discarded                              %10d\n", ctx->nFilteredLocalEnd);
	}

	if (ctx->nLowCovALn)
	{
		printf("low tip coverage of <%2d discarded                %10d\n", ctx->removeLowCoverageOverlaps, ctx->nLowCovALn);
	}

	if (ctx->nSymDiscard)
	{
		printf("symmetrically remove discarded overlaps          %10d\n", ctx->nSymDiscard);
	}

	if (ctx->nCovFiltReads)
	{
		printf("coverage filtered reads %4d -> overlaps %10lld -> bases %10lld\n", ctx->nCovFiltReads, ctx->nCovFiltOverlaps, ctx->nCovFiltBases);
	}

#endif

	if (ctx->trim)
	{
		trim_close(ctx->trim);
		free(ctx->trim);
	}

	free(ctx->le_lbins);
	free(ctx->le_rbins);

	free(ctx->rm_bins);
	free(ctx->rm_anno);
	free(ctx->rm_data);

	if (ctx->discardedBreads)
		free(ctx->discardedBreads);
	if (ctx->discardedAreadList)
		free(ctx->discardedAreadList);

	if (ctx->removeMultiMappers > 1)
		free(ctx->cover_multi_mapper);

	if (ctx->nMaxProperChains >= 0)
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

static int cmp_ovls_qual(const void* a, const void* b)
{
	Overlap* o1 = *(Overlap**) a;
	Overlap* o2 = *(Overlap**) b;

	return ((100 - (o1->path.diffs * 100.0 / (o1->path.aepos - o1->path.abpos))) * 10 - (100 - (o2->path.diffs * 100.0 / (o2->path.aepos - o2->path.abpos))) * 10);
}

static void removeWorstAlignments(FilterContext* ctx, Overlap* ovl, int novl)
{
	int i;
	Overlap** ovl_sort = (Overlap**) malloc(sizeof(Overlap*) * novl);

	int numIncomingReads, numLeavingReads;
	int cumOverallBases;

	numIncomingReads = numLeavingReads = cumOverallBases = 0;

	int trimABeg, trimAEnd;

	trimABeg = 0;
	trimAEnd = DB_READ_LEN(ctx->db, ovl->aread);

	if (ctx->trackTrim)
		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trimABeg, &trimAEnd);

	int aTrimLen = trimAEnd - trimABeg;
	int aQuarterTrimLen = aTrimLen / 4;

	for (i = 0; i < novl; i++)
	{
		ovl_sort[i] = ovl + i;

		if (ovl_sort[i]->flags & OVL_DISCARD)
			continue;

		if (ovl_sort[i]->path.abpos <= trimABeg)
			numIncomingReads += 1;

		if (ovl_sort[i]->path.aepos >= trimAEnd)
			numLeavingReads += 1;

		cumOverallBases += ovl_sort[i]->path.aepos - ovl_sort[i]->path.abpos;
	}

	printf("Coverage[%d]: beg,end [%3d, %3d] avgCov %.2f\n", ovl->aread, numIncomingReads, numLeavingReads,
			cumOverallBases*1.0/aTrimLen);

	qsort(ovl_sort, novl, sizeof(Overlap*), cmp_ovls_qual);

	// todo hard coded
	int MinTipCov = MIN(3, ctx->removeLowCoverageOverlaps);
	int MaxRemovedAlnBasesPerc = ctx->remUpToXPercAln;
	float maxQV = (ctx->fMaxDiffs > 0) ? ctx->fMaxDiffs : 28;
	int removeAlnBases = 0;

	int numRemovedIncomingReads, numRemovedLeavingReads;
	numRemovedIncomingReads = numRemovedLeavingReads = 0;

	for (i = 0; i < novl; i++)
	{
		Overlap* so = ovl_sort[i];
		if (so->flags & OVL_DISCARD)
			continue;

		int err = (int) (so->path.diffs * 100.0 / (so->path.aepos - so->path.abpos));

		if(so->path.abpos <= trimABeg)
			numRemovedIncomingReads++;

		if(so->path.aepos >= trimAEnd)
			numRemovedLeavingReads++;

		removeAlnBases +=  so->path.aepos - so->path.abpos;

		if (removeAlnBases*100.0/cumOverallBases < MaxRemovedAlnBasesPerc && numIncomingReads - numRemovedIncomingReads > MinTipCov
				&& numLeavingReads - numRemovedLeavingReads > MinTipCov)
		{
				so->flags |= OVL_DISCARD | OVL_DIFF;
				ctx->nFilteredDiffs += 1;
				if(ctx->nVerbose)
				printf("DISCARD %d%% bad overlaps [%d, %d] a[%d, %d] b [%d, %d] %c ODIF %d\n", MaxRemovedAlnBasesPerc, so->aread, so->bread,
						so->path.abpos, so->path.aepos, so->path.bbpos, so->path.bepos, (so->flags & OVL_COMP) ? 'C' : 'N', err);
		}
		else
		{
			if(ctx->nVerbose)
				printf("DO NOT DISCARD %d%% bad overlaps [%d, %d] a[%d, %d] b [%d, %d] %c ODIF %d\n", MaxRemovedAlnBasesPerc, so->aread, so->bread,
					so->path.abpos, so->path.aepos, so->path.bbpos, so->path.bepos, (so->flags & OVL_COMP) ? 'C' : 'N', err);
			break;
		}

		if(err < maxQV)
		{
			if(ctx->nVerbose)
				printf("STOP reached maxqv of %d [%d, %d] a[%d, %d] b [%d, %d] %c\n",err, so->aread, so->bread,
					so->path.abpos, so->path.aepos, so->path.bbpos, so->path.bepos, (so->flags & OVL_COMP) ? 'C' : 'N');
			break;
		}
	}

	free(ovl_sort);
}

static void filterByCoverage(FilterContext* ctx, Overlap* ovl, int novl)
{

	int j, k, l;
	j = k = 0;

	int64 cov;
	int64 bases = 0;

	bzero(ctx->cov_read_active, DB_READ_MAXLEN(ctx->db));

	int trimABeg, trimAEnd;
	int trimBBeg, trimBEnd;

	trimABeg = 0;
	trimAEnd = DB_READ_LEN(ctx->db, ovl->aread);

	if (ctx->trackTrim)
		get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trimABeg, &trimAEnd);

	while (j < novl)
	{
		while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
		{
			k++;
		}

		if (ovl[j].aread == ovl[j].bread)
		{
			j = k + 1;
			continue;
		}

		trimBBeg = 0;
		trimBEnd = DB_READ_LEN(ctx->db, ovl[j].bread);

		if (ctx->trackTrim)
			get_trim(ctx->db, ctx->trackTrim, ovl[j].bread, &trimBBeg, &trimBEnd);

		chain(ctx, ovl + j, k - j + 1, trimABeg, trimAEnd);

		if (ctx->curChains > 0)
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
			if (MIN(bestChain->ovls[0]->path.abpos - trimABeg, bestChain->ovls[0]->path.bbpos) - trimBBeg < 1000)
				properBeg = 1;
			// check for proper end
			if (MIN(trimAEnd - bestChain->ovls[bestChain->novl - 1]->path.aepos, trimBEnd - bestChain->ovls[bestChain->novl - 1]->path.bepos) < 1000)
				properEnd = 1;

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
						if (gap > 1000)
						{
							gapBasesInA = -1;
							break;
						}
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
						if (gap > 1000)
						{
							gapBasesInB = -1;
							break;
						}
						gapBasesInB += gap;
					}
					overlapBases -= MAX(its_a, its_b);
				}
			}

			// if there is a proper chain between A and B reads, then discard all overlaps between A and B for the repcomp step, (otherwise do repcomp)
			if (properBeg && properEnd && itsBasesInA >= 0 && itsBasesInB >= 0 && gapBasesInA >= 0 && gapBasesInB >= 0
					&& overlapBases * 0.3 > MAX(gapBasesInA, gapBasesInB) && overlapBases >= ctx->nMinAlnLength)
			{
				for (l = 0; l < bestChain->novl; l++)
				{
					bases += bestChain->ovls[l]->path.aepos - bestChain->ovls[l]->path.abpos;

					memset(ctx->cov_read_active + bestChain->ovls[l]->path.abpos, 1, bestChain->ovls[l]->path.aepos - bestChain->ovls[l]->path.abpos);
				}
			}
			// reset chain and ovl counter
			for (l = 0; l < ctx->curChains; l++)
			{
				int m;
				// reset OVL_TEMP flag ovl
				for (m = j; m < k - j + 1; m++)
				{
					if (ovl[m].flags & OVL_TEMP)
						ovl[m].flags &= ~OVL_TEMP;
				}
				ctx->ovlChains[l].novl = 0;
			}
			ctx->curChains = 0;
		}
		else
		{
			int properBeg = 0;
			int properEnd = 0;

			// check for proper begin
			if (MIN(ovl->path.abpos - trimABeg, ovl->path.bbpos - trimBBeg) < 1000)
				properBeg = 1;

			// check for proper end
			if (MIN(trimAEnd - ovl->path.aepos, trimBEnd - ovl->path.bepos) < 1000)
				properEnd = 1;

			if (properBeg && properEnd && (ovl->path.aepos - ovl->path.abpos) > ctx->nMinAlnLength)
			{
				bases += ovl->path.aepos - ovl->path.abpos;

				memset(ctx->cov_read_active + ovl->path.abpos, 1, ovl->path.aepos - ovl->path.abpos);
			}
		}
		j = k + 1;
	}

	int active = 0;
	for (j = trimABeg; j < trimAEnd; j++)
	{
		active += ctx->cov_read_active[j];
	}

	if (active > 0)
	{
		cov = bases / active;
	}
	else
	{
		cov = 0;
	}

	if (cov >= ctx->lowCoverageFilter && cov <= ctx->hghCoverageFilter)
	{
		for (j = 0; j < novl; j++)
		{
			ovl[j].flags |= OVL_DISCARD;
			ctx->nCovFiltBases += ovl[j].path.aepos - ovl[j].path.abpos;
		}
		ctx->nCovFiltReads++;
		ctx->nCovFiltOverlaps += novl;
		printf("COV FILT OUT READ %d (novl: %d) cov: %lld range [%d, %d] active: %d\n", ovl->aread, novl, cov, ctx->lowCoverageFilter, ctx->hghCoverageFilter,
				active);

	}
	else
		printf("COV READ %d (novl: %d) cov: %lld range [%d, %d] active: %d\n", ovl->aread, novl, cov, ctx->lowCoverageFilter, ctx->hghCoverageFilter, active);
}

static int filter_handler(void* _ctx, Overlap* ovl, int novl)
{
	FilterContext* ctx = (FilterContext*) _ctx;
	int j;

	if (ctx->downsample)
	{
		int max_rid = (int) (ctx->downsample * DB_NREADS(ctx->db) / 100.0);

		for (j = 0; j < novl; j++)
		{
			Overlap* o = ovl + j;
			if (o->aread > max_rid || o->bread > max_rid)
			{
				o->flags |= OVL_DISCARD;
			}
		}
	}

	if (ctx->trim)
	{
		for (j = 0; j < novl; j++)
		{
			trim_overlap(ctx->trim, ovl + j);
		}
	}

	if (ctx->lowCoverageFilter < ctx->hghCoverageFilter)
		filterByCoverage(ctx, ovl, novl);

	if (ctx->stitch >= 0)
	{
		int k;
		j = k = 0;

		while (j < novl)
		{
			while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
			{
				k++;
			}

			ctx->nStitched += stitch(ovl + j, k - j + 1, ctx->stitch, ctx->stitch_aggressively);

			j = k + 1;
		}
	}

	if (ctx->nMaxProperChains >= 0)
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
			while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
			{
				k++;
			}

			chain(ctx, ovl + j, k - j + 1, trimABeg, trimAEnd);
			{

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

			j = k + 1;
		}
	}

	// set filter flags
	for (j = 0; j < novl; j++)
	{
		ovl[j].flags |= filter(ctx, ovl + j);
	}

	if (ctx->removeMultiMappers)
	{
		int foundMultiMapper = 0;

		if (ctx->removeMultiMappers > 1)
		{
			bzero(ctx->cover_multi_mapper, DB_READ_MAXLEN(ctx->db));

			if (ctx->trackRepeat)
			{
				track_anno* ranno = (track_anno*) (ctx->trackRepeat->anno);
				track_data* rdata = (track_data*) (ctx->trackRepeat->data);

				track_anno ob = ranno[ovl->aread] / sizeof(track_data);
				track_anno oe = ranno[ovl->aread + 1] / sizeof(track_data);

				while (ob < oe)
				{
					memset(ctx->cover_multi_mapper + rdata[ob], 1, rdata[ob + 1] - rdata[ob]);
					ob += 2;
				}
			}
		}

		int k, l, m;
		j = k = 0;

		int trimABeg, trimAEnd;
		int trimBBeg, trimBEnd;

		int numMultiMapper = 0;
		while (j < novl)
		{
			while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
			{
				k++;
			}

			if (k > j)
			{
				foundMultiMapper = 0;

				for (l = j; l <= k && !foundMultiMapper; l++)
				{
					Overlap *o1 = ovl + l;

					if (o1->flags & OVL_DISCARD)
						continue;

					trimABeg = 0;
					trimAEnd = DB_READ_LEN(ctx->db, o1->aread);

					if (ctx->trackTrim)
						get_trim(ctx->db, ctx->trackTrim, o1->aread, &trimABeg, &trimAEnd);

					trimBBeg = 0;
					trimBEnd = DB_READ_LEN(ctx->db, o1->bread);

					if (ctx->trackTrim)
						get_trim(ctx->db, ctx->trackTrim, o1->bread, &trimBBeg, &trimBEnd);

					for (m = l + 1; m <= k && !foundMultiMapper; m++)
					{
						Overlap *o2 = ovl + m;

						if (o2->flags & OVL_DISCARD)
							continue;

						if (((o1->path.abpos == trimABeg || o1->path.bbpos == trimBBeg) && (o1->path.aepos == trimAEnd || o1->path.bepos == trimBEnd))
								&& ((o2->path.abpos == trimABeg || o2->path.bbpos == trimBBeg) && (o2->path.aepos == trimAEnd || o2->path.bepos == trimBEnd))
								&& ((abs(o1->path.abpos - o2->path.abpos) > ctx->twidth) && (abs(o1->path.aepos - o2->path.aepos) > ctx->twidth)))
						{
							foundMultiMapper = 1;
						}
					}
				}

				if (foundMultiMapper)
				{
					ctx->nMultiMapper++;
					numMultiMapper++;
					for (l = j; l <= k; l++)
					{
						if (!(ovl[l].flags & OVL_DISCARD))
						{
#ifdef VERBOSE
							printf("remove multi mapper: %d vs %d [%d, %d] %c [%d, %d]\n", ovl[l].aread, ovl[l].bread, ovl[l].path.abpos, ovl[l].path.aepos,
									(ovl[l].flags & OVL_COMP) ? 'c' : 'n', ovl[l].path.bbpos, ovl[l].path.bepos);
#endif
							ctx->nMultiMapperBases += ovl[l].path.aepos - ovl[l].path.abpos;
							ovl[l].flags |= OVL_DISCARD;
						}

						if (ctx->removeMultiMappers > 1)
						{
							memset(ctx->cover_multi_mapper + ovl[l].path.abpos, 1, ovl[l].path.aepos - ovl[l].path.abpos);
						}

					}
				}
			}
			j = k + 1;
		}
		if (ctx->removeMultiMappers > 1 && numMultiMapper)
		{
			for (j = 0; j < novl; j++)
			{
				Overlap *o = ovl + j;

				if (o->flags & OVL_DISCARD)
					continue;

				int sumMultMapBases = 0;
				for (k = o->path.abpos; k < o->path.aepos; k++)
				{
					if ((int) (ctx->cover_multi_mapper[k]))
						sumMultMapBases++;
					else
						break;
				}

				int anchorBases = MAX(500, ctx->nMinNonRepeatBases);

				if (sumMultMapBases + anchorBases >= (o->path.aepos - o->path.abpos))
				{
#ifdef VERBOSE
					printf("remove ovl, falls into multi mapper interval: %d vs %d [%d, %d] %c [%d, %d]\n", o->aread, o->bread, o->path.abpos, o->path.aepos,
							(o->flags & OVL_COMP) ? 'c' : 'n', o->path.bbpos, o->path.bepos);
#endif
					o->flags |= OVL_DISCARD;
					ctx->nMultiMapper++;
					ctx->nMultiMapperBases += o->path.aepos - o->path.abpos;
				}

			}
		}
	}

	// detect weird overlaps and discard them, i.e. those that are obviously wrong but couldn't be removed with repeat annotation or local alignment filter
	// TODO extract + create own function
	// TODO perform tests
	if (ctx->removeLowCoverageOverlaps)
	{
		int trimBeg, trimEnd;
		trimBeg = 0;
		trimEnd = DB_READ_LEN(ctx->db, ovl->aread);

		if (ctx->trackTrim)
			get_trim(ctx->db, ctx->trackTrim, ovl->aread, &trimBeg, &trimEnd);

		int entercov = 0;
		int leavecov = 0;

		for (j = 0; j < novl; j++)
		{
			Overlap* ovl_j = ovl + j;

			if (ovl_j->flags & OVL_DISCARD)
				continue;

			if (ovl_j->path.abpos <= trimBeg)
				entercov++;

			if (ovl_j->path.aepos >= trimEnd)
				leavecov++;

//			if (ovl_j->aread == 142600)
//			{
//				printf("check %d vs %d inA [%d, %d] inB [%d, %d] trimA [%d, %d] >> %d  --- %d << \n", ovl_j->aread, ovl_j->bread, ovl_j->path.abpos,
//						ovl_j->path.aepos, ovl_j->path.bbpos, ovl_j->path.bepos, trimBeg, trimEnd, entercov, leavecov);
//			}

		}

		// actually if this happens the whole aread should be discarded,
		// because it leaves regions in the aread that have no overlaps at all (depends on -u flag, number of unaligned bases, which is 0 by default)
		if (entercov && entercov < ctx->removeLowCoverageOverlaps)
		{
//			printf("found entercov %d < %d: aread %d: remove", entercov, ctx->removeLowCoverageOverlaps, ovl->aread);
			for (j = 0; j < novl; j++)
			{
				Overlap* ovl_j = ovl + j;

				if (ovl_j->flags & OVL_DISCARD)
					continue;

				if (ovl_j->path.abpos <= trimBeg)
				{
//					printf(" %d", ovl_j->bread);
					ovl_j->flags |= OVL_DISCARD;
					ctx->nLowCovALn++;
					if (ctx->fileOutDiscardedOverlaps)
					{
						fprintf(ctx->fileOutDiscardedOverlaps, "%d %d\n", ovl_j->aread, ovl_j->bread);
					}
				}
			}
//				printf("\n");
		}

		if (leavecov && leavecov < ctx->removeLowCoverageOverlaps)
		{
//			printf("found leavecov %d < %d: aread %d remove", leavecov, ctx->removeLowCoverageOverlaps, ovl->aread);
			for (j = 0; j < novl; j++)
			{
				Overlap* ovl_j = ovl + j;

				if (ovl_j->flags & OVL_DISCARD)
					continue;

				if (ovl_j->path.aepos >= trimEnd)
				{
//					printf(" %d", ovl_j->bread);
					ovl_j->flags |= OVL_DISCARD;
					ctx->nLowCovALn++;
					if (ctx->fileOutDiscardedOverlaps)
					{
						fprintf(ctx->fileOutDiscardedOverlaps, "%d %d\n", ovl_j->aread, ovl_j->bread);
					}
				}
			}
//			printf("\n");
		}
	}

	// rescue chain overlaps if a single local alignment has at least -n anchor bases
	if (ctx->nMinNonRepeatBasesChain != -1)
	{
		int k;
		j = k = 0;

		int rescComp, rescNorm;

		while (j < novl)
		{
			rescComp = rescNorm = 0;

			if (!(ovl[j].flags & OVL_DISCARD))
			{
				if (ovl[j].flags & OVL_COMP)
					rescComp++;
				else
					rescNorm++;
			}

			while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
			{
				if (!(ovl[k + 1].flags & OVL_DISCARD))
				{
					if (ovl[k + 1].flags & OVL_COMP)
						rescComp++;
					else
						rescNorm++;
				}
				k++;
			}

			if (rescComp + rescNorm < k - j + 1)
			{
				int l;
				for (l = j; l <= k; ++l)
				{
					if (ovl[l].flags & OVL_DISCARD)
					{
						if (((ovl[l].flags & OVL_COMP) && rescComp) || (!(ovl[l].flags & OVL_COMP) && rescNorm))
						{
							ovl[l].flags &= ~OVL_DISCARD;
							ovl[l].flags |= OVL_MODULE;

							ctx->nRepeatOvlsKept++;
						}
					}
				}
			}
			j = k + 1;
		}
	}

	if (ctx->remUpToXPercAln)
	{
		removeWorstAlignments(ctx, ovl, novl);
	}

	// find repeat modules and rescue overlaps

	if (ctx->rm_cov != -1)
	{
		find_repeat_modules(ctx, ovl, novl);
	}

	// evaluate remove flags
	if (ctx->removeFlags != 0)
		removeOvls(ctx, ovl, novl, ctx->removeFlags);

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

	if (ctx->includeReadFlag)
	{
		for (j = 0; j < novl; j++)
		{
			if (!(reads[ovl[j].aread].flags & READ_KEEP) && !(reads[ovl[j].bread].flags & READ_KEEP))
			{
				ovl[j].flags |= OVL_DISCARD;
			}
		}
	}

	if (ctx->fileSpanningReads)
	{
		write_spanning_reads(ctx, ovl->aread);
	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "[-vpLqTwWZ] [-dnolRsSumMfyYzNZ <int>] [-rt <track>] [-xPIaA <file>] <db> <overlaps_in> <overlaps_out>\n");

	fprintf(stderr, "options: -v ... verbose\n");
	fprintf(stderr, "         -d ... max divergence allowed [0,100]\n");
	fprintf(stderr, "         -n ... min number of non-repeat bases\n");
	fprintf(stderr, "         -o ... min overlap length\n");
	fprintf(stderr, "         -l ... min read length\n");
	fprintf(stderr, "         -p ... purge discarded overlaps\n");
	fprintf(stderr, "         -s ... stitch (%d)\n", DEF_ARG_S);
	fprintf(stderr, "         -S ... stitch aggressively (%d)\n", DEF_ARG_S);
	fprintf(stderr, "         -r ... repeat track name (%s)\n", DEF_ARG_R);
	fprintf(stderr, "         -t ... trim track name (%s)\n", DEF_ARG_T);
	fprintf(stderr, "         -T ... trim overlaps (%d)\n", DEF_ARG_TT);
	fprintf(stderr, "         -u ... max number of unaligned bases\n");
	fprintf(stderr, "         -R ... remove (multiple -R possible) ... \n");
	fprintf(stderr, "                0 stitched overlaps, i.e. overlaps with invalid trace points\n");
	fprintf(stderr, "                1 module overlaps\n");
	fprintf(stderr, "                2 trace points\n");
	fprintf(stderr, "                3 non-identity overlaps\n");
	fprintf(stderr, "                4 remove B-read repeat overlaps, if a proper overlap between A and B exists \n");
	fprintf(stderr, "                5 identity overlaps\n");
	fprintf(stderr, "         -L ... two pass processing with read caching\n");
	fprintf(stderr, "experimental features:\n");
	fprintf(stderr, "         -f ... percentage of overlaps to keep (downsampling)\n");
	fprintf(stderr, "         -m ... resolve repeat modules, pass coverage as argument\n");
	fprintf(stderr, "         -M ... -m + more aggressive module detection\n");
	fprintf(stderr, "         -x ... exclude read ids found in file\n");
	fprintf(stderr, "         -I ... include read ids found in file, all other are excluded\n");
	fprintf(stderr, "         -P ... write read ids of repeat spanners to file\n");
	fprintf(stderr, "         -z ... drop entering/leaving alignments if number is below -z <int>\n");
	fprintf(stderr,
			"         -y ... merge repeats if they are closer then -y bases apart (if distance > 100, then smaller repeats (< 100) usually from DBdust are ignored)\n");
	fprintf(stderr, "         -Y ... merge repeats with start/end position of read if repeat interval starts/ends with fewer then -Y\n");
	fprintf(stderr, "         -a ... write discarded overlaps that may not symmetrically removed to file\n");
	fprintf(stderr, "         -A ... read file of discarded overlaps and remove them symmetrically\n");
	fprintf(stderr, "         -w ... remove multi-mapper \n");
	fprintf(stderr, "         -W ... -w + remove overlaps that fall into multi-mapping alignment intervals\n");
	fprintf(stderr, "         -N ... min number of non-repeat bases for a whole alignment chain\n");
	fprintf(stderr, "         -C ... keep only the best -C alignment chains for all local alignemnts of reads and A vs B\n");
	fprintf(stderr, "                0 only THE best chain is kept, all others are discarded\n");
	fprintf(stderr, "                1 only THE best chain in forward as well as in reverse complement orientation are kept\n");
	fprintf(stderr, "                otherwise keep all chains\n");
	fprintf(stderr, "         -c ... discard overlaps that don't show given coverage range. Intention get rid of contamination by coverage\n");
	fprintf(stderr, "                e.g. -c -1000  ... discard overlaps that where A-read has lower than 1000x coverage \n");
	fprintf(stderr, "                     -c  1000  ... discard overlaps that where A-read has higher than 1000x coverage \n");
	fprintf(stderr, "         -Z ... remove at most -Z percent of the worst alignments. Set -d INT to avoid loss of good alignments. Set -z INT to avoid loss of contiguity !\n");
	fprintf(stderr, "                This option was included to get rid of low coverage repeats or random alignments, that clearly get separated by diff scores\n");
}

static int opt_repeat_count(int argc, char** argv, char opt)
{
	int i;
	int count = 0;
	for (i = 1; i < argc; i++)
	{
		char* arg = argv[i];

		if (*arg == '-')
		{
			arg += 1;

			while (*arg == opt)
			{
				count += 1;
				arg += 1;
			}

			if (count)
			{
				argv[i][2] = '\0';
				break;
			}
		}
	}

	return count;
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

	char* pathSpanningReads = NULL;
	char* pathExcludeReads = NULL;
	char* pathIncludeReads = NULL;
	char* pathOutDiscardedOvls = NULL;
	char* pathInDiscardedOvls = NULL;
	char* pcTrackRepeats = DEF_ARG_R;
	char* arg_trimTrack = DEF_ARG_T;
	int arg_purge = 0;

	fctx.fileSpanningReads = NULL;
	fctx.fMaxDiffs = -1;
	fctx.nMaxUnalignedBases = -1;
	fctx.nMinAlnLength = -1;
	fctx.nMinNonRepeatBases = -1;
	fctx.nMinNonRepeatBasesChain = -1;
	fctx.nMinReadLength = -1;
	fctx.nVerbose = 0;
	fctx.stitch = DEF_ARG_S;
	fctx.rm_cov = -1;
	fctx.rm_aggressive = 0;
	fctx.useRLoader = 0;
	fctx.do_trim = DEF_ARG_TT;
	fctx.rm_mode = 0;
	fctx.downsample = 0;
	fctx.includeReadFlag = 0;
	fctx.removeLowCoverageOverlaps = 0;

	fctx.stitch_aggressively = 0;
	fctx.removeFlags = 0;
	fctx.rm_merge = 50;
	fctx.rp_mergeTips = 0;
	fctx.discardedAreadList = NULL;
	fctx.discardedBreads = NULL;
	fctx.removeMultiMappers = 0;
	fctx.nMaxProperChains = -1;
	fctx.lowCoverageFilter = -1;
	fctx.hghCoverageFilter = -1;
	fctx.cov_read_active = NULL;
	fctx.remUpToXPercAln = 0;

	int c;

	fctx.rm_mode = opt_repeat_count(argc, argv, 'm');
	if (fctx.rm_mode == 0)
	{
		fctx.rm_mode = opt_repeat_count(argc, argv, 'M');
	}

	opterr = 0;
	while ((c = getopt(argc, argv, "TvLpwWy:z:d:n:o:l:R:s:S:u:m:M:r:t:P:x:f:I:Y:a:A:N:C:c:Z:")) != -1)
	{
		switch (c)
		{
		case 'f':
			fctx.downsample = atoi(optarg);
			break;

		case 'T':
			fctx.do_trim = 1;
			break;

		case 'x':
			pathExcludeReads = optarg;
			break;

		case 'Z':
			fctx.remUpToXPercAln = atoi(optarg);
			break;

		case 'I':
			pathIncludeReads = optarg;
			break;

		case 'P':
			pathSpanningReads = optarg;
			break;

		case 'a':
			pathOutDiscardedOvls = optarg;
			break;

		case 'c':
		{
			int tmp = atoi(optarg);
			if (tmp == 0)
			{
				fprintf(stderr, "[ERROR]: Unsupported coverage filer -c [%d].\n", tmp);
				usage();
				exit(1);
			}
			if (tmp < 0)
			{
				fctx.lowCoverageFilter = 0;
				fctx.hghCoverageFilter = abs(tmp);
			}
			else
			{
				fctx.lowCoverageFilter = tmp;
				fctx.hghCoverageFilter = INT_MAX;
			}
		}
			break;

		case 'A':
			pathInDiscardedOvls = optarg;
			break;

		case 'L':
			fctx.useRLoader = 1;
			break;

		case 'w':
			fctx.removeMultiMappers = 1;
			break;

		case 'W':
			fctx.removeMultiMappers += 2;
			break;

		case 'z':
			fctx.removeLowCoverageOverlaps = atoi(optarg);
			break;

		case 'y':
			fctx.rm_merge = atoi(optarg);
			break;

		case 'Y':
			fctx.rp_mergeTips = atoi(optarg);
			break;

		case 'M':
			fctx.rm_aggressive = 1;

			// fall through

		case 'm':
			fctx.rm_cov = atoi(optarg);
			break;

		case 'R':
		{
			int flag = atoi(optarg);
			if (flag < 0 || flag > 5)
			{
				fprintf(stderr, "[ERROR]: Unsupported -R [%d] option.\n", flag);
				usage();
				exit(1);
			}
			fctx.removeFlags |= (1 << flag);
		}
			break;

		case 'S':
			fctx.stitch_aggressively = 1;

			// fall through

		case 's':
			fctx.stitch = atoi(optarg);
			break;

		case 'v':
			fctx.nVerbose = 1;
			break;

		case 'p':
			arg_purge = 1;
			break;

		case 'd':
			fctx.fMaxDiffs = atof(optarg) / 100.0;
			break;

		case 'o':
			fctx.nMinAlnLength = atoi(optarg);
			break;

		case 'l':
			fctx.nMinReadLength = atoi(optarg);
			break;

		case 'u':
			fctx.nMaxUnalignedBases = atoi(optarg);
			break;

		case 'N':
			fctx.nMinNonRepeatBasesChain = 1;

			// fall through

		case 'n':
			fctx.nMinNonRepeatBases = atoi(optarg);
			break;

		case 'r':
			pcTrackRepeats = optarg;
			break;

		case 't':
			arg_trimTrack = optarg;
			break;

		case 'C':
			fctx.nMaxProperChains = atoi(optarg);
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

	if (fctx.lowCoverageFilter < fctx.hghCoverageFilter)
		fctx.cov_read_active = (char*) malloc(DB_READ_MAXLEN(fctx.db));

	int i;
	for (i = 0; i < DB_NREADS(&db); i++)
	{
		db.reads[i].flags = READ_NONE;
	}

	if (fctx.downsample)
	{
		if (fctx.downsample > 99 || fctx.downsample < 0)
		{
			fprintf(stderr, "invalid downsampling factor %d. must be in [1, 99]\n", fctx.downsample);
			exit(1);
		}
	}

	if (pathSpanningReads)
	{
		fctx.fileSpanningReads = fopen(pathSpanningReads, "w");

		if (fctx.fileSpanningReads == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathSpanningReads);
			exit(1);
		}
	}

	if (fctx.nMinNonRepeatBases != -1 || fctx.fileSpanningReads)
	{
		fctx.trackRepeat = track_load(&db, pcTrackRepeats);

		if (!fctx.trackRepeat)
		{
			fprintf(stderr, "could not load track %s\n", pcTrackRepeats);
			exit(1);
		}
	}

	fctx.trackTrim = track_load(&db, arg_trimTrack);

	if (!fctx.trackTrim)
	{
		fprintf(stderr, "could not load track %s\n", arg_trimTrack);
		// exit( 1 );
	}

	if (pathExcludeReads && pathIncludeReads)
	{
		fprintf(stderr, "-x and -I cannot be used at the same time!!!");
		exit(1);
	}

	if (pathExcludeReads)
	{
		FILE* fileIn = fopen(pathExcludeReads, "r");

		if (fileIn == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathExcludeReads);
			exit(1);
		}

		int* values;
		int nvalues;

		fread_integers(fileIn, &values, &nvalues);

		printf("excluding %d reads\n", nvalues);

		for (i = 0; i < nvalues; i++)
		{
			db.reads[values[i]].flags = READ_DISCARD;
		}

		free(values);
		fclose(fileIn);
	}
	if (pathIncludeReads)
	{
		FILE* fileIn = fopen(pathIncludeReads, "r");

		if (fileIn == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathIncludeReads);
			exit(1);
		}

		int* values;
		int nvalues;

		fread_integers(fileIn, &values, &nvalues);

		printf("including %d reads\n", nvalues);

		for (i = 0; i < nvalues; i++)
		{
			db.reads[values[i]].flags = READ_KEEP;
		}

		fctx.includeReadFlag = 1;
		free(values);
		fclose(fileIn);
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

	if (pathInDiscardedOvls)
	{
		FILE* fileIn = fopen(pathInDiscardedOvls, "r");

		if (fileIn == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathInDiscardedOvls);
			exit(1);
		}

		// todo check if memory consumption is reasonable for larger projects
		int ** discardedAreadList = (int **) calloc((DB_NREADS(fctx.db) + 1), sizeof(int *));
		assert(discardedAreadList);

		int aread, bread;
		int prevAread, prevBread;
		prevAread = -2;
		prevBread = -1;

		int numAreads = 0;
		int numBreads = 0;

		int line = 0;
		// search for exact number of areads and breads, to avoid a bread reallocation
//		printf("check for proper file format\n");
		while (fscanf(fileIn, "%d %d\n", &aread, &bread) == 2)
		{
			if (aread != prevAread)
				numAreads++;

			numBreads++;

			assert(prevAread <= aread);
			if (prevAread == aread)
				assert(prevBread < bread);

			assert(aread < DB_NREADS(fctx.db) + 1);

			line++;

//			printf("line %5d: prev [%7d, %7d] cur [%7d, %7d] #a %5d #b %5d\n", line, prevAread, prevBread, aread, bread, numAreads, numBreads);
			prevAread = aread;
			prevBread = bread;
		}
//		printf("#areads: %6d, #breads %6d\n", numAreads, numBreads);
		fseek(fileIn, 0L, SEEK_SET);

		int * discardedBreads = (int*) malloc(sizeof(int) * (numBreads + numAreads + 10));
		assert(discardedBreads);

		prevAread = -2;
		prevBread = -1;
		int curBreadIdx = 0;
//		printf("put overlaps into data structure\n");
		while (fscanf(fileIn, "%d %d\n", &aread, &bread) == 2)
		{
			assert(prevAread <= aread);
			if (prevAread == aread)
				assert(prevBread < bread);

			if (prevAread != aread)
			{
				if (prevAread > 0)
					discardedBreads[curBreadIdx++] = -1;

				discardedAreadList[aread] = discardedBreads + curBreadIdx;
			}
			discardedBreads[curBreadIdx++] = bread;

//			printf("line %5d: prev [%7d, %7d] cur [%7d, %7d] #a %5d #b %5d \n", line, prevAread, prevBread, aread, bread, numAreads, numBreads);
			prevAread = aread;
			prevBread = bread;
		}
		discardedBreads[curBreadIdx++] = -1;

//		printf("output data structure:\n");
//		int i;
//		for(i=0; i<DB_NREADS(fctx.db)+1; i++)
//		{
//			printf("i %d\n",i);
//			if(discardedAreadList[i])
//			{
//				printf("discard ovls from aread %7d:", i);
//				int j=0;
//				while(discardedAreadList[i][j] > -1)
//				{
//					printf(" %7d", discardedAreadList[i][j++]);
//				}
//				printf("\n");
//			}
//		}
		fctx.discardedAreadList = discardedAreadList;
		fctx.discardedBreads = discardedBreads;
		fclose(fileIn);
	}

	// passes

	if (fctx.useRLoader)
	{
		fctx.rl = rl_init(&db, 1);

		pctx = pass_init(fileOvlIn, NULL);

		pctx->data = &fctx;
		pctx->split_b = 1;
		pctx->load_trace = 0;

		pass(pctx, loader_handler);
		rl_load_added(fctx.rl);
		pass_free(pctx);
	}

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

	if (fctx.useRLoader)
	{
		rl_free(fctx.rl);
	}

	Close_DB(&db);

	if (fctx.fileSpanningReads)
	{
		fclose(fctx.fileSpanningReads);
	}

	if (fctx.fileOutDiscardedOverlaps)
	{
		fclose(fctx.fileOutDiscardedOverlaps);
	}

	fclose(fileOvlOut);
	fclose(fileOvlIn);

	return 0;
}
