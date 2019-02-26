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

#define REMOVE_STITCH_OVL ( 1 << 0 )
#define REMOVE_MOD_OVL ( 1 << 1 )
#define REMOVE_TP ( 1 << 2 )
#define REMOVE_NONID_OVL ( 1 << 3 )
#define REMOVE_SPECREP_OVL ( 1 << 4 )
#define REMOVE_ID_OVL ( 1 << 5 )

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
	int nMinNonRepeatBases, nMinReadLength;
	int nVerbose;
	int stitch;
	int stitch_aggressively;
	int rm_cov;        // repeat modules, coverage
	int rm_aggressive; // -M
	int do_trim;
	int downsample;
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
	HITS_TRACK* trackDust;

	int* r2bin;
	int max_r2bin;

	int useRLoader;
	TRIM* trim;
	Read_Loader* rl;

	ovl_header_twidth twidth;

	FILE* fileOutDiscardedOverlaps;
	int ** discardedAreadList;
	int * discardedBreads;

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

				if ((ovlk->flags & ignore_mask) || (ovli->flags & OVL_COMP) != (ovlk->flags & OVL_COMP))
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

		if(tmp)
		{
			*cumBases += tmp;
			if(*largest < (rEnd - rBeg))
				*largest = (rEnd - rBeg);
		}

		rb += 2;
	}
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
			b = repeats_anno[ovl->aread] / sizeof(track_data);
			e = repeats_anno[ovl->aread + 1] / sizeof(track_data);

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

	if (ctx->nVerbose)
		printf("Coverage[%d]: beg,end [%3d, %3d] avgCov %.2f\n", ovl->aread, numIncomingReads, numLeavingReads, cumOverallBases * 1.0 / aTrimLen);

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

		int err = (int) (so->path.diffs * 100.0 / MIN(so->path.aepos - so->path.abpos, so->path.bepos - so->path.bbpos));

		if (so->path.abpos <= trimABeg)
			numRemovedIncomingReads++;

		if (so->path.aepos >= trimAEnd)
			numRemovedLeavingReads++;

		removeAlnBases += so->path.aepos - so->path.abpos;

		if (removeAlnBases * 100.0 / cumOverallBases < MaxRemovedAlnBasesPerc && numIncomingReads - numRemovedIncomingReads >= MinTipCov
				&& numLeavingReads - numRemovedLeavingReads >= MinTipCov)
		{
			so->flags |= OVL_DISCARD | OVL_DIFF;
			ctx->nFilteredDiffs += 1;
			if (ctx->nVerbose)
				printf("DISCARD %d%% bad overlaps [%d, %d] a[%d, %d] b [%d, %d] %c ODIF %d\n", MaxRemovedAlnBasesPerc, so->aread, so->bread, so->path.abpos,
						so->path.aepos, so->path.bbpos, so->path.bepos, (so->flags & OVL_COMP) ? 'C' : 'N', err);
		}
		else
		{
			if (ctx->nVerbose)
				printf("DO NOT DISCARD %d%% bad overlaps [%d, %d] a[%d, %d] b [%d, %d] %c ODIF %d\n", MaxRemovedAlnBasesPerc, so->aread, so->bread, so->path.abpos,
						so->path.aepos, so->path.bbpos, so->path.bepos, (so->flags & OVL_COMP) ? 'C' : 'N', err);
			break;
		}

		if (err < maxQV)
		{
			if (ctx->nVerbose)
				printf("STOP reached maxqv of %d [%d, %d] a[%d, %d] b [%d, %d] %c\n", err, so->aread, so->bread, so->path.abpos, so->path.aepos, so->path.bbpos,
						so->path.bepos, (so->flags & OVL_COMP) ? 'C' : 'N');
			break;
		}
	}

	free(ovl_sort);
}

#define ANCHOR_INVALID 	(1 << 0)
#define ANCHOR_TRIM 		(1 << 1)
#define ANCHOR_LOWCOMP 	(1 << 2)

typedef struct
{

	int beg;
	int end;
	int flag;
} anchorItv;

//<0 The element pointed by p1 goes before the element pointed by p2
//0  The element pointed by p1 is equivalent to the element pointed by p2
//>0 The element pointed by p1 goes after the element pointed by p2


static int cmp_aIvl(const void *a, const void *b)
{
	anchorItv * a1 = (anchorItv*)a;
	anchorItv * a2 = (anchorItv*)b;

	if(a1->flag & ANCHOR_INVALID)
	{
		return 1;
	}
	if(a2->flag & ANCHOR_INVALID)
	{
		return -1;
	}
	return a1->beg - a2->beg;
}

static void analyzeRepeatIntervals(FilterContext *ctx, int aread)
{
	int trim_ab, trim_ae;
	int trim_alen;

	int arlen = DB_READ_LEN(ctx->db, aread);

	trim_ab = 0;
	trim_ae = arlen;

	if (ctx->trackTrim)
	{
		get_trim(ctx->db, ctx->trackTrim, aread, &trim_ab, &trim_ae);
		trim_alen = trim_ae - trim_ab;
	}
	else
	{
		trim_ab = 0;
		trim_ae = arlen;
		trim_alen = arlen;
	}

	int WINDOW   = 500;
	int MAXMERGE = 3000;
	int i, b, e;

	track_anno* repeats_anno = ctx->trackRepeat->anno;
	track_data* repeats_data = ctx->trackRepeat->data;

	b = repeats_anno[aread] / sizeof(track_data);
	e = repeats_anno[aread + 1] / sizeof(track_data);

	int numIntervals = (e - b + 1) + 4;
	anchorItv *uniqIntervals = malloc(sizeof(anchorItv) * numIntervals);
	bzero(uniqIntervals, sizeof(anchorItv) * numIntervals);
	int curItv = 0;

	if (b < e)
	{
		int rb1, rb2;
		int re1, re2;

		rb1 = repeats_data[b];
		re1 = repeats_data[b + 1];

		if (rb1 > 0)
		{
			uniqIntervals[curItv].beg = 0;
			uniqIntervals[curItv].end = rb1;
			curItv++;
		}

		b += 2;
		while (b < e)
		{
			rb2 = repeats_data[b];
			re2 = repeats_data[b + 1];

			uniqIntervals[curItv].beg = re1;
			uniqIntervals[curItv].end = rb2;
			curItv++;

			rb1 = rb2;
			re1 = re2;
			b += 2;
		}

		if (re1 < arlen)
		{
			uniqIntervals[curItv].beg = re1;
			uniqIntervals[curItv].end = arlen;
			curItv++;
		}

		// update unique intervals based on trim track
		if (trim_ab > 0 || trim_ae < arlen)
		{
			for (i = 0; i < numIntervals; i++)
			{
				anchorItv *a = uniqIntervals + i;

				if (trim_ab >= a->end)
				{
					a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
				}
				else if (trim_ab > a->beg)
				{
					a->beg = trim_ab;
				}

				if (a->beg >= trim_ae)
				{
					a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
				}
				else if (a->end > trim_ae)
				{
					a->end = trim_ae;
				}
			}
		}

		// update unique intervals based on low complexity and tandem repeat
		// todo hardcoded values !!!
		int predust, dust, postdust, longestDust, longestDustl, longestDustr;
		for (i = 0; i < curItv; i++)
		{
			anchorItv *a = uniqIntervals + i;

			if (a->flag & ANCHOR_INVALID)
				continue;

			if(a->end - a->beg > MAXMERGE)
				continue;

			getRepeatBasesFromInterval(ctx->trackDust, aread, a->beg, a->end, &dust, &longestDust);

			if(dust * 100.0 / (a->end - a->beg) > 50.0)
			{
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
			}
			else if(dust * 100.0 / (a->end - a->beg) > 20.0 && longestDust > 100)
			{
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
			}
			else if ((a->end - a->beg) < 100 && longestDust > 29)
			{
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
			}
			else // check if neighboring repeats end in low complexity interval
			{
				getRepeatBasesFromInterval(ctx->trackDust, aread, MAX(0, a->beg - WINDOW), a->beg, &predust, &longestDustl);
				getRepeatBasesFromInterval(ctx->trackDust, aread, a->end, MIN(a->end + WINDOW, arlen),&postdust, &longestDustr);

				if (predust > 200 || longestDustl > 100 || postdust > 200 || longestDustr > 100)
				{
					a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
				}
			}

			// merge tips if required
			if(ctx->rp_mergeTips)
			{
				for (i = 0; i < curItv; i++)
				{
					anchorItv *a = uniqIntervals + i;

					if (a->flag & ANCHOR_INVALID)
						continue;

					if(a->end < trim_ab + ctx->rp_mergeTips)
					{
						a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
					}
					else if(a->beg < trim_ab + ctx->rp_mergeTips)
					{
						a->beg = trim_ab + ctx->rp_mergeTips;
					}
					else
					{
						break;
					}
				}

				for (i = curItv-1; i >= 0; --i)
				{
					anchorItv *a = uniqIntervals + i;

					if (a->flag & ANCHOR_INVALID)
						continue;

					if(a->beg > trim_ae - ctx->rp_mergeTips)
					{
						a->flag |= (ANCHOR_TRIM | ANCHOR_INVALID);
					}
					else if(a->end > trim_ae - ctx->rp_mergeTips)
					{
						a->end = trim_ae + ctx->rp_mergeTips;
					}
					else
					{
						break;
					}
				}
			}

			printf("#LC %d %d %d f%d PRE %d %d %.2f DUST %d %d %.2f post %d %d %.2f SUM %d %d %.2f\n", aread, a->beg, a->end, a->flag, predust, a->beg - MAX(0, a->beg - WINDOW),
					predust * 100.0 / (a->beg - MAX(0, a->beg - WINDOW)), dust, a->end - a->beg, dust * 100.0 / (a->end - a->beg), postdust,
					MIN(a->end + WINDOW, arlen) - a->end, postdust * 100.0 / (MIN(a->end + WINDOW, arlen) - a->end), predust + dust + postdust,
					(a->beg - MAX(0, a->beg - WINDOW)) + (a->end - a->beg) + (MIN(a->end + WINDOW, arlen) - a->end),
					(predust + dust + postdust) * 100.0 / ((a->beg - MAX(0, a->beg - WINDOW)) + (a->end - a->beg) + (MIN(a->end + WINDOW, arlen) - a->end)));
		}
	}
	else // add full read inteval as uniq range
	{
		uniqIntervals[0].beg = trim_ab;
		uniqIntervals[0].end = trim_ae;
		curItv++;

	}

	repeats_anno = ctx->trackDust->anno;
	repeats_data = ctx->trackDust->data;

	b = repeats_anno[aread] / sizeof(track_data);
	e = repeats_anno[aread + 1] / sizeof(track_data);

	int rb, re;

	// update unique anchors with all low complexity intervals !!!
	int c = curItv;
	for (i = 0; i < c; i++)
	{
		anchorItv *a = uniqIntervals + i;

		if (a->flag & ANCHOR_INVALID)
			continue;

		printf("check valid unique region %d, %d\n", a->beg, a->end);
		while (b<e)
		{
			rb = repeats_data[b];
			re = repeats_data[b + 1];

			printf("check dust region [%d, %d]\n", rb, re);
			if(rb > a->end)
			{
				printf("dust behind unique region %d, %d\n", a->beg, a->end);
				break;
			}

			if(re < a->beg)
			{
				printf("dust before unique region %d, %d\n", a->beg, a->end);
				b+=2;
				continue;
			}

			// dust fully covers unique part
			if(rb <= a->beg && re >= a->end)
			{
				printf("dust fully covers unique region %d, %d\n", a->beg, a->end);
				a->flag |= (ANCHOR_LOWCOMP | ANCHOR_INVALID);
				break;
			}

			// dust aligns left with unique part
			if(rb <= a->beg)
			{
				a->beg = re;
			}
			// dust aligns with right unique part
			if(re >= a->end)
			{
				a->end = rb;
			}
			// dust splits uniq part, i.e. make unique part invalid an append splits to the end of uniqueIntervals
			printf("dust %d,%d splits unique range %d, %d\n", rb, re, a->beg, a->end);
			a->beg = re;

			printf("curItv %d >= numIntervals %d\n", curItv, numIntervals);
			if(curItv >= numIntervals)
			{
				numIntervals = 1.2*numIntervals + 10;
				uniqIntervals = (anchorItv*)realloc(uniqIntervals, sizeof(anchorItv)*numIntervals);
				bzero(uniqIntervals + curItv, sizeof(anchorItv)*(numIntervals-curItv));
			}
			printf("curItv %d >= numIntervals %d\n", curItv, numIntervals);
			uniqIntervals[curItv].beg = uniqIntervals[i].beg;
			uniqIntervals[curItv].end = rb;
			curItv++;

			b+=2;
		}
	}

	qsort(uniqIntervals, curItv, sizeof(anchorItv), cmp_aIvl);

	// report final unique anchors:
	int anchorbases = 0;
	printf("#anchors %d",aread);
	for (i = 0; i < curItv; i++)
	{
		anchorItv *a = uniqIntervals + i;
		if(a->flag & ANCHOR_INVALID)
			break;

		printf(" %d-%d-%d", a->beg, a->end, a->flag);
		anchorbases += a->end-a->beg;
	}
	curItv = i;
	printf(" sum n%d b%d\n", i, anchorbases);

	free(uniqIntervals);

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

	if (ctx->trackDust)
	{
		analyzeRepeatIntervals(ctx, ovl->aread);
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
			k++;
			j = k;
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
	fprintf(stderr, "[-vpLqTwWZ] [-dnolRsSumMfyYzZ <int>] [-rtD <track>] [-xPIaA <file>] <db> <overlaps_in> <overlaps_out>\n");

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
	fprintf(stderr,
			"         -Z ... remove at most -Z percent of the worst alignments. Set -d INT to avoid loss of good alignments. Set -z INT to avoid loss of contiguity !\n");
	fprintf(stderr,
			"                This option was included to get rid of low coverage repeats or random alignments, that clearly get separated by diff scores\n");
	fprintf(stderr, "         -D ... read low complexity track (dust or tan_dust)\n");
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
	char* arg_dustTrack = NULL;

	int arg_purge = 0;

	fctx.fileSpanningReads = NULL;
	fctx.fMaxDiffs = -1;
	fctx.nMaxUnalignedBases = -1;
	fctx.nMinAlnLength = -1;
	fctx.nMinNonRepeatBases = -1;
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
	fctx.remUpToXPercAln = 0;

	int c;

	fctx.rm_mode = opt_repeat_count(argc, argv, 'm');
	if (fctx.rm_mode == 0)
	{
		fctx.rm_mode = opt_repeat_count(argc, argv, 'M');
	}

	opterr = 0;
	while ((c = getopt(argc, argv, "TvLpwWy:z:d:n:o:l:R:s:S:u:m:M:r:t:P:x:f:I:Y:a:A:Z:D:")) != -1)
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

		case 'n':
			fctx.nMinNonRepeatBases = atoi(optarg);
			break;

		case 'r':
			pcTrackRepeats = optarg;
			break;

		case 't':
			arg_trimTrack = optarg;
			break;

		case 'D':
			arg_dustTrack = optarg;
			break;

		default:
			fprintf(stderr, "[ERROR] unknown option -%c\n", optopt);
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

	if (arg_dustTrack)
	{
		fctx.trackDust = track_load(&db, arg_dustTrack);

		if (!fctx.trackDust)
		{
			fprintf(stderr, "could not load track %s\n", arg_dustTrack);
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
			if (values[i] < 0 || values[i] >= DB_NREADS(&db))
			{
				fprintf(stderr, "[WARNING] LAfilter: excluding read %d not possible! Must be in range: [0, %d]", values[i], DB_NREADS(&db) - 1);
				continue;
			}
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
			if (values[i] < 0 || values[i] >= DB_NREADS(&db))
			{
				fprintf(stderr, "[WARNING] LAfilter: including read %d not possible! Must be in range: [0, %d]", values[i], DB_NREADS(&db) - 1);
				continue;
			}
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
