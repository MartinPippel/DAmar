/*******************************************************************************************
 *
 *  A->B overlaps can be split into multiple A->B records due to bad regions in either one
 *  of the reads, causing the overlapper to stop aligning. Here we look for those records
 *  (with a maximum gap of -f), join the two overlaps into a single record, discard the
 *  superfluous one, and (re-)compute pass-through points & diffs for the segments surrounding
 *  the break point.
 *
 *  Author  :  MARVEL Team
 *
 *  Date    :  November 2014
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <sys/param.h>

#include "lib/read_loader.h"
#include "lib/pass.h"
#include "lib/oflags.h"
#include "lib/colors.h"
#include "lib/utils.h"
#include "lib/tracks.h"

#include "db/DB.h"
#include "dalign/align.h"

// command line defaults

#define DEF_ARG_F   40
#define DEF_ARG_P    0
#define DEF_ARG_A  500
#define DEF_ARG_M   10
#define DEF_ARG_S  200
#define DEF_ARG_L  500
#define DEF_ARG_D  0.3

#define MAX_STITCH_ALIGN_ERROR  0.57

// switches

#define VERBOSE
#undef DEBUG_STITCH

// macros

#define OVL_STRAND(ovl) ( ( (ovl)->flags & OVL_COMP) ? 'c' : 'n' )

// context

typedef struct
{
	HITS_DB* db;                    // database

	int64 stitched;                 // how many stitch operations have occurred
	int64 stitchedLC;                 // how many stitch low complexity operations have occurred
	int fuzz;                       // max gap between the A/B overlaps
	int verbose;
	int minLength;
	int anchorBases;
	int mergeDistance;
	int minLCLenToMerge;
	float gapLenDifference;
	int useRLoader;

	ovl_header_twidth twidth;       // trace point spacing
	size_t tbytes;                  // bytes used for trace point storage

	ovl_trace* trace;               // new trace buffer for all B reads in the stitch_handler
	int tcur;                       // current position in the buffer
	int tmax;                       // size of the buffer

	Work_Data* align_work;          // working storage for the alignment module
	Alignment align;                // alignment and path record for computing the
	Path path;                      // alignment in the gap region

	Read_Loader* rl;

	HITS_TRACK* lowCompTrack;
	HITS_TRACK* repeatTrack;

} StitchContext;

// externals for getopt()

extern char *optarg;
extern int optind, opterr, optopt;

static void create_pass_through_points(StitchContext* ctx, ovl_trace* trace)
{
	int twidth = ctx->twidth;
	int a = ctx->align.path->abpos;
	int b = ctx->align.path->bbpos;
	int p, t;
	int diffs = 0;
	int matches = 0;

	int aprev = a;
	int bprev = b;

	int tcur = 0;

	for (t = 0; t < ctx->align.path->tlen; t++)
	{
		if ((p = ((ovl_trace*) (ctx->align.path->trace))[t]) < 0)
		{
			p = -p - 1;
			while (a < p)
			{
				if (ctx->align.aseq[a] != ctx->align.bseq[b])
					diffs++;
				else
					matches++;

				a += 1;
				b += 1;

				if (a % twidth == 0)
				{
					trace[tcur++] = diffs;
					trace[tcur++] = b - bprev;

					aprev = a;
					bprev = b;

#ifdef DEBUG_STITCH
					printf(" a(%4dx%4d %3d %3d %3d)", a, b, diffs, matches, tcur);
#endif
					diffs = matches = 0;
				}
			}

			diffs++;
			b += 1;
		}
		else
		{
			p--;

			while (b < p)
			{
				if (ctx->align.aseq[a] != ctx->align.bseq[b])
					diffs++;
				else
					matches++;

				a += 1;
				b += 1;

				if (a % twidth == 0)
				{
					trace[tcur++] = diffs;
					trace[tcur++] = b - bprev;

					aprev = a;
					bprev = b;

#ifdef DEBUG_STITCH
					printf(" b(%4dx%4d %3d %3d %3d)", a, b, diffs, matches, tcur);
#endif
					diffs = matches = 0;
				}
			}

			diffs++;
			a += 1;

			if (a % twidth == 0)
			{
				trace[tcur++] = diffs;
				trace[tcur++] = b - bprev;

				aprev = a;
				bprev = b;

#ifdef DEBUG_STITCH
				printf(" c(%4dx%4d %3d %3d %3d)", a, b, diffs, matches, tcur);
#endif
				diffs = matches = 0;
			}
		}
	}

	p = ctx->align.path->aepos;
	while (a < p)
	{
		if (ctx->align.aseq[a] != ctx->align.bseq[b])
			diffs++;
		else
			matches++;

		a += 1;
		b += 1;

		if (a % twidth == 0 && a != ctx->align.path->aepos)
		{
			trace[tcur++] = diffs;
			trace[tcur++] = b - bprev;

			aprev = a;
			bprev = b;

#ifdef DEBUG_STITCH
			printf(" d(%4dx%4d %3d %3d %3d)", a, b, diffs, matches, tcur);
#endif
			diffs = matches = 0;
		}
	}

	if (a != aprev)
	{
		trace[tcur++] = diffs;
		trace[tcur++] = b - bprev;

#ifdef DEBUG_STITCH
		printf(" e(%4dx%4d %3d %3d %3d)", a, b, diffs, matches, tcur);
#endif
	}
	else
	{
		trace[tcur - 1] += b - bprev;
	}

#ifdef DEBUG_STITCH
	printf("\n");
#endif
}

/*
 ensures that there is enough space in StitchContext.trace for
 an additional needed values.
 */
static void ensure_trace(StitchContext* sctx, int needed)
{
	// space needed

	fflush(stdout);
	if (needed + sctx->tcur >= sctx->tmax)
	{
		// void* traceb_old = sctx->trace;
		// void* tracee_old = sctx->trace + sctx->tmax;

		int tmax = (needed + sctx->tcur) * 2 + 1000;
		void* trace = realloc(sctx->trace, sizeof(ovl_trace) * tmax);

		/*
		 // re-adjust trace pointers from the overlap records
		 int i;
		 for (i = 0; i < novl; i++)
		 {
		 if (ovl[i].path.trace >= traceb_old && ovl[i].path.trace <= tracee_old)
		 {
		 printf("adjust %d\n", i);

		 ovl[i].path.trace = trace + (ovl[i].path.trace - traceb_old);
		 }
		 }
		 */

		sctx->tmax = tmax;
		sctx->trace = trace;
	}
}

/*
 initialises the StitchContext
 */
static void stitch_pre(PassContext* pctx, StitchContext* sctx)
{
#ifdef VERBOSE
	printf(ANSI_COLOR_GREEN "PASS stitching" ANSI_COLOR_RESET "\n");
#endif

	sctx->tbytes = pctx->tbytes;
	sctx->twidth = pctx->twidth;
	sctx->align_work = New_Work_Data();

	sctx->tmax = 2 * (( DB_READ_MAXLEN( sctx->db ) + pctx->twidth) / pctx->twidth);
	sctx->trace = (ovl_trace*) malloc(sizeof(ovl_trace) * sctx->tmax);

	sctx->align.path = &(sctx->path);
	sctx->align.aseq = New_Read_Buffer(sctx->db);
	sctx->align.bseq = New_Read_Buffer(sctx->db);
}

/*
 cleanup the StitchContext
 */
static void stitch_post(PassContext* pctx, StitchContext* sctx)
{
#ifdef VERBOSE
	printf("stitched %lld out of %lld overlaps\n", sctx->stitched, pctx->novl);

	if (sctx->lowCompTrack)
		printf("stitched low complexity %lld out of %lld overlaps\n", sctx->stitchedLC, pctx->novl);
#endif

	Free_Work_Data(sctx->align_work);
	free(sctx->align.aseq - 1);
	free(sctx->align.bseq - 1);

	free(sctx->trace);
}

static int getAnchorBases(StitchContext *ctx, int read, int from, int to)
{
	if (from >= to)
		return 0;

	track_anno* rep_anno;
	track_data* rep_data;

	// if available use repeat track to analyze anchor bases
	if (ctx->repeatTrack)
	{
		rep_anno = ctx->repeatTrack->anno;
		rep_data = ctx->repeatTrack->data;
	}
	else
	{
		if (ctx->lowCompTrack == NULL)
			return (to - from);

		rep_anno = ctx->lowCompTrack->anno;
		rep_data = ctx->lowCompTrack->data;
	}

	track_anno rb, re;

	int repBases = 0;
	int rBeg, rEnd;

	if (read < 0 || read >= DB_NREADS(ctx->db))
	{
		fprintf(stderr, "[ERROR] - getAnchorBases readID: %d out of bounds [0, %d]\n", read, DB_NREADS(ctx->db) - 1);
		fflush(stderr);
		exit(1);
	}

	// repeat bases in a-read
	rb = rep_anno[read] / sizeof(track_data);
	re = rep_anno[read + 1] / sizeof(track_data);

	while (rb < re)
	{
		rBeg = rep_data[rb];
		rEnd = rep_data[rb + 1];

		repBases += intersect(from, to, rBeg, rEnd);

		rb += 2;
	}

//    printf("getRepeatBasesFromInterval: %d %d %d --> %d", readID, beg, end, repBases);
	return MAX((to - from) - repBases, 0);
}

static int isLowComplexityBreak(StitchContext *ctx, Overlap *ovl1, Overlap *ovl2)
{

	// todo check if anchor bases are really unique bases, i.e. include repeat mask, other low complexity bases
	// todo do we need to check the alignment coordinates, anchor bases ,etc. according to bread too?
	if (ctx->lowCompTrack)
	{
		assert(ovl1->aread == ovl2->aread);
		assert(ovl1->bread == ovl2->bread);
		assert((ovl1->flags & OVL_COMP) == (ovl2->flags & OVL_COMP));

		track_anno* tanno = (track_anno*) (ctx->lowCompTrack->anno);
		track_data* tdata = (track_data*) (ctx->lowCompTrack->data);

		track_anno ob = tanno[ovl1->aread];
		track_anno oe = tanno[ovl1->aread + 1];

		if (ob < oe)
		{
			ob /= sizeof(track_data);
			oe /= sizeof(track_data);

			// try to merge close by low complexity intervals
			int merge_dist = ctx->mergeDistance;
			int anchor_dist = ctx->anchorBases;
			int minLCLenToMerge = ctx->minLCLenToMerge;

			int b = tdata[ob];
			int e = tdata[ob + 1];

//			printf("LAS %d vs %d %c a1[%d, %d] a2[%d, %d] --- b1[%d, %d] b2[%d, %d]\n", ovl1->aread, ovl1->bread, (ovl1->flags & OVL_COMP) ? 'C' : 'N',
//					ovl1->path.abpos, ovl1->path.aepos, ovl2->path.abpos, ovl2->path.aepos, ovl1->path.bbpos, ovl1->path.bepos, ovl2->path.bbpos, ovl2->path.bepos);

			ob += 2;

			while (ob < oe)
			{
//				printf("b %d, e: %d\n", tdata[ob], tdata[ob + 1]);
				if (e + merge_dist >= tdata[ob])
				{
					if (tdata[ob + 1] - tdata[ob] >= minLCLenToMerge)
					{
						e = tdata[ob + 1];
//						printf("LAstitch LC - MERGE: [%d, %d]\n", b, e);
					}
					else
					{
//						printf("LAstitch LC - IGNORE TINY REPEAT: [%d, %d] len %d < %d\n", tdata[ob], tdata[ob + 1], tdata[ob + 1] - tdata[ob],
//								minLCLenToMerge);
					}
					ob += 2;
					continue;
				}
				else
				{
					// check overlap pattern
					// 								                 XXXXX_LOW_COMPLEX_XXXXXX
					//      A read ----------------------------------------------------------
					// 		B read              --------------------------
					//  		B read              	        -------------------------------------

					if (e - b >= ctx->minLength) // LOW COMPLEXITY MUST BVE AT LEAST 500 BASES LONG ?
					{
						int validIn = 0;
						int validOut = 0;

//						printf("LAstitch LC - CHECK OVL: [%d, %d]\n", b, e);

						if (getAnchorBases(ctx, ovl1->aread, ovl1->path.abpos, b) >= anchor_dist && ovl1->path.aepos + 200 > b && ovl1->path.aepos - 25 < e)
							validIn = 1;
						if (getAnchorBases(ctx, ovl2->aread, e, ovl2->path.aepos) >= anchor_dist && ovl2->path.abpos + 25 > b && ovl2->path.abpos - 200 < e)
							validOut = 1;

						if (validIn && validOut)
						{
							// check gap length differences in both reads
							int gLenA = ovl2->path.abpos - ovl1->path.aepos;
							int gLenB = ovl2->path.bbpos - ovl1->path.bepos;

//							printf("glenA %d, glenB %d", gLenA, gLenB);
							int valid = 0;
							double diff = 1.0;
							// case 1: gap in both reads - bad region in both reads
							if (gLenA > 0 && gLenB > 0)
							{
								if (gLenA < gLenB)
								{
									diff = 1.0 - gLenA * 1.0 / gLenB;
									valid = 1;
								}
								else
								{
									diff = 1.0 - gLenB * 1.0 / gLenA;
									valid = 2;
								}
							}
							// case 2: overhang in both reads - why?
							else if (gLenA < 0 && gLenB < 0)
							{
								if (gLenA < gLenB)
								{
									diff = 1.0 - gLenB * 1.0 / gLenA;
									valid = 3;
								}
								else
								{
									diff = 1.0 - gLenA * 1.0 / gLenB;
									valid = 4;
								}

							}
							// case 3: one TE is larger then the other
							else
							{
								if (gLenA < 0)
								{
									gLenB = gLenB + abs(gLenA);
									diff = 1.0 - abs(gLenA) * 1.0 / gLenB;
									valid = 5;
								}
								else
								{
									gLenA = gLenA + abs(gLenB);
									diff = 1.0 - abs(gLenB) * 1.0 / gLenA;
									valid = 6;

								}
							}

//							printf("valid? %d, diff %.5f glen(%d, %d)", valid, diff, gLenA, gLenB);

							if (diff < ctx->gapLenDifference)
								return 1;

							return 0;
						}
					}

					b = tdata[ob];
					e = tdata[ob + 1];
//					printf("LAstitch LC - ASSIGN NEW: [%d, %d]\n", b, e);
					ob += 2;
				}
			}

			// check remaining lowComplexity interval (either there is only one, or its the last from the while loop)
			if (e - b >= ctx->minLength)
			{
//				printf("LAstitch LC - LAST CHECK OVL: [%d, %d]\n", b, e);
				int validIn = 0;
				int validOut = 0;

				if (getAnchorBases(ctx, ovl1->aread, ovl1->path.abpos, b) >= anchor_dist && ovl1->path.aepos + 200 > b && ovl1->path.aepos - 25 < e)
					validIn = 1;
				if (getAnchorBases(ctx, ovl2->aread, e, ovl2->path.aepos) >= anchor_dist && ovl2->path.abpos + 25 > b && ovl2->path.abpos - 200 < e)
					validOut = 1;

				if (validIn && validOut)
				{
					// check gap length differences in both reads
					int gLenA = ovl2->path.abpos - ovl1->path.aepos;
					int gLenB = ovl2->path.bbpos - ovl1->path.bepos;

//					printf("glenA %d, glenB %d", gLenA, gLenB);
					int valid = 0;
					double diff = 1.0;
					// case 1: gap in both reads - bad region in both reads
					if (gLenA > 0 && gLenB > 0)
					{
						if (gLenA < gLenB)
						{
							diff = 1.0 - gLenA * 1.0 / gLenB;
							valid = 1;
						}
						else
						{
							diff = 1.0 - gLenB * 1.0 / gLenA;
							valid = 2;
						}
					}
					// case 2: overhang in both reads - why?
					else if (gLenA < 0 && gLenB < 0)
					{
						if (gLenA < gLenB)
						{
							diff = 1.0 - gLenB * 1.0 / gLenA;
							valid = 3;
						}
						else
						{
							diff = 1.0 - gLenA * 1.0 / gLenB;
							valid = 4;
						}

					}
					// case 3: one TE is larger then the other
					else
					{
						if (gLenA < 0)
						{
							gLenB = gLenB + abs(gLenA);
							diff = 1.0 - abs(gLenA) * 1.0 / gLenB;
							valid = 5;
						}
						else
						{
							gLenA = gLenA + abs(gLenB);
							diff = 1.0 - abs(gLenB) * 1.0 / gLenB;
							valid = 6;

						}
					}

//					printf("valid? %d, diff %.5f glen(%d, %d)\n", valid, diff, gLenA, gLenB);

					if (diff < ctx->gapLenDifference)
						return 1;

					return 0;
				}
			}
		}
	}
	return 0;
}

//  Produce the concatentation of path1, path2, and path3 where they are known to meet at
//    the the ends of path2 which was produced by Compute-Alignment. Place this result in
//    a big growing buffer, that gets reset when Bridge is called with path1 = NULL.

static void Bridge(StitchContext *ctx, Path *path1, Path *path3)
{
	int k, k1, k2;
	int len, diff;
	ovl_trace *trace;

	Path *path2 = ctx->align.path;

	k1 = 2 * ((path2->abpos / ctx->twidth) - (path1->abpos / ctx->twidth));
	if (path2->aepos == path3->aepos)
		k2 = path3->tlen;
	else
		k2 = 2 * ((path2->aepos / ctx->twidth) - (path3->abpos / ctx->twidth));

	len = k1 + path2->tlen + (path3->tlen - k2);
#ifdef DEBUG_BRIDGE
	printf("k1: %d, k2: %d, len: %d\n", k1, k2, len);
#endif
	int tcur = ctx->tcur;
	trace = ctx->trace + tcur;

	diff = 0;
	len = 0;
	if (k1 > 0)
	{
		ovl_trace *t = (ovl_trace*) (path1->trace);
		for (k = 0; k < k1; k += 2)
		{
			trace[len++] = t[k];
			trace[len++] = t[k + 1];
			diff += t[k];
#ifdef DEBUG_BRIDGE
			printf("path1 %d %d %d %d\n", k, path1->abpos / ctx->twidth * ctx->twidth + len / 2 * ctx->twidth, t[k], t[k + 1]);
#endif
		}
	}
	if (path2->tlen > 0)
	{
		ovl_trace *t = (ovl_trace *) (path2->trace);
		for (k = 0; k < path2->tlen; k += 2)
		{
			trace[len++] = t[k];
			trace[len++] = t[k + 1];
			diff += t[k];
#ifdef DEBUG_BRIDGE
			printf("path2 %d %d %d %d\n", k, path1->abpos / ctx->twidth * ctx->twidth + len / 2 * ctx->twidth, t[k], t[k + 1]);
#endif
		}
	}
	if (k2 < path3->tlen)
	{
		ovl_trace *t = (ovl_trace*) (path3->trace);
		for (k = k2; k < path3->tlen; k += 2)
		{
			trace[len++] = t[k];
			trace[len++] = t[k + 1];
			diff += t[k];
#ifdef DEBUG_BRIDGE
			 printf("path3 %d %d %d %d\n", k, path1->abpos / ctx->twidth * ctx->twidth + len / 2 * ctx->twidth, t[k], t[k + 1]);
#endif
		}
	}

	path1->aepos = path3->aepos;
	path1->bepos = path3->bepos;
	path1->diffs = diff;
	path1->trace = (void *) (ctx->trace + tcur);
	path1->tlen = len;

	ctx->tcur += len;
}

static inline int MapToTPAbove(Path *path, int *x, int isA, int MR_tspace)
{
	ovl_trace *trace = (ovl_trace*) path->trace;
	int a, b, i;

	a = (path->abpos / MR_tspace) * MR_tspace;
	b = path->bbpos;
	for (i = 1; i < path->tlen; i += 2)
	{
		a += MR_tspace;
		b += trace[i];
		if (a > path->aepos)
			a = path->aepos;
		if (isA)
		{
			if (a >= *x)
			{
				*x = a;
				return (b);
			}
		}
		else
		{
			if (b >= *x)
			{
				*x = b;
				return (a);
			}
		}
	}
	if (isA)
	{
		*x = a;
		return (b);
	}
	else
	{
		*x = b;
		return (a);
	}
}

static inline int MapToTPBelow(Path *path, int *x, int isA, int MR_tspace)
{
	ovl_trace *trace = (ovl_trace*) path->trace;
	int a, b, i;

	a = ((path->aepos + (MR_tspace - 1)) / MR_tspace) * MR_tspace;
	b = path->bepos;
	for (i = path->tlen - 1; i >= 0; i -= 2)
	{
		a -= MR_tspace;
		b -= trace[i];
		if (a < path->abpos)
			a = path->abpos;
		if (isA)
		{
			if (a <= *x)
			{
				*x = a;
				return (b);
			}
		}
		else
		{
			if (b <= *x)
			{
				*x = b;
				return (a);
			}
		}
	}
	if (isA)
	{
		*x = a;
		return (b);
	}
	else
	{
		*x = b;
		return (a);
	}
}

static int Check_Bridge(Path *path, int MR_tspace)
{
	uint16 *trace = (uint16 *) path->trace;
	int i;
	int cumDiff = 0;
	int cntBadSegments = 0;

	int apos = path->abpos;
	int bpos = path->bbpos;

	for (i = 0; i < path->tlen; i += 2)
	{
		apos = (apos / MR_tspace) * MR_tspace + MR_tspace;
		bpos += trace[i + 1];

		if (MR_tspace <= TRACE_XOVR && (trace[i] > 250 || trace[i + 1] > 250))
			return (1);

		cumDiff += trace[i];

		if (1.0 * trace[i] / trace[i + 1] > MAX_STITCH_ALIGN_ERROR && i + 2 < path->tlen) // ignore last segment
		{
#ifdef DEBUG_BRIDGE
			printf("found bad segment at (a: %d -> b: %d) erate: %f\n", apos, bpos, 1.0 * trace[i] / trace[i + 1]);
#endif
			cntBadSegments++;
		}

	}

	if (cumDiff != path->diffs || cntBadSegments > 0)
		return (1);

	return (0);
}

static int Compute_Bridge_Path(StitchContext *ctx, Overlap *ovl1, Overlap *ovl2)
{
	Path *apath;
	int ain, aout;
	int bin, bout;
	int err;
	int i, j, p;
	ovl_trace *trk;

	Read_Loader* rl = ctx->rl;

	assert(ovl1->aread == ovl2->aread);
	assert(ovl1->bread == ovl2->bread);
	assert((ovl1->flags & OVL_COMP) == (ovl2->flags & OVL_COMP));

	if (ctx->useRLoader)
	{
		rl_load_read(rl, ovl1->aread, ctx->align.aseq, 0);
		rl_load_read(rl, ovl1->bread, ctx->align.bseq, 0);
	}
	else
	{
		Load_Read(ctx->db, ovl1->aread, ctx->align.aseq, 0);
		Load_Read(ctx->db, ovl1->bread, ctx->align.bseq, 0);
	}

	Path *path1 = &(ovl1->path);
	Path *path2 = &(ovl2->path);

	apath = ctx->align.path;

	ain = path2->abpos;
	aout = path1->aepos;

	bin = MapToTPBelow(path1, &ain, 1, ctx->twidth);
	bout = MapToTPAbove(path2, &aout, 1, ctx->twidth);

#ifdef DEBUG_BRIDGE
	printf("\n  Tangle [%5d..%5d] vs [%5d..%5d] \n", path1->abpos, path1->aepos, path2->abpos, path2->aepos);
	printf("         [%5d..%5d] vs [%5d..%5d] \n", path1->bbpos, path1->bepos, path2->bbpos, path2->bepos);
	printf("      (%d,%d) to (%d,%d)\n", ain, bin, aout, bout);
	fflush(stdout);
#endif

	for (i = 0; bin > ovl2->path.bbpos && i < 10; i++)
	{
		ain = (ain - ctx->twidth) / ctx->twidth * ctx->twidth;
		aout = (aout + ctx->twidth) / ctx->twidth * ctx->twidth;

		if (ain < ovl1->path.abpos || aout > ovl2->path.aepos)
			return 1;

		bin = MapToTPBelow(path1, &ain, 1, ctx->twidth);
		bout = MapToTPAbove(path2, &aout, 1, ctx->twidth);

#ifdef DEBUG_BRIDGE
		printf("\n  Tangle [%5d..%5d] vs [%5d..%5d] \n", path1->abpos, path1->aepos, path2->abpos, path2->aepos);
		printf("         [%5d..%5d] vs [%5d..%5d] \n", path1->bbpos, path1->bepos, path2->bbpos, path2->bepos);
		printf("      (%d,%d) to (%d,%d)\n", ain, bin, aout, bout);
		fflush(stdout);
#endif
	}

	if (bout < bin)
		return 1;

	apath->abpos = ain - 2 * ctx->twidth;
	apath->aepos = aout + 2 * ctx->twidth;
	apath->bbpos = MapToTPBelow(path1, &(apath->abpos), 1, ctx->twidth);
	apath->bepos = MapToTPAbove(path2, &(apath->aepos), 1, ctx->twidth);

	if ((ovl1->flags & OVL_COMP))
	{
		Complement_Seq(ctx->align.bseq, DB_READ_LEN(ctx->db, ovl1->bread));
		int p = apath->bbpos;
		apath->bbpos = DB_READ_LEN(ctx->db, ovl1->bread) - apath->bepos;
		apath->bepos = DB_READ_LEN(ctx->db, ovl1->bread) - p;
	}

#ifdef DEBUG_BRIDGE
	printf("align path:  (%d,%d) to (%d,%d)\n", apath->abpos, apath->bbpos, apath->aepos, apath->bepos);
	fflush(stdout);
#endif

	if (Compute_Alignment(&(ctx->align), ctx->align_work, DIFF_TRACE, ctx->twidth))
		printf("Compute_Trace_ALL failed\n");

#ifdef DEBUG_BRIDGE
	printf("TRACE diffs: %d, len %d\n", ctx->align.path->diffs, ctx->align.path->tlen);
	for (j = 0; j < ctx->align.path->tlen; j += 2)
	{
		printf(" %d %d\n", ((ovl_trace*) (ctx->align.path->trace))[j], ((ovl_trace*) (ctx->align.path->trace))[j + 1]);
	}
	printf("\n");
	fflush(stdout);
#endif

	trk = (ovl_trace *) apath->trace;
	if (ovl1->flags & OVL_COMP)
	{
		j = apath->tlen - 2;
		i = 0;
		while (i < j)
		{
			p = trk[i];
			trk[i] = trk[j];
			trk[j] = p;
			p = trk[i + 1];
			trk[i + 1] = trk[j + 1];
			trk[j + 1] = p;
			i += 2;
			j -= 2;
		}

		int p = apath->bbpos;
		apath->bbpos = DB_READ_LEN(ctx->db, ovl1->bread) - apath->bepos;
		apath->bepos = DB_READ_LEN(ctx->db, ovl1->bread) - p;
	}

	bin = apath->bbpos;
	bout = apath->bepos;
	err = apath->diffs;

	p = 2 * (ain / ctx->twidth - apath->abpos / ctx->twidth);
	for (i = 0; i < p; i += 2)
	{
		bin += trk[i + 1];
		err -= trk[i];
	}

	p = 2 * (apath->aepos / ctx->twidth - aout / ctx->twidth);
	for (i = ctx->align.path->tlen, p = i - p; i > p; i -= 2)
	{
		bout -= trk[i - 1];
		err -= trk[i - 2];
	}

#ifdef DEBUG_BRIDGE
	printf("      (%d,%d) to (%d,%d)\n", ain, bin, aout, bout);
	printf("  Box %d vs %d -> %d %d%%\n", aout - ain, bout - bin, err, (200 * err) / ((aout - ain) + (bout - bin)));
	fflush(stdout);
#endif

	return 0;
}

static int stitch_handler2(void* _ctx, Overlap* ovl, int novl)
{
	StitchContext* ctx = (StitchContext*) _ctx;

	if (novl < 2)
	{
		return 1;
	}

	if (ovl->aread == ovl->bread)
	{
		return 1;
	}

	int fuzz = ctx->fuzz;

	int i, j, k;
	int ab2, ae1, ae2;
	int bb2, be1, be2;
	int ab1, bb1;

	ctx->tcur = 0;

	int tsum = 0;
	for (i = 0; i < novl; i++)
	{
		tsum += ovl[i].path.tlen;
	}

	ensure_trace(ctx, tsum * 4);

	for (i = 0; i < novl; i++)
	{
		if (ovl[i].flags & OVL_DISCARD)
			continue;

		// b = ovl[i].bread;

		ab1 = ovl[i].path.abpos;
		ae1 = ovl[i].path.aepos;

		bb1 = ovl[i].path.bbpos;
		be1 = ovl[i].path.bepos;

		for (k = i + 1; k < novl; k++)
		{
			if ((ovl[k].flags & OVL_DISCARD) || (ovl[i].flags & OVL_COMP) != (ovl[k].flags & OVL_COMP))
			{
				continue;
			}

			ab2 = ovl[k].path.abpos;
			ae2 = ovl[k].path.aepos;

			bb2 = ovl[k].path.bbpos;
			be2 = ovl[k].path.bepos;

			if (ab1 > ab2 || ae1 > ae2 || bb1 > bb2 || be1 > be2)
				continue;

			// check always for anchor bases
			// todo check bread as well???
			if (getAnchorBases(ctx, ovl[i].aread, ovl[i].path.abpos, ovl[i].path.aepos) < ctx->anchorBases ||
					getAnchorBases(ctx, ovl[k].aread, ovl[k].path.abpos, ovl[k].path.aepos) < ctx->anchorBases)
				continue;

			// todo check also bread?
			int lowCompBreak = isLowComplexityBreak(ctx, ovl + i, ovl + k);

#ifdef DEBUG_STITCH
			if (lowCompBreak)
			{
				printf("ovl %d vs %d a len %d blen %d ratio: %f\n", ovl[i].aread, ovl[i].bread, abs(ae1 - ab2), abs(be1 - bb2),
						fabs(1 - abs(ae1 - ab2) / (1.0 * abs(be1 - bb2))));
			}
#endif

			if ((abs(ae1 - ab2) < fuzz && abs(be1 - bb2) < fuzz && abs((ae1 - ab2) - (be1 - bb2)) < fuzz)          			 // added 2016-02-17
			|| lowCompBreak)	 // added 2018-02-14
			{
				// try to stitch

#ifdef DEBUG_STITCH
				printf("overlap i %d vs %d [%d, %d] [%d, %d] %c\n", ovl[i].aread, ovl[i].bread, ovl[i].path.abpos, ovl[i].path.aepos, ovl[i].path.bbpos,
						ovl[i].path.bepos, ovl[i].flags & OVL_COMP ? 'c' : 'n');
#endif
				int cntBadSegments = 0;
				int apos = ovl[i].path.abpos;
				int bpos = ovl[i].path.bbpos;
				ovl_trace* trace = ovl[i].path.trace;
				for (j = 0; j < ovl[i].path.tlen; j += 2)
				{
					apos = (apos / ctx->twidth) * ctx->twidth + ctx->twidth;
					bpos += trace[j + 1];

#ifdef DEBUG_STITCH
					if (j > 0 && j % 10 == 0)
						printf("\n");
					printf("(%3d, %3d (%.2f) %5d, %5d) ", trace[j + 1], trace[j], 100.0 * trace[j] / trace[j + 1], apos, bpos);
#endif
					if (1.0 * ctx->trace[j] / ctx->trace[j + 1] > MAX_STITCH_ALIGN_ERROR && j + 2 < ctx->tcur) // ignore last segment
					{
#ifdef DEBUG_STITCH
						printf("found bad segment at %d erate: %f\n", apos, 1.0 * trace[j] / trace[j + 1]);
#endif
						cntBadSegments++;
					}
				}
#ifdef DEBUG_STITCH
				printf("\n");

				printf("overlap k %d vs %d [%d, %d] [%d, %d] %c\n", ovl[k].aread, ovl[k].bread, ovl[k].path.abpos, ovl[k].path.aepos, ovl[k].path.bbpos,
						ovl[k].path.bepos, ovl[k].flags & OVL_COMP ? 'c' : 'n');
#endif
				cntBadSegments = 0;
				apos = ovl[k].path.abpos;
				bpos = ovl[k].path.bbpos;
				trace = ovl[k].path.trace;
				for (j = 0; j < ovl[k].path.tlen; j += 2)
				{
					apos = (apos / ctx->twidth) * ctx->twidth + ctx->twidth;
					bpos += trace[j + 1];

#ifdef DEBUG_STITCH
					if (j > 0 && j % 10 == 0)
						printf("\n");
					printf("(%3d, %3d (%.2f) %5d, %5d) ", trace[j + 1], trace[j], 100.0 * trace[j] / trace[j + 1], apos, bpos);
#endif
					if (1.0 * ctx->trace[j] / ctx->trace[j + 1] > MAX_STITCH_ALIGN_ERROR && j + 2 < ctx->tcur) // ignore last segment
					{
#ifdef DEBUG_STITCH
						printf("found bad segment at %d erate: %f\n", apos, 1.0 * trace[j] / trace[j + 1]);
#endif
						cntBadSegments++;
					}
				}
#ifdef DEBUG_STITCH
				printf("\n\n");
#endif

				int tcur_start = ctx->tcur;

				if (Compute_Bridge_Path(ctx, ovl + i, ovl + k))
					continue;

				if (Check_Bridge(ctx->align.path, ctx->twidth))
				{
#ifdef DEBUG_BRIDGE
					printf("Check_Bridge FAILED: %d vs %d [%d, %d] [%d,%d] ---- [%d, %d] [%d,%d]\n", ovl[i].aread, ovl[i].bread, ovl[i].path.abpos, ovl[i].path.aepos,
							ovl[i].path.bbpos, ovl[i].path.bepos, ovl[k].path.abpos, ovl[k].path.aepos, ovl[k].path.bbpos, ovl[k].path.bepos);
#endif
					continue;
				}

				Bridge(ctx, &(ovl[i].path), &(ovl[k].path));

				cntBadSegments = 0;
				apos = ovl[i].path.abpos;
				bpos = ovl[i].path.bbpos;
				for (j = tcur_start; j < ctx->tcur; j += 2)
				{
					apos = (apos / ctx->twidth) * ctx->twidth + ctx->twidth;
					bpos += ctx->trace[j + 1];

#ifdef DEBUG_STITCH
					if (j > 0 && j % 10 == 0)
					printf("\n");
					printf("(%3d, %3d (%.2f) %5d, %5d) ", ctx->trace[j + 1], ctx->trace[j], 100.0 * ctx->trace[j] / ctx->trace[j + 1], apos, bpos);
#endif
					if (1.0 * ctx->trace[j] / ctx->trace[j + 1] > MAX_STITCH_ALIGN_ERROR && j + 2 < ctx->tcur) // ignore last segment
					{
#ifdef DEBUG_STITCH
						printf("found bad segment at %d erate: %f\n", apos, 1.0 * ctx->trace[j] / ctx->trace[j + 1]);
#endif
						cntBadSegments++;
					}
				}
#ifdef DEBUG_STITCH
				printf("\n\n");
#endif
				assert(bpos == ovl[k].path.bepos);

				ovl[k].flags |= OVL_DISCARD | OVL_STITCH;
				ae1 = ae2;
				be1 = be2;

				if (lowCompBreak)
					ctx->stitchedLC++;
				else
					ctx->stitched++;
			}
		}
	}
	return 1;
}

/*
 called for each distinct A->B overlap pair
 */
static int stitch_handler(void* _ctx, Overlap* ovl, int novl)
{
	StitchContext* ctx = (StitchContext*) _ctx;

	Read_Loader* rl = ctx->rl;

	if (novl < 2)
	{
		return 1;
	}

	if (ovl->aread == ovl->bread)
	{
		return 1;
	}

	int fuzz = ctx->fuzz;

	int i, j, k;
	int ab2, ae1, ae2;
	int bb2, be1, be2;
	int ab1, bb1;

	ctx->tcur = 0;

	int tsum = 0;
	for (i = 0; i < novl; i++)
	{
		tsum += ovl[i].path.tlen;
	}

	ensure_trace(ctx, tsum * 4);

	for (i = 0; i < novl; i++)
	{
		if (ovl[i].flags & OVL_DISCARD)
			continue;

		// b = ovl[i].bread;

		ab1 = ovl[i].path.abpos;
		ae1 = ovl[i].path.aepos;

		bb1 = ovl[i].path.bbpos;
		be1 = ovl[i].path.bepos;

		for (k = i + 1; k < novl; k++)
		{
			if ((ovl[k].flags & OVL_DISCARD) || (ovl[i].flags & OVL_COMP) != (ovl[k].flags & OVL_COMP))
			{
				continue;
			}

			ab2 = ovl[k].path.abpos;
			ae2 = ovl[k].path.aepos;

			bb2 = ovl[k].path.bbpos;
			be2 = ovl[k].path.bepos;

			// todo check also bread?
			int lowCompBreak = isLowComplexityBreak(ctx, ovl + i, ovl + k);

#ifdef DEBUG_STITCH
			if (lowCompBreak)
			{
				printf("ovl %d vs %d a len %d blen %d ratio: %f\n", ovl[i].aread, ovl[i].bread, abs(ae1 - ab2), abs(be1 - bb2),
						fabs(1 - abs(ae1 - ab2) / (1.0 * abs(be1 - bb2))));
			}
#endif

			if ((abs(ae1 - ab2) < fuzz && abs(be1 - bb2) < fuzz && abs((ae1 - ab2) - (be1 - bb2)) < fuzz)          			 // added 2016-02-17
			|| lowCompBreak)	 // added 2018-02-14
			{
				int segb = (MIN(ae1, ab2) - (ctx->twidth / 2)) / ctx->twidth;
				int sege = (MAX(ae1, ab2) + (ctx->twidth / 2)) / ctx->twidth;

				{
#ifdef DEBUG_STITCH
					printf("TRACE diffs: %d, len %d\n", ovl[i].path.diffs, ovl[i].path.tlen);
#endif
					int apos = ovl[i].path.abpos;
					int bpos = ovl[i].path.bbpos;
					for (j = 0; j < ovl[i].path.tlen; j += 2)
					{
						apos = (apos / ctx->twidth) * ctx->twidth + ctx->twidth;
						bpos += ((ovl_trace*) (ovl[i].path.trace))[j + 1];
#ifdef DEBUG_STITCH
						printf(" %d %d, %d %d\n", ((ovl_trace*) (ovl[i].path.trace))[j + 1], ((ovl_trace*) (ovl[i].path.trace))[j], apos, bpos);
#endif
					}

#ifdef DEBUG_STITCH
					printf("TRACE diffs: %d, len %d\n", ovl[k].path.diffs, ovl[k].path.tlen);
#endif
					apos = ovl[k].path.abpos;
					bpos = ovl[k].path.bbpos;
					for (j = 0; j < ovl[k].path.tlen; j += 2)
					{
						apos = (apos / ctx->twidth) * ctx->twidth + ctx->twidth;
						bpos += ((ovl_trace*) (ovl[k].path.trace))[j + 1];
#ifdef DEBUG_STITCH
						printf(" %d %d, %d %d\n", ((ovl_trace*) (ovl[k].path.trace))[j + 1], ((ovl_trace*) (ovl[k].path.trace))[j], apos, bpos);
#endif
					}
				}

#ifdef DEBUG_STITCH
				printf("[DEBUG] segB %3d segE %3d\n", segb, sege);
#endif
				//TODO experimental - increase segb,sege interval by 1 trace space to the left and right (if possible)
				{
					int shift = 1;
					segb = MAX(ab1 / ctx->twidth, segb - shift);
					sege = MIN(ae2 / ctx->twidth, sege + shift);
				}
#ifdef DEBUG_STITCH
				printf("[DEBUG] update segB %3d segE %3d\n", segb, sege);
#endif

				ovl_trace* tracei = ovl[i].path.trace;
				ovl_trace* tracek = ovl[k].path.trace;
				int align_bb, align_be, seg, apos, bpos;

				assert(segb <= sege);

				int tcur = ctx->tcur;
				int tcur_start = tcur;

				if (ctx->verbose)
				{
					if (lowCompBreak)
						printf("STITCH LC ");
					else
						printf("STITCH    ");

					printf("%8d %2d @ %5d..%5d -> %8d @ %5d..%5d %c\n", ovl[i].aread, i, ab1, ae1, ovl[i].bread, bb1, be1, OVL_STRAND(ovl + i));

					printf("                %2d @ %5d..%5d -> %8d @ %5d..%5d %c\n", k, ab2, ae2, ovl[k].bread, bb2, be2, OVL_STRAND(ovl + k));
				}

#ifdef DEBUG_STITCH
				char* color1 = (ae1 > ab2) ? ANSI_COLOR_RED : "";
				char* color2 = (be1 > bb2) ? ANSI_COLOR_RED : "";

				printf("[DFEBUG] STITCH %8d %2d @ %5d..%s%5d" ANSI_COLOR_RESET " -> %8d @ %5d..%s%5d" ANSI_COLOR_RESET " %c", ovl[i].aread, i, ab1, color1, ae1,
						ovl[i].bread, bb1, color2, be1, OVL_STRAND(ovl + i));

				apos = ovl[i].path.abpos;
				bpos = ovl[i].path.bbpos;
				for (j = 0; j < ovl[i].path.tlen; j += 2)
				{
					if (j == ovl[i].path.tlen - 2)
					{
						apos = ovl[i].path.aepos;
					}
					else
					{
						apos = (apos / ctx->twidth + 1) * ctx->twidth;
					}

					bpos += tracei[j + 1];

					if (j >= ovl[i].path.tlen - 6)
					printf(" (%3d, %3d, %5d, %5d)", tracei[j + 1], tracei[j], apos, bpos);
				}
				printf("\n");

				printf("                %2d @ %s%5d" ANSI_COLOR_RESET "..%5d -> %8d @ %s%5d" ANSI_COLOR_RESET "..%5d %c", k, color1, ab2, ae2, ovl[k].bread,
						color2, bb2, be2, OVL_STRAND(ovl + k));

				apos = ovl[k].path.abpos;
				bpos = ovl[k].path.bbpos;
				for (j = 0; j < ovl[k].path.tlen; j += 2)
				{
					if (j == ovl[i].path.tlen - 2)
					{
						apos = ovl[i].path.aepos;
					}
					else
					{
						apos = (apos / ctx->twidth + 1) * ctx->twidth;
					}
					bpos += tracek[j + 1];

					if (j < 6)
					printf(" (%3d, %3d, %5d, %5d)", tracek[j + 1], tracek[j], apos, bpos);
				}
				printf("\n");

#endif

				bpos = ovl[i].path.bbpos;
				for (seg = ovl[i].path.abpos / ctx->twidth, j = 0; j < ovl[i].path.tlen && seg < segb; j += 2, seg++)
				{
////					ctx->trace[tcur++] = tracei[j];
////					ctx->trace[tcur++] = tracei[j + 1];

					bpos += tracei[j + 1];
				}

				align_bb = bpos;
				align_be = 0;

				for (j = segb; j <= sege; j++)
				{
//					ctx->trace[tcur++] = 0;
//					ctx->trace[tcur++] = 0;
				}

				bpos = ovl[k].path.bbpos;

				for (seg = ovl[k].path.abpos / ctx->twidth, j = 0; j < ovl[k].path.tlen; j += 2, seg++)
				{
					if (seg == sege)
					{
						align_be = bpos + tracek[j + 1];
					}
					else if (seg > sege)
					{
//						ctx->trace[tcur++] = tracek[j];
//						ctx->trace[tcur++] = tracek[j + 1];
					}

					bpos += tracek[j + 1];
				}

				if ((align_bb >= align_be)
//						||
//						(fabs(1.0 - ((((sege + 1) * ctx->twidth) - (segb * ctx->twidth)) / (1.0 *(align_be - align_bb)))) > ctx->gapLenDifference + 0.05)
				)
				{
#ifdef DEBUG_STITCH
					printf("Stitching failed due to gap length differences: %f > %f\n",
							fabs(1.0 - ((((sege + 1) * ctx->twidth) - (segb * ctx->twidth)) / (1.0 * (align_be - align_bb)))), ctx->gapLenDifference + 0.05);
#endif
					continue;
				}

#ifdef DEBUG_STITCH
				printf("[DEBUG] %3d..%3d %5d..%5d %5d..%5d\n", segb, sege, segb * ctx->twidth,
						MIN((sege + 1) * ctx->twidth, DB_READ_LEN(ctx->db, ovl[i].aread)), align_bb, align_be);
#endif

				if (ctx->useRLoader)
				{
					rl_load_read(rl, ovl[i].aread, ctx->align.aseq, 0);
					rl_load_read(rl, ovl[i].bread, ctx->align.bseq, 0);
				}
				else
				{
					Load_Read(ctx->db, ovl[i].aread, ctx->align.aseq, 0);
					Load_Read(ctx->db, ovl[i].bread, ctx->align.bseq, 0);
				}

				if ((ovl[i].flags & OVL_COMP))
				{
					Complement_Seq(ctx->align.bseq, DB_READ_LEN(ctx->db, ovl[i].bread));
				}

				ctx->align.alen = DB_READ_LEN(ctx->db, ovl[i].aread);
				ctx->align.blen = DB_READ_LEN(ctx->db, ovl[i].bread);

				ctx->align.path->abpos = MIN(ae1, ab2) - 2 * ctx->twidth; //MAX(segb * ctx->twidth, ab1);
				ctx->align.path->aepos = MAX(ae1, ab2) + 2 * ctx->twidth;
				ctx->align.path->bbpos = MapToTPBelow(&(ovl[i].path), &(ctx->align.path->abpos), 1, ctx->twidth);
				ctx->align.path->bepos = MapToTPAbove(&(ovl[k].path), &(ctx->align.path->aepos), 1, ctx->twidth);

//				ctx->align.path->diffs = (ctx->align.path->aepos - ctx->align.path->abpos) + (ctx->align.path->bepos - ctx->align.path->bbpos);

#ifdef DEBUG_STITCH
				printf("Compute_Alignment: %d - %d and %d - %d\n", ctx->align.path->abpos, ctx->align.path->aepos, ctx->align.path->bbpos, ctx->align.path->bepos);
#endif

				if (Compute_Alignment(&(ctx->align), ctx->align_work, DIFF_TRACE, ctx->twidth))
					printf("Compute_Trace_ALL failed\n");

#ifdef DEBUG_STITCH
				printf("TRACE diffs: %d, len %d\n", ctx->align.path->diffs, ctx->align.path->tlen);
				for (j = 0; j < ctx->align.path->tlen; j += 2)
				{
					printf(" %d %d\n", ((ovl_trace*) (ctx->align.path->trace))[j], ((ovl_trace*) (ctx->align.path->trace))[j + 1]);
				}
				printf("\n");
				fflush(stdout);
#endif

//				int falk = tcur_start + 2 * (segb - ovl[i].path.abpos / ctx->twidth);
//				int test = ctx->align.path->abpos;
//				for (j = 0; j < ctx->align.path->tlen; j+=2)
//				{
//					ctx->trace[falk++] = ((ovl_trace*) (ctx->align.path->trace))[j];
//					ctx->trace[falk++] = ((ovl_trace*) (ctx->align.path->trace))[j+1];
//					test+=ctx->twidth;
//				}

				//create_pass_through_points(ctx, ctx->trace + tcur_start + 2 * (segb - ovl[i].path.abpos / ctx->twidth));
				Bridge(ctx, &(ovl[i].path), &(ovl[k].path));

#ifdef DEBUG_STITCH
				printf("TRACE");
				for (j = 0; j < ctx->align.path->tlen; j++)
				{
					printf(" %d", ((ovl_trace*) (ctx->align.path->trace))[j]);
				}
				printf("\n");
#endif

				int cntBadSegments = 0;
				apos = ovl[i].path.abpos;
				bpos = ovl[i].path.bbpos;
				for (j = tcur_start; j < ctx->tcur; j += 2)
				{
					apos = (apos / ctx->twidth) * ctx->twidth + ctx->twidth;
					bpos += ctx->trace[j + 1];

#ifdef DEBUG_STITCH
					if (j > 0 && j % 10 == 0)
					printf("\n");
					printf("(%3d, %3d (%.2f) %5d, %5d) ", ctx->trace[j + 1], ctx->trace[j], 100.0 * ctx->trace[j] / ctx->trace[j + 1], apos, bpos);
#endif
					if (1.0 * ctx->trace[j] / ctx->trace[j + 1] > MAX_STITCH_ALIGN_ERROR && j + 2 < tcur) // ignore last segment
					{
#ifdef DEBUG_STITCH
						printf("found bad segment at %d erate: %f\n", apos, 1.0 * ctx->trace[j] / ctx->trace[j + 1]);
#endif
						cntBadSegments++;
					}
				}
#ifdef DEBUG_STITCH
				printf("\n\n");
#endif
				assert(bpos == ovl[k].path.bepos);

				// trace points can overflow if we stitch across a big gap
				int bOverflow = 0;
				if (ctx->tbytes == sizeof(uint8))
				{
					for (j = 0; j < tcur; j += 2)
					{
						if (ctx->trace[j + 1] > 255)
						{
							bOverflow = 1;
							break;
						}
					}
				}

				if (!bOverflow && cntBadSegments == 0)
				{
					ctx->tcur = tcur;

					ovl[i].path.aepos = ae1 = ae2;
					ovl[i].path.bepos = be1 = be2;
					ovl[i].path.diffs += ovl[k].path.diffs;

					ovl[i].path.trace = ctx->trace + tcur_start;
					ovl[i].path.tlen = ctx->tcur - tcur_start;

					ovl[k].flags |= OVL_DISCARD | OVL_STITCH;

					if (lowCompBreak)
						ctx->stitchedLC++;
					else
						ctx->stitched++;
				}
				else
				{
					if(ctx->verbose)
						printf("overflow recognized or bad #segements %d > 0 \n", cntBadSegments);
				}
			}

		}
	}

	return 1;
}

static int loader_handler(void* _ctx, Overlap* ovl, int novl)
{
	StitchContext* ctx = (StitchContext*) _ctx;
	Read_Loader* rl = ctx->rl;

	if (novl < 2)
	{
		return 1;
	}

	int fuzz = ctx->fuzz;

	int i, k; // , b;
	int ab2, ae1; // , ae2;
	int bb2, be1; // , be2;
	// int ab1, bb1;

	for (i = 0; i < novl; i++)
	{
		if (ovl[i].flags & OVL_DISCARD)
			continue;

		// b = ovl[i].bread;

		// ab1 = ovl[i].path.abpos;
		ae1 = ovl[i].path.aepos;

		// bb1 = ovl[i].path.bbpos;
		be1 = ovl[i].path.bepos;

		for (k = i + 1; k < novl; k++)
		{
			if ((ovl[k].flags & OVL_DISCARD) || (ovl[i].flags & OVL_COMP) != (ovl[k].flags & OVL_COMP))
			{
				continue;
			}

			ab2 = ovl[k].path.abpos;
			// ae2 = ovl[k].path.aepos;

			bb2 = ovl[k].path.bbpos;
			// be2 = ovl[k].path.bepos;

			int lowCompBreak = 0;

			if (ctx->lowCompTrack)
				lowCompBreak = isLowComplexityBreak(ctx, ovl + i, ovl + k);

			if ((abs(ae1 - ab2) < fuzz && abs(be1 - bb2) < fuzz) || lowCompBreak)
			{
				rl_add(rl, ovl[i].aread);
				rl_add(rl, ovl[i].bread);
			}
		}
	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "usage: [-pvL] [-amlf <int>] [-rt <track>] <db> <ovl.in> <ovl.out>\n");
	fprintf(stderr, "options: -v ... verbose\n");
	fprintf(stderr, "         -f ... fuzzing for stitch (%d)\n", DEF_ARG_F);
	fprintf(stderr, "         -L ... two pass processing with read caching\n");
	fprintf(stderr, "         -p ... purge discarded overlaps\n");
	fprintf(stderr, "\nEXPERIMENTAL\n");
	fprintf(stderr, "                try to stitch overlaps, that break in low complexity regions. Only works if trace spacing is >= 126\n");
	fprintf(stderr, "         -t ... low complexity (LC) or tandem repeat (TR) interval track (e.g. dust, tan)\n");
	fprintf(stderr, "         -r ... repeat track - is used to evaluate unique anchor bases\n");
	fprintf(stderr, "         -a ... anchor bases left and right of LC/TR interval (default: %d)\n", DEF_ARG_A);
	fprintf(stderr, "         -m ... merge LC/TR intervals that are closer then m bases (default: %d)\n", DEF_ARG_M);
	fprintf(stderr, "         -l ... consider only LC/TR intervals of length greater than -l bases (default: %d)\n", DEF_ARG_L);
	fprintf(stderr, "         -M ... merge only LC/TR intervals that are at least -M bases long (default: %d)\n", DEF_ARG_S);
	fprintf(stderr,
			"         -d ... gap difference in percent (default: %.2f). Stitch only overlaps where gap size difference in A-read and B-read is lower then -m %%\n",
			DEF_ARG_D);
}

int main(int argc, char* argv[])
{
	HITS_DB db;
	StitchContext sctx;
	PassContext* pctx;
	FILE* fileOvlIn;
	FILE* fileOvlOut;

	bzero(&sctx, sizeof(StitchContext));
	sctx.fuzz = DEF_ARG_F;
	sctx.anchorBases = DEF_ARG_A;
	sctx.mergeDistance = DEF_ARG_M;
	sctx.minLength = DEF_ARG_L;
	sctx.minLCLenToMerge = DEF_ARG_S;
	sctx.gapLenDifference = DEF_ARG_D;
	sctx.verbose = 0;
	sctx.db = &db;
	sctx.useRLoader = 0;

	char* arg_lowCompTrack = NULL;
	char* arg_repeatTrack = NULL;

	// process arguments

	int arg_purge = DEF_ARG_P;

	opterr = 0;

	int c;
	while ((c = getopt(argc, argv, "Lpvf:t:a:m:l:d:r:M:")) != -1)
	{
		switch (c)
		{
		case 'p':
			arg_purge = 1;
			break;

		case 'L':
			sctx.useRLoader = 1;
			break;

		case 'v':
			sctx.verbose = 1;
			break;

		case 'f':
			sctx.fuzz = atoi(optarg);
			break;

		case 'M':
			sctx.minLCLenToMerge = atoi(optarg);
			break;

		case 'a':
			sctx.anchorBases = atoi(optarg);
			break;

		case 'm':
			sctx.mergeDistance = atoi(optarg);
			break;

		case 'l':
			sctx.minLength = atoi(optarg);
			break;

		case 'd':
			sctx.gapLenDifference = atof(optarg);
			break;

		case 't':
			arg_lowCompTrack = optarg;
			break;

		case 'r':
			arg_repeatTrack = optarg;
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
		fprintf(stderr, "could not open '%s'\n", pcPathOverlapsIn);
		exit(1);
	}

	if ((fileOvlOut = fopen(pcPathOverlapsOut, "w")) == NULL)
	{
		fprintf(stderr, "could not open '%s'\n", pcPathOverlapsOut);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "could not open database '%s'\n", pcPathReadsIn);
		exit(1);
	}

	if (sctx.fuzz < 0)
	{
		fprintf(stderr, "invalid fuzzing value of %d\n", sctx.fuzz);
		exit(1);
	}

	if (arg_lowCompTrack != NULL)
	{
		sctx.lowCompTrack = track_load(&db, arg_lowCompTrack);

		if (!sctx.lowCompTrack)
		{
			fprintf(stderr, "could not load track %s\n", arg_lowCompTrack);
			exit(1);
		}
	}

	if (arg_repeatTrack != NULL)
	{
		sctx.repeatTrack = track_load(&db, arg_repeatTrack);

		if (!sctx.repeatTrack)
		{
			fprintf(stderr, "could not load track %s\n", arg_repeatTrack);
			exit(1);
		}
	}

	if (sctx.useRLoader)
	{
		// collect read ids for loading

		sctx.rl = rl_init(&db, 1);

		pctx = pass_init(fileOvlIn, NULL);

		pctx->data = &sctx;
		pctx->split_b = 1;
		pctx->load_trace = 0;

		pass(pctx, loader_handler);

		rl_load_added(sctx.rl);

		pass_free(pctx);
	}

	// process overlaps

	pctx = pass_init(fileOvlIn, fileOvlOut);

	pctx->data = &sctx;

	pctx->split_b = 1;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->write_overlaps = 1;
	pctx->purge_discarded = arg_purge;

	stitch_pre(pctx, &sctx);

	pass(pctx, stitch_handler2);

	stitch_post(pctx, &sctx);

	// cleanup

	if (sctx.useRLoader)
	{
		rl_free(sctx.rl);
	}

	Close_DB(&db);

	pass_free(pctx);

	fclose(fileOvlIn);
	fclose(fileOvlOut);

	return 0;
}

