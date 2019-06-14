/*
 repairs gaps and weak regions based on a read's overlaps
 and produces a new set of sequences

 Author: MARVEL Team
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <sys/param.h>

#include "lib/colors.h"
#include "lib/tracks.h"
#include "lib/pass.h"
#include "lib/oflags.h"
#include "lib/utils.h"
#include "lib.ext/types.h"
#include "lib.ext/bitarr.h"

#include "db/DB.h"
#include "dalign/align.h"

// defaults

#define DEF_ARG_X    	 	1000
#define DEF_ARG_G     	 500   // don't patch gaps above a certain size
#define DEF_ARG_Q      		28   // low quality cutoff
#define DEF_ARG_B      		 7   // border coverage
#define DEF_ARG_C       8000   // maximum chimer length
#define DEF_ARG_F 		     0   // allow fuzzy LAS chain begin/end

// settings

#define FASTA_WIDTH			 	60   // wrapping for result fasta files

#define MIN_INT_LEN      	 5   // min length of an adjusted track interval

#define MIN_SPAN      	 400   // only alignments with at least MIN_SPAN bases left and right of a segment are considering as support
#define MIN_SPAN_REPEAT	 800   // in repeats: only alignments with at least MIN_SPAN_REPEAT bases left and right of a segment are considering as support


// toggles

#undef DEBUG
#undef DEBUG_INTERVAL_ADJUSTMENT
#define DEBUG_FLIP
#undef DEBUG_CHAIN
#undef DEBUG_CHIMER

#define VERBOSE

// for getopt()

extern char* optarg;
extern int optind, opterr, optopt;

// context

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

	FILE* fileFastaOut;
	FILE* fileQvOut;
	FILE* fileTrimOut;

	// arguments

	int minlen;
	int lowq;
	int maxgap;
	int trim;
	int maxspanners;
	int minsupport;
	int fixChimers;
	int maxChimerLen;
	int minChimerBorderCov;
	int fuzzyChain;

	char *qName;
	HITS_TRACK* qtrack;
	char *trimName;
	HITS_TRACK* trimtrack;
	char *repeatName;
	HITS_TRACK* repeattrack;

	HITS_TRACK** convertTracks;
	int curctracks;
	int maxctracks;

	uint64_t stats_bases_before;      // bases before patching
	uint64_t stats_bases_after;       // bases after patching

	int discardChimericReads;

	// stats

	int num_flips;
	int num_gaps;
	int num_chimers;

	char* reada;
	char* readb;
	char* read_patched;

	char** qva;
	char** qvb;
	char** qv_patched;

	int* apatches;

	Chain *ovlChains;
	int curChains;
	int maxChains;

} FixContext;

// information on a gap/weak region

typedef struct
{
	int ab;         // a begin
	int ae;         // a end

	int bb;         // b begin
	int be;         // b end

	int diff;       // quality
	int b;          // b read id
	int support;    // how many reads support the gap
	// int span;       // reads spanning the gap

	int comp;       // complement sequence when writing
} Gap;

static void fix_pre(PassContext* pctx, FixContext* fctx)
{
#ifdef VERBOSE
	printf(ANSI_COLOR_GREEN "PASS fix" ANSI_COLOR_RESET "\n");
#endif

	fctx->twidth = pctx->twidth;

	if (!(fctx->qtrack = track_load(fctx->db, fctx->qName)))
	{
		fprintf(stderr, "failed to open track %s\n", fctx->qName);
		exit(1);
	}

	if (fctx->trimName)
	{
		if (!(fctx->trimtrack = track_load(fctx->db, fctx->trimName)))
		{
			fprintf(stderr, "failed to open track %s\n", fctx->trimName);
			exit(1);
		}
	}
	else
	{
		fctx->trimtrack = NULL;
	}

	if (fctx->repeatName)
	{
		if (!(fctx->repeattrack = track_load(fctx->db, fctx->repeatName)))
		{
			fprintf(stderr, "failed to open track %s\n", fctx->repeatName);
			exit(1);
		}
	}
	else
	{
		fctx->repeattrack = NULL;
	}

	int maxlen = fctx->db->maxlen;

	fctx->reada = New_Read_Buffer(fctx->db);
	fctx->readb = New_Read_Buffer(fctx->db);
	fctx->read_patched = malloc(maxlen * 2 + 4);

	if (fctx->fileQvOut)
	{
		fctx->qva = New_QV_Buffer(fctx->db);
		fctx->qvb = New_QV_Buffer(fctx->db);
		fctx->qv_patched = malloc(sizeof(char*) * NUM_QV_STREAMS);

		char* qvs = malloc(maxlen * 2 * NUM_QV_STREAMS);
		int i;
		for (i = 0; i < NUM_QV_STREAMS; i++)
		{
			fctx->qv_patched[i] = qvs + i * maxlen * 2;
		}
	}

	fctx->apatches = malloc((maxlen / pctx->twidth + 1) * 3 * sizeof(int));


	fctx->curChains = 0;
	fctx->maxChains = 5;
	fctx->ovlChains = (Chain*) malloc(sizeof(Chain) * fctx->maxChains);
	bzero(fctx->ovlChains, sizeof(Chain) * fctx->maxChains);
}

static void fix_post(PassContext* pctx, FixContext* fctx)
{
#ifdef VERBOSE
	printf("gaps:    %5d\n", fctx->num_gaps);
	printf("flips:   %5d\n", fctx->num_flips);
	printf("chimers: %5d\n", fctx->num_chimers);
	printf("replaced %'" PRIu64 " with %'" PRIu64 " bases\n", fctx->stats_bases_before, fctx->stats_bases_after);
#endif

	UNUSED(pctx);

	free(fctx->reada - 1);
	free(fctx->readb - 1);
	free(fctx->read_patched);

	if (fctx->fileQvOut)
	{
		Free_QV_Buffer(fctx->qva);
		Free_QV_Buffer(fctx->qvb);
		Free_QV_Buffer(fctx->qv_patched);
	}

	free(fctx->apatches);

	int i;
	for (i = 0; i < fctx->maxChains; i++)
	{
		Chain *chain = fctx->ovlChains + i;
		if (chain)
			free(chain->ovls);
		else
			break;
	}
	free(fctx->ovlChains);
}

static int cmp_gaps(const void* x, const void* y)
{
	Gap* a = (Gap*) x;
	Gap* b = (Gap*) y;

	int cmp = a->ab - b->ab;

	if (cmp == 0)
	{
		cmp = a->ae - b->ae;
	}

	if (cmp == 0)
	{
		cmp = a->diff - b->diff;
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

static int getRepeatCount(FixContext* ctx, int readID, int beg, int end)
{
	if (ctx->repeattrack == NULL)
		return 0;

	track_anno* rep_anno = ctx->repeattrack->anno;
	track_data* rep_data = ctx->repeattrack->data;

	track_anno rb, re;

	int repBases = 0;
	int rBeg, rEnd;

	if (readID < 0 || readID >= DB_NREADS(ctx->db))
	{
		fprintf(stderr, "[ERROR] - isInRepeat readID: %d out of bounds [0, %d]\n", readID, DB_NREADS(ctx->db) - 1);
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

static int contained(int ab, int ae, int bb, int be)
{
	if (ab >= bb && ae <= be)
	{
		return 1;
	}

	return 0;
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

static int getRepeatBases(FixContext *ctx, Overlap *ovl, int read)
{
	if (ctx->repeattrack == NULL)
	{
		return 0;
	}

	assert(ovl->aread == read || ovl->bread == read);

	int bLen = ovl->path.bepos - ovl->path.bbpos;

	// get repeats track
	track_anno* rep_anno = ctx->repeattrack->anno;
	track_data* rep_data = ctx->repeattrack->data;

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


static int spanningChain(FixContext*ctx, Overlap *ovls, int n, int chimerBeg, int chimerEnd, int trim_ab, int trim_ae)
{
	/// TODO hard coded
	int MIN_OVL_LOOKAHEAD = 2000;
	int MAX_OVL_LOOKAHEAD = 10000;
	int STRIDE_OVL_LOOKAHEAD = 2000;
	int MIN_ANCHOR = 800;
	int MIN_CHAIN_LEN = 3000;
	int FUZZY = ctx->fuzzyChain;

	if (n < 2)
	{
		if (!(ovls->flags & OVL_DISCARD) && ovls->path.abpos + MIN_ANCHOR < chimerBeg && ovls->path.aepos - MIN_ANCHOR > chimerEnd && ovls->path.aepos - ovls->path.abpos > MIN_CHAIN_LEN)
		{
#ifdef DEBUG_CHAIN
			printf("found spnning chain %d vs %d\n", ovls->aread, ovls->bread);
#endif
			return 1;
		}
		return 0;
	}

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

#ifdef DEBUG_CHAIN
	printf("find detect already excluded overlaps\n");
#endif
	{
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = ovls + i;

			if (ovl_i->flags & (OVL_CONT | OVL_DISCARD))
			{
				nremain--;
			}
		}
	}
#ifdef DEBUG_CHAIN
	printf("nremain %d\n", nremain);
#endif

	if(nremain == 0)
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

			if (ovl_i->flags & (OVL_CONT | OVL_DISCARD))
				continue;

			for (j = i + 1; j < n; j++)
			{
				Overlap *ovl_j = ovls + j;

				if (ovl_j->flags & (OVL_CONT | OVL_DISCARD))
					continue;

				if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos)
						&& contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
				{
					nremain--;
					ovl_j->flags |= (OVL_CONT | OVL_DISCARD);
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

	// cleanup all flags
	for (i=0; i<n; i++)
	{
		ovls[i].flags &= ~(OVL_DISCARD | OVL_TEMP | OVL_CONT);
	}

	// check if chain spans putative chimer
	Chain* bestChain = ctx->ovlChains;

	int trim_bb, trim_be;

	if (ctx->trimtrack)
	{
		get_trim(ctx->db, ctx->trimtrack, ovls->bread, &trim_bb, &trim_be);
	}
	else
	{
		trim_bb = 0;
		trim_be = DB_READ_LEN(ctx->db, ovls->bread);
	}

	if(bestChain->ovls[0]->path.abpos + MIN_ANCHOR > chimerBeg)
		return 0;

	if(bestChain->ovls[bestChain->novl-1]->path.aepos - MIN_ANCHOR < chimerEnd)
		return 0;

	int chainLen = 0;
	for(i=0; i<bestChain->novl; i++)
		chainLen+=bestChain->ovls[i]->path.aepos - bestChain->ovls[i]->path.abpos;

	if(chainLen < MIN_CHAIN_LEN)
		return 0;

	// chain must be proper

	// todo
	// 1. allow not fully proper chain? i.e. that start/end near the tips of the trim intervals????
	// 2. set a minimum overlaps length to avoid 'sparse chains' (only induced by repeats) ??
	// 3. set a minimum of non-repeat bases for a chain ???		--> DONE

	if(bestChain->ovls[0]->flags & OVL_COMP)
	{
		if(((bestChain->ovls[0]->path.abpos - FUZZY < trim_ab) || ( (blen - bestChain->ovls[0]->path.bepos) - FUZZY < trim_bb)) &&
				((bestChain->ovls[bestChain->novl-1]->path.aepos + FUZZY >= trim_ae) || ( (blen - bestChain->ovls[bestChain->novl-1]->path.bbpos) + FUZZY >= trim_be)))
		{
	#ifdef DEBUG_CHAIN
			printf("found spanning chain %d vs %d\n", ovls->aread, ovls->bread);
	#endif
			return 1;
		}
	}
	else
	{
		if(((bestChain->ovls[0]->path.abpos - FUZZY < trim_ab) || (bestChain->ovls[0]->path.bbpos - FUZZY < trim_bb)) &&
				((bestChain->ovls[bestChain->novl-1]->path.aepos + FUZZY >= trim_ae) || (bestChain->ovls[bestChain->novl-1]->path.bepos + FUZZY >= trim_be)))
		{
	#ifdef DEBUG_CHAIN
			printf("found spanning chain %d vs %d\n", ovls->aread, ovls->bread);
	#endif
			return 1;
		}
	}
	return 0;
}

static int oChainIntervalIntersection(FixContext *fctx, Overlap* ovls, int novl, int b, int e, int trim_ab, int trim_ae)
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

static int bReadIntersectionInterval(FixContext *fctx, Overlap* ovls, int novl, int b, int e)
{
	bit* visited = ba_new(DB_NREADS(fctx->db));

	int NON_REP_BASES = 10;

	int i;
	int pre = 0;
	int post= 0;
	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;

		if (ovl->path.abpos + MIN(b/2, 1000) < b && ((ovl->path.aepos - ovl->path.abpos) - getRepeatCount(fctx, ovl->aread, ovl->path.abpos, ovl->path.aepos) >= NON_REP_BASES))
		{
			ba_assign(visited, ovl->bread, TRUE);
			++pre;
		}
	}

	int intersection = 0;

	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;
		if (ovl->path.aepos - MIN((DB_READ_LEN(fctx->db, ovl->aread) - e) / 2, 1000) > e && ((ovl->path.aepos - ovl->path.abpos) - getRepeatCount(fctx, ovl->aread, ovl->path.abpos, ovl->path.aepos) >= NON_REP_BASES))
		{
			if(ba_value(visited, ovl->bread))
			{
				ba_assign(visited, ovl->bread, FALSE);
				intersection++;
			}
			++post;
		}
	}
	free(visited);

	// if there are almost no overlaps before or after chimer, then chimer detection is not reliable!
//	if(pre <= 3  || post <= 3)
//		return 99;

	return intersection;
}

static int spanners_interval(FixContext *fctx, Overlap* ovls, int novl, int b, int e)
{
	int span = 0;
	int minSpan = (e-b-getRepeatCount(fctx, ovls->aread, b, e) < 100) ? MIN_SPAN_REPEAT : MIN_SPAN;

	int i;
	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;

		if (ovl->path.abpos < b - minSpan && ovl->path.aepos > e + minSpan)
		{
			span++;
		}
	}

#ifdef DEBUG
	if (ovls->aread == 218393)
		printf("spanners_interval %d minspan %d [%d ,%d] #span: %d\n", ovls->aread, minSpan, b, e, span);
#endif

	return span;
}

static int spanners_point(FixContext *fctx, Overlap* ovls, int novl, int p, int max)
{
	int span = 0;
	int minSpan = getRepeatCount(fctx, ovls->aread, p, p + 1) ? MIN_SPAN_REPEAT : MIN_SPAN;

	int i;
	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;

		if (ovl->path.abpos < p - minSpan && ovl->path.aepos > p + minSpan)
		{
			span++;

			if (span > max)
			{
				break;
			}
		}
	}

#ifdef DEBUG
	if (ovls->aread == 218393)
		printf("spanners_point %d minspan %d point %d #span: %d max %d\n", ovls->aread, minSpan, p, span, max);
#endif

	return (span > max);
}

static int filter_chimers(FixContext* fctx, Overlap* ovls, int novl, int* trim_b, int* trim_e)
{
	int *alnBeg = (int *) malloc((DB_READ_MAXLEN(fctx->db) / fctx->twidth + 1) * sizeof(int));
	memset(alnBeg, 0, (DB_READ_MAXLEN(fctx->db) / fctx->twidth + 1) * sizeof(int));

	int *alnEnd = (int *) malloc((DB_READ_MAXLEN(fctx->db) / fctx->twidth + 1) * sizeof(int));
	memset(alnEnd, 0, (DB_READ_MAXLEN(fctx->db) / fctx->twidth + 1) * sizeof(int));

#ifdef DEBUG_CHIMER
	printf("filter_chimers aread %d #ovls %d trim [%d, %d]\n", ovls->aread, novl, *trim_b, *trim_e);
#endif
	int i;

	for (i = 1; i < novl - 1; i++)
	{
		Overlap *opre = ovls + i - 1;
		Overlap *ocur = ovls + i;
		Overlap *opas = ovls + i + 1;

		int prematEnd = (ocur->path.bepos + 1000 < DB_READ_LEN(fctx->db, ocur->bread)) && (ocur->path.aepos + 1000 < *trim_e) && (ocur->path.aepos - 100 > *trim_b);
		int prematBeg = (ocur->path.bbpos > 1000) && (ocur->path.abpos > *trim_b + 1000) && (ocur->path.abpos + 100 < *trim_e);

//		// check weird case first
		if (prematBeg && prematEnd && (opre->bread != ocur->bread) && (ocur->bread != opas->bread )) // most probably a repeat
			continue;

		if ((ocur->bread != opas->bread ) && prematEnd)
		{
			alnEnd[(ocur->path.aepos / fctx->twidth)] += 1;
#ifdef DEBUG_CHIMER
			printf("found end %d vs %d a[%d, %d] b[%d, %d] ends[%d]=%d\n", ocur->aread, ocur->bread, ocur->path.abpos, ocur->path.aepos, ocur->path.bbpos,
					ocur->path.bepos, (ocur->path.aepos / fctx->twidth), alnEnd[(ocur->path.aepos / fctx->twidth)]);
#endif
		}

		else if((opre->bread != ocur->bread) && prematBeg)
		{
			alnBeg[(ocur->path.abpos / fctx->twidth)] += 1;
#ifdef DEBUG_CHIMER
			printf("found beg %d vs %d a[%d, %d] b[%d, %d] beg[%d]=%d\n", ocur->aread, ocur->bread, ocur->path.abpos, ocur->path.aepos, ocur->path.bbpos,
									ocur->path.bepos, (ocur->path.abpos / fctx->twidth), alnBeg[(ocur->path.abpos / fctx->twidth)]);
#endif
		}
	}

	int chimBeg = -1;
	int chimEnd = -1;

	int numChim = 0;
	int remainBeg = *trim_b;
	int remainEnd = *trim_e;

	int j;
	i = j = 0;

	for (i = remainBeg / fctx->twidth; i <= remainEnd / fctx->twidth; i++)
	{
		if (alnEnd[i] >= fctx->minChimerBorderCov)
		{
			chimBeg = i;
		}
		else if (alnEnd[i] + alnEnd[i + 1] >= fctx->minChimerBorderCov)
		{
			chimBeg = ++i;
		}

		if (chimBeg < 0)
			continue;

		for (j = remainBeg / fctx->twidth; j <= remainEnd / fctx->twidth; j++)
		{
			if (alnBeg[j] >= fctx->minChimerBorderCov)
			{
				chimEnd = j;
			}

			else if (alnBeg[j] + alnBeg[j + 1] >= fctx->minChimerBorderCov)
			{
				chimEnd = ++j;
			}

			if (chimEnd < 0)
				continue;
#ifdef DEBUG_CHIMER
			printf("check putative chimer indexes [%d, %d]\n", chimBeg, chimEnd);
#endif
			chimBeg = i * fctx->twidth;
			chimEnd *= fctx->twidth;

			if (chimBeg > chimEnd)
			{
				int tmp = chimBeg;
				chimBeg = chimEnd;
				chimEnd = tmp;
			}

#ifdef DEBUG_CHIMER
			printf("check putative chimer interval [%d, %d]\n", chimBeg, chimEnd);
#endif

			if((chimEnd - chimBeg) > fctx->maxChimerLen)
			{
#ifdef DEBUG_CHIMER
				printf("chimer interval too large [%d, %d] len %d-- IGNORED\n", chimBeg, chimEnd, (chimEnd - chimBeg));
#endif
				chimEnd = -1;
				continue;

			}

			int tmpRep = getRepeatCount(fctx, ovls->aread, chimBeg, chimEnd);
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

			int crossingChains = oChainIntervalIntersection(fctx, ovls, novl, chimBeg, chimEnd, remainBeg, remainEnd);
#ifdef DEBUG_CHIMER
			printf("intersectBreads %d\n", crossingChains);
#endif
			if (crossingChains < 1)
			{
#ifdef DEBUG_CHIMER
				printf("check putative chimer (spanner < 1 true) interval [%d, %d]\n", chimBeg, chimEnd);
#endif
				if (numChim == 0)
				{
					if (chimBeg - *trim_b > *trim_e - chimEnd)
					{
						remainEnd = chimBeg;
					}
					else
					{
						remainBeg = chimEnd;
					}
				}
				else if (chimBeg < remainEnd)
				{
					if (chimBeg - remainBeg > remainEnd - chimEnd)
					{
						remainEnd = chimBeg;
					}
					else
					{
						remainBeg = chimEnd;
					}
				}
				++numChim;
#ifdef DEBUG_CHIMER
				printf("AREAD %d FOUND %d. CHIMER [%d, %d] best remaining part [%d, %d]\n", ovls->aread, numChim, chimBeg, chimEnd, remainBeg, remainEnd);
#endif
				chimEnd = -1;
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
		*trim_b = remainBeg;
		*trim_e = remainEnd;

		if(fctx->discardChimericReads)
		{
			int tmp = *trim_b;
			*trim_b = *trim_e;
			*trim_e = tmp;
		}
#ifdef DEBUG_FLIP
		printf("FINALCHIMER %d [%d, %d] remaining part [%d, %d]\n", ovls->aread, chimBeg, chimEnd, *trim_b, *trim_e);
#endif
	}

	free(alnBeg);
	free(alnEnd);
	return numChim;
}

static int filter_flips(FixContext* fctx, Overlap* ovls, int novl, int* trim_b, int* trim_e)
{
	int self_n = 0;
	int self_c = 0;
	int aread = ovls->aread;
	int trimmed = 0;

	int b = -1;
	int e = -1;

	int i;
	for (i = 0; i < novl; i++)
	{
		Overlap* ovl = ovls + i;
		int bread = ovl->bread;

		if (aread < bread)
		{
			if (b != -1)
			{
				e = i;
			}

			break;
		}
		else if (aread == bread)
		{
			if (b == -1)
			{
				b = i;
			}

			if (ovl->flags & OVL_COMP)
			{
				self_c++;
			}
			else
			{
				self_n++;
			}
		}
	}

	if (self_c == 0)
	{
		return trimmed;
	}
#ifdef DEBUG
	printf("%7d %3d %3d [%4d..%4d]\n", aread, self_c, self_n, b, e);
#endif
	int alen = DB_READ_LEN(fctx->db, aread);

	for (i = b; i < e; i++)
	{
		Overlap* ovl = ovls + i;

		if (ovl->flags & OVL_COMP)
		{
			int ab = ovl->path.abpos;
			int ae = ovl->path.aepos;

			int ab_c = alen - ovl->path.bepos;
			int ae_c = alen - ovl->path.bbpos;

			if (intersect(ab, ae, ab_c, ae_c))
			{
#ifdef DEBUG
				printf("  -> crosses diagonal %5d..%5d x %5d..%5d\n", ab, ae, ab_c, ae_c);
#endif
				ovl_trace* trace = ovl->path.trace;

				int sab = ovl->path.abpos;
				int sae = (sab / fctx->twidth + 1) * fctx->twidth;
				int sbb = ovl->path.bbpos;
				int sbe = sbb + trace[1];
#ifdef DEBUG
				printf("  -> sab %d sae %d sbb %d sbe %d\n", sab, sae, sbb, sbe);
#endif
				int j;
				for (j = 2; j < ovl->path.tlen - 2; j += 2)
				{
					//if (intersect(sab, sae, alen - sbe, alen - sbb)) // && spanners(ovls, novl, sab, sae) <= 1 )
					if (intersect(sab, sae + 1, alen - sbe, alen - sbb - 1)) // && spanners(ovls, novl, sab, sae) <= 1 ) // modified 25.06.18 change interval from [ ) to [ ]
					{
#ifdef DEBUG_FLIP
						printf("%8d CROSS @ %5d..%5d x %5d..%5d\n", aread, sab, sae, alen - sbe, alen - sbb);
#endif
						trimmed = 1;

						if (*trim_b < sab && sae < *trim_e)
						{
							if (sab - *trim_b < *trim_e - sae)
							{
								*trim_b = sae;
							}
							else
							{
								*trim_e = sab;
							}
						}
					}

					sab = sae;
					sae += fctx->twidth;

					sbb = sbe;
					sbe += trace[j + 1];
				}

				sae = ovl->path.aepos;
				sbe = ovl->path.bepos;
			}
		}
	}

	i = b;
	while (i < e - 1)
	{
		Overlap* ovl = ovls + i;
		Overlap* ovl2 = ovls + i + 1;

		if ((ovl->flags & OVL_COMP) && (ovl2->flags & OVL_COMP))
		{
			int ab = ovl->path.aepos;
			int ae = ovl2->path.abpos;
			int ab_c = alen - ovl2->path.bbpos;
			int ae_c = alen - ovl->path.bepos;

			if (intersect(ab, ae, ab_c, ae_c) && spanners_interval(fctx, ovls, novl, ab, ae) <= 1)
			{
#ifdef DEBUG_FLIP
				printf("%8d GAP   @ %5d..%5d x %5d..%5d\n", aread, ab, ae, ab_c, ae_c);
#endif
				trimmed = 1;

				int mid = (ab + ae) / 2;

				if (*trim_b < mid && mid < *trim_e)
				{
					if (mid - *trim_b < *trim_e - mid)
					{
						*trim_b = mid;
					}
					else
					{
						*trim_e = mid;
					}
				}
			}
		}

		i += 1;
	}

	return trimmed;
}

static int fix_process(void* _ctx, Overlap* ovl, int novl)
{
// #warning "REMOVE ME"
//    if ( ovl->aread != 166815 ) return 1;

	FixContext* fctx = (FixContext*) _ctx;

// #warning "REMOVE ME"
	// if ( DB_READ_LEN(fctx->db, ovl->aread) < 20000 ) return 1;

	int maxgap = fctx->maxgap;
	int lowq = fctx->lowq;
	int twidth = fctx->twidth;
	int maxspanners = fctx->maxspanners;
	int minsupport = fctx->minsupport;

	int dcur = 0;
	int dmax = 1000;
	Gap* data = malloc(sizeof(Gap) * dmax);

	track_anno* qanno = fctx->qtrack->anno;
	track_data* qdata = fctx->qtrack->data;

	track_data* qa = qdata + (qanno[ovl->aread] / sizeof(track_data));

	// get trim offsets and skip reads that get trimmed away

	int trim_ab, trim_ae;
	if (fctx->trimtrack)
	{
		get_trim(fctx->db, fctx->trimtrack, ovl->aread, &trim_ab, &trim_ae);
	}
	else
	{
		trim_ab = 0;
		trim_ae = DB_READ_LEN(fctx->db, ovl->aread);
	}

	if (trim_ab >= trim_ae)
	{
		free(data);
		return 1;
	}

	// locate missed adapters
	int flips_trim_b = trim_ab;
	int flips_trim_e = trim_ae;
	if (filter_flips(fctx, ovl, novl, &flips_trim_b, &flips_trim_e))
	{
		fctx->num_flips += 1;
	}

	trim_ab = MAX(flips_trim_b, trim_ab);
	trim_ae = MIN(flips_trim_e, trim_ae);

	if (trim_ae - trim_ab < fctx->minlen)
	{
		printf("skip read %d len %d < min Len %d due to missed adapter\n", ovl->aread, trim_ae - trim_ab, fctx->minlen);
		free(data);
		return 1;
	}

	// locate chimers
	if(fctx->fixChimers)
	{
		int chimer_trim_b = trim_ab;
		int chimer_trim_e = trim_ae;
		if (filter_chimers(fctx, ovl, novl, &chimer_trim_b, &chimer_trim_e))
		{
			fctx->num_chimers += 1;
		}

		trim_ab = MAX(chimer_trim_b, trim_ab);
		trim_ae = MIN(chimer_trim_e, trim_ae);

		if (trim_ae - trim_ab < fctx->minlen)
		{
			printf("skip read %d len %d < min Len %d due to chimer \n", ovl->aread, trim_ae - trim_ab, fctx->minlen);
			free(data);
			return 1;
		}
	}

	if(fctx->fileTrimOut)
		fprintf(fctx->fileTrimOut,"%d %d %d\n" , ovl->aread, trim_ab, trim_ae);

	// sanity check tracks

	int alen = DB_READ_LEN(fctx->db, ovl->aread);
	int nsegments = (alen + fctx->twidth - 1) / twidth;

	int ob = qanno[ovl->aread] / sizeof(track_data);
	int oe = qanno[ovl->aread + 1] / sizeof(track_data);

	if (oe - ob != nsegments)
	{
		fprintf(stderr, "read %d expected %d Q track entries, found %d\n", ovl->aread, nsegments, oe - ob);
		exit(1);
	}

	if (trim_ab < 0 || trim_ab > alen || trim_ab > trim_ae || trim_ae > alen)
	{
		fprintf(stderr, "trim interval %d..%d outside read length %d\n", trim_ab, trim_ae, alen);
		exit(1);
	}

	// locate breaks in A and move outwards to the next segment boundary

	int i;
	Overlap* ocur = ovl + 0;
	Overlap* oprev = NULL;
	for (i = 1; i < novl; i++)
	{
		oprev = ocur;
		ocur = ovl + i;

		if (oprev->bread == ocur->bread && oprev->path.aepos < ocur->path.abpos && (oprev->flags & OVL_COMP) == (ocur->flags & OVL_COMP))
		{
			if (dcur >= dmax)
			{
				dmax = dmax * 1.2 + 1000;
				data = realloc(data, sizeof(Gap) * dmax);
			}

			ovl_trace* trace_left = ovl[i - 1].path.trace;
			ovl_trace* trace_right = ovl[i].path.trace;

			int ab = (ovl[i - 1].path.aepos - 1) / twidth;
			int ae = ovl[i].path.abpos / twidth + 1;

			int j = ovl[i - 1].path.tlen - 1;

			int bb = ovl[i - 1].path.bepos - trace_left[j];
			int be = ovl[i].path.bbpos + trace_right[1];

			/*
			 while (qa[ab-1] > lowq)
			 {
			 j -= 2;
			 ab--;

			 bb -= trace_left[j];
			 }

			 j = 1;
			 while (qa[ae+1] > lowq)
			 {
			 j += 2;
			 ae++;

			 be += trace_right[j];
			 }
			 */

			if (bb >= be)
			{
				continue;
			}

			if (ovl[i].flags & OVL_COMP)
			{
				int t = bb;
				int blen = DB_READ_LEN(fctx->db, ovl[i].bread);

				bb = blen - be;
				be = blen - t;
			}

			int weak_b = 0;
			track_data* qb = qdata + (qanno[ovl[i].bread] / sizeof(track_data));
			int beg = bb / twidth;
			int end = be / twidth + 1;
			int q = 0;
			j = beg;
			while (j < end)
			{
				if (qb[j] == 0)
				{
					weak_b = 1;
				}

				q += qb[j];
				j++;
			}

			if (weak_b)
			{
				continue;
			}

			if (spanners_point(fctx, ovl, novl, ab * twidth, maxspanners) && spanners_point(fctx, ovl, novl, ae * twidth, maxspanners))
			{
				continue;
			}

			int q_total = 0;
			int q_zero = 0;
			for (j = ab + 1; j < ae - 1; j++)
			{
				q_total += 1;

				if (qa[j] == 0)
				{
					q_zero += 1;
				}
			}

			if (q_total - q_zero > 1) // TODO --- hardcoded
			{
				continue;
			}

			if ((ae - ab) * twidth * 10 < (be - bb)) // TODO --- hardcoded
			{
				continue;
			}

#ifdef DEBUG

			printf("A %7d %5d..%5d %5d..%5d -> ", ovl->aread, ovl[i - 1].path.abpos, ovl[i - 1].path.aepos, ovl[i].path.abpos, ovl[i].path.aepos);
			printf("B %7d %5d..%5d %5d..%5d | ", ovl[i].bread, ovl[i - 1].path.bbpos, ovl[i - 1].path.bepos, ovl[i].path.bbpos, ovl[i].path.bepos);
			printf("%d %d %d %d %d %d %d\n", ab, ae, q, beg, end, bb, be);

#endif

			// gap due to potential weak region in A

			data[dcur].ab = ab * twidth;
			data[dcur].ae = ae * twidth;
			data[dcur].b = ovl[i].bread;
			data[dcur].bb = bb;
			data[dcur].be = be;
			data[dcur].support = 1;
			data[dcur].diff = fctx->twidth * 1.0 * q / (be - bb);
			data[dcur].comp = (ovl[i].flags & OVL_COMP);

#ifdef DEBUG
			printf("OVL %d..%d -> %d..%d\n", ovl[i - 1].path.abpos, ovl[i - 1].path.aepos, ovl[i].path.abpos, ovl[i].path.aepos);
#endif
			dcur++;
		}
	}

	qsort(data, dcur, sizeof(Gap), cmp_gaps);

	int j = 0;

	// merge breaks located at the same position in A

	for (i = 0; i < dcur; i++)
	{
		if (data[i].support == -1)
		{
			continue;
		}

		if (maxgap != -1 && (data[i].ae - data[i].ab >= maxgap || abs(data[i].be - data[i].bb) >= maxgap))
		{
			data[i].support = -1;
			continue;
		}

		for (j = i + 1; j < dcur && data[i].ab == data[j].ab && data[i].ae == data[j].ae; j++)
		{
			if (data[j].support == -1)
			{
				continue;
			}

			if (abs((data[j].be - data[j].bb) - (data[i].be - data[i].bb)) < 50)
			{
				data[i].support += 1;
				data[j].support = -1;
			}
		}
	}

	// merge overlapping breaks

	for (i = 0; i < dcur; i++)
	{
		if (data[i].support == -1)
		{
			continue;
		}

		for (j = i + 1; j < dcur && data[i].ae > data[j].ab && data[i].ab < data[j].ae; j++)
		{
			if (data[j].support == -1)
			{
				continue;
			}

			// if ( data[i].ae - data[i].ab > data[j].ae - data[j].ab )
			if (data[i].support > data[j].support)
			{
				data[i].support += data[j].support;
				data[j].support = -1;
			}
			else
			{
				data[j].support += data[i].support;
				data[i].support = -1;
				break;
			}
		}
	}

	// filter breaks with not enough support (# B reads) or no accompanying Q drop in A
	/*
	 for (j = 0; j < dcur; j++)
	 {
	 if ( data[j].support != -1 && spanners_point(ovl, novl, data[j].ab) > maxspanners
	 && spanners_point(ovl, novl, data[j].ae) > maxspanners )
	 {
	 data[j].support = -1;
	 }
	 }
	 */

	j = 0;
	for (i = 0; i < dcur; i++)
	{
		if (data[i].support < minsupport)
		{
			continue;
		}

		int bad_q = 0;
		int k;
		for (k = data[i].ab / twidth; k < data[i].ae / twidth; k++)
		{
			if (qa[k] == 0 || qa[k] >= lowq)
			{
				bad_q = 1;
				break;
			}
		}

		if (!bad_q)
		{
			continue;
		}

		data[j] = data[i];
		j++;
	}

	dcur = j;

	// scan for bad regions ~1k from both ends

	int seg_first = trim_ab / twidth;
	int seg_last = trim_ae / twidth;

	while (qa[seg_first] == 0)
	{
		seg_first++;
	}

	while (qa[seg_last - 1] == 0)
	{
		seg_last--;
	}

	for (i = seg_first; i < seg_last; i++)
	{
		// TODO ---------

		// if (i == seg_first + BEND_SEGMENTS)
		// {
		//    i = MAX(i, seg_last - BEND_SEGMENTS - 1);
		// }

		if (qa[i] != 0 && qa[i] < lowq)
		{
			continue;
		}

		int ab = i * twidth;
		int ae = (i + 1) * twidth;

		// already covered by a break interval
		int contained = 0;
		for (j = 0; j < dcur; j++)
		{
			if (data[j].support != -1 && data[j].ab <= ab && data[j].ae >= ae)
			{
				contained = 1;
				break;
			}
		}

		if (contained)
		{
#ifdef DEBUG
			printf("%d @ %d contained\n", i, qa[i]);
#endif
			continue;
		}

		// spanners & reads starting/stopping
		int span = 0;
		int border = 0;

		float q_min = 0;
		int ovl_q_min = -1;
		int bb_q_min, be_q_min;

		for (j = 0; j < novl; j++)
		{
			// if (ovl[j].flags & OVL_COMP) continue;

			if (ovl[j].path.abpos + fctx->twidth <= ab && ovl[j].path.aepos - fctx->twidth >= ae)     // TODO --- hardcoded
			{

				if (dcur >= dmax)
				{
					dmax = dmax * 1.2 + 1000;
					data = realloc(data, sizeof(Gap) * dmax);
				}

				// locate replacement segment(s) in B

				int bb, be;

				bb = -1;
				be = ovl[j].path.bbpos;
				ovl_trace* tracerep = ovl[j].path.trace;

				int apos = ovl[j].path.abpos;
				int k = 0;
				while (apos <= ab)
				{
					apos = (apos / twidth + 1) * twidth;

					bb = be;
					be += tracerep[k + 1];
					k += 2;
				}

				assert(bb != -1);

				if (ovl[j].flags & OVL_COMP)
				{
					int t = bb;
					int blen = DB_READ_LEN(fctx->db, ovl[j].bread);

					bb = blen - be;
					be = blen - t;
				}

				// get Q in B read
				track_data* qb = qdata + (qanno[ovl[j].bread] / sizeof(track_data));
				int beg = bb / twidth;
				int end = be / twidth;
				int q = 0;
				k = beg;
				while (k < end)
				{
					if (qb[k] == 0)
					{
						q = 0;
						break;
					}

					q += qb[k];
					k++;
				}

				if (q == 0)
				{
					continue;
				}

				// printf("%d..%d %d..%d\n", bb, be, beg, end);

				float q_new = 1.0 * q / (end - beg);

				if (ovl_q_min == -1 || q_new < q_min)
				{
					bb_q_min = bb;
					be_q_min = be;
					q_min = q_new;
					ovl_q_min = j;
				}
#ifdef DEBUG
				printf("%c %8d %5d..%5d %3d..%3d %5.2f\n", ovl_q_min == j ? '*' : ' ', ovl[j].bread, bb, be, beg, end, q_new);
#endif
				span++;
			}

			if ((ovl[j].path.abpos >= ab && ovl[j].path.abpos <= ae) || (ovl[j].path.aepos >= ab && ovl[j].path.aepos <= ae))
			{
				border++;
			}
		}

		// nothing spans the bad region or nothing starts/stops there
		if (ovl_q_min == -1) // || border == 0)
		{
			continue;
		}

		// locate replacement segment in B
		Overlap* ovlrep = ovl + ovl_q_min;

#ifdef DEBUG
		printf("in %d bad Q %d @ %d..%d SPAN %d BORDER %d\n", ovl->aread, qa[i], ab, ae, span, border);

		printf("  -> using B %d..%d x %d..%d @ %.2f\n", ovlrep->path.abpos, ovlrep->path.aepos, ovlrep->path.bbpos, ovlrep->path.bepos, q_min);
#endif

		data[dcur].bb = bb_q_min;
		data[dcur].be = be_q_min;
		data[dcur].b = ovlrep->bread;
		// data[dcur].span = span;
		data[dcur].support = border;
		data[dcur].ab = ab;
		data[dcur].ae = ae;
		data[dcur].diff = q_min;
		data[dcur].comp = (ovlrep->flags & OVL_COMP);

		dcur++;
	}

	// no problems in read

	if (dcur == 0)
	{
		Load_Read(fctx->db, ovl->aread, fctx->reada, 1);

		if (trim_ae - trim_ab >= fctx->minlen)
		{
			fprintf(fctx->fileFastaOut, ">trimmed_%d source=%d", ovl->aread, ovl->aread);

			for (i = 0; i < fctx->curctracks; i++)
			{
				track_anno* anno = fctx->convertTracks[i]->anno;
				track_data* data = fctx->convertTracks[i]->data;
				track_anno ob = anno[ovl->aread] / sizeof(track_data);
				track_anno oe = anno[ovl->aread + 1] / sizeof(track_data);
				char* track = fctx->convertTracks[i]->name;

				int first = 1;
				int beg, end;
				for (; ob < oe; ob += 2)
				{
					beg = data[ob] - trim_ab;
					end = data[ob + 1] - trim_ab;

					// check trim begin
					if (end < 0)
					{
						continue;
					}

					if (beg < 0)
					{
						beg = 0;
					}

					// check trim end
					if (beg > trim_ae - trim_ab)
					{
						break;
					}

					if (end > trim_ae - trim_ab)
					{
						end = (trim_ae - trim_ab);
					}

					if (first)
					{
						fprintf(fctx->fileFastaOut, " %s=", track);
					}
					else
					{
						fprintf(fctx->fileFastaOut, ",");
					}

					fprintf(fctx->fileFastaOut, "%d,%d", beg, end);

					first = 0;
				}

			}

			fprintf(fctx->fileFastaOut, "\n");

			wrap_write(fctx->fileFastaOut, fctx->reada + trim_ab, trim_ae - trim_ab, FASTA_WIDTH);

			if (fctx->fileQvOut)
			{
				Load_QVentry(fctx->db, ovl->aread, fctx->qva, 1);

				fprintf(fctx->fileQvOut, "@fixed/%d_%d source=%d\n", 0, trim_ae - trim_ab, ovl->aread);

				for (i = 0; i < NUM_QV_STREAMS; i++)
				{
					fprintf(fctx->fileQvOut, "%.*s\n", trim_ae - trim_ab, fctx->qva[i] + trim_ab);
				}
			}
		}

		// cleanup
		free(data);

		return 1;
	}

	qsort(data, dcur, sizeof(Gap), cmp_gaps);

	// count reads that span the break
	/*
	 for (i = 0; i < novl; i++)
	 {
	 for (j = 0; j < dcur; j++)
	 {
	 if (ovl[i].path.abpos + 100 < data[j].ab && ovl[i].path.aepos - 100 > data[j].ae)       // TODO --- hardcoded
	 {
	 data[j].span += 1;
	 }
	 }
	 }
	 */

	// calculate new read length and patch segments
	Load_Read(fctx->db, ovl->aread, fctx->reada, 1);

	if (fctx->fileQvOut)
	{
		Load_QVentry(fctx->db, ovl->aread, fctx->qva, 1);
	}

	char* read = fctx->read_patched;
	char** qv = fctx->qv_patched;
	int rlen = 0;

	int ab = trim_ab;
	int ae;

#ifdef DEBUG
	printf("A %7d TRIM %5d..%5d\n", ovl->aread, trim_ab, trim_ae);
#endif

	int* apatches = fctx->apatches;
	int napatches = 0;

	for (i = 0; i < dcur; i++)
	{
		if (trim_ab > data[i].ab)
		{
			ab = data[i].ae;
			continue;
		}

		if (trim_ae < data[i].ae)
		{
			// ae = data[i].ae;
			break;
		}

		ae = data[i].ab;

		if (trim_ab < ae && trim_ab > ab)
		{
			ab = trim_ab;
		}

		// A[ab..ae]

		assert(ab <= ae);

		if (ab < ae)
		{
#ifdef DEBUG
			printf("A %7d %5d..%5d\n", ovl->aread, ab, ae);
#endif
			apatches[napatches] = ab;
			apatches[napatches + 1] = ae;
			apatches[napatches + 2] = rlen;
			napatches += 3;

			if (fctx->fileQvOut)
			{
				for (j = 0; j < NUM_QV_STREAMS; j++)
				{
					memcpy(qv[j] + rlen, fctx->qva[j] + ab, ae - ab);
				}
			}

			memcpy(read + rlen, fctx->reada + ab, ae - ab);
			rlen += ae - ab;
		}

		ab = data[i].ae;

		// B[bb..be]

		fctx->num_gaps += 1;

		int bb = data[i].bb;
		int be = data[i].be;

		fctx->stats_bases_before += data[i].ae - data[i].ab;
		fctx->stats_bases_after += data[i].be - data[i].bb;

		if (fctx->fileQvOut)
		{
			Load_QVentry(fctx->db, data[i].b, fctx->qvb, 1);

			for (j = 0; j < NUM_QV_STREAMS; j++)
			{
				if (data[i].comp)
				{
					rev(fctx->qvb[j] + bb, be - bb);
				}

				memcpy(qv[j] + rlen, fctx->qvb[j] + bb, be - bb);
			}
		}

		Load_Read(fctx->db, data[i].b, fctx->readb, 1);

		if (data[i].comp)
		{
			revcomp(fctx->readb + bb, be - bb);
		}

		memcpy(read + rlen, fctx->readb + bb, be - bb);
		rlen += be - bb;

#ifdef DEBUG
		printf("B %7d %5d..%5d (%6d) @ DIFF %3d SUPPORT %3d", data[i].b, bb, be, be - bb, data[i].diff, data[i].support);

		printf("    Q");
		for (j = data[i].ab / twidth; j < data[i].ae / twidth; j++)
		{
			printf(" %2d", qa[j]);
		}
		printf("\n");
#endif

	}

	ae = trim_ae;

	if (ab < ae)
	{
		apatches[napatches] = ab;
		apatches[napatches + 1] = ae;
		apatches[napatches + 2] = rlen;
		napatches += 3;

		if (fctx->fileQvOut)
		{
			for (j = 0; j < NUM_QV_STREAMS; j++)
			{
				memcpy(qv[j] + rlen, fctx->qva[j] + ab, ae - ab);
			}
		}

		memcpy(read + rlen, fctx->reada + ab, ae - ab);
		rlen += ae - ab;

#ifdef DEBUG
		printf("A %7d %5d..%5d\n", ovl->aread, ab, ae);
#endif
	}

#ifdef DEBUG
	printf("A %7d RLEN %5d -> %5d\n", ovl->aread, DB_READ_LEN(fctx->db, ovl->aread), rlen);
#endif

	// write patched sequence

	if (rlen >= fctx->minlen)
	{
		fprintf(fctx->fileFastaOut, ">fixed_%d source=%d", ovl->aread, ovl->aread);

#ifdef DEBUG_INTERVAL_ADJUSTMENT
		printf("\n\n");
		for (i = 0; i < napatches; i += 3)
		{
			printf("A-PATCH %5d..%5d -> %5d\n", apatches[i], apatches[i + 1], apatches[i + 2]);
		}
#endif

		// for each track
		for (i = 0; i < fctx->curctracks; i++)
		{
			track_anno* anno = fctx->convertTracks[i]->anno;
			track_data* data = fctx->convertTracks[i]->data;
			track_anno ob = anno[ovl->aread] / sizeof(track_data);
			track_anno oe = anno[ovl->aread + 1] / sizeof(track_data);
			char* track = fctx->convertTracks[i]->name;

			// adjust intervals if present
			if (ob < oe)
			{
				int first = 1;

				while (ob < oe)
				{
					int ib = data[ob];
					int ie = data[ob + 1];
					int ib_adj = -1;
					int ie_adj = -1;

					if (ie < apatches[0] || ib > apatches[napatches - 2])
					{
#ifdef DEBUG_INTERVAL_ADJUSTMENT
						printf("INTRVL  %5d..%5d -> OUTSIDE\n", ib, ie);
#endif

						ob += 2;
						continue;
					}

					for (j = 0; j < napatches; j += 3)
					{
						if (ib_adj == -1)
						{
							if (ib < apatches[j + 1])
							{
								ib_adj = MAX(ib, apatches[j]);
								ib_adj = apatches[j + 2] + (ib_adj - apatches[j]);
							}
						}

						if (ie_adj == -1)
						{
							if (ie <= apatches[j + 1])
							{
								if (ie < apatches[j] && j > 0)
								{
									ie_adj = apatches[j - 2];
									ie_adj = apatches[j - 1] + (ie_adj - apatches[j - 3]);

									break;
								}
								else if (ie > apatches[j])
								{
									ie_adj = ie;
									ie_adj = apatches[j + 2] + (ie_adj - apatches[j]);

									break;
								}
							}
						}
					}

					if (ie_adj - ib_adj > MIN_INT_LEN)
					{
#ifdef DEBUG_INTERVAL_ADJUSTMENT
						printf("INTRVL  %5d..%5d -> %5d..%5d\n", ib, ie, ib_adj, ie_adj);
#endif

						if (!first)
						{
							fprintf(fctx->fileFastaOut, ",");
						}
						else
						{
							fprintf(fctx->fileFastaOut, " %s=", track);
						}

						// sanity check
						if (ib_adj < 0 || ib_adj > rlen || ib_adj > ie_adj || ie_adj > rlen)
						{
							fprintf(stderr, "adjust interval %d..%d outside read length %d\n", ib_adj, ie_adj, rlen);
							exit(1);
						}

						fprintf(fctx->fileFastaOut, "%d,%d", ib_adj, ie_adj);

						first = 0;
					}
					else
					{
#ifdef DEBUG_INTERVAL_ADJUSTMENT
						printf("INTRVL  %5d..%5d -> SKIP\n", ib, ie);
#endif
					}

					ob += 2;
				}
			}
		}

		fprintf(fctx->fileFastaOut, "\n");

		wrap_write(fctx->fileFastaOut, read, rlen, FASTA_WIDTH);

		if (fctx->fileQvOut)
		{
			fprintf(fctx->fileQvOut, "@fixed/%d_%d source=%d\n", 0, rlen, ovl->aread);

			for (j = 0; j < NUM_QV_STREAMS; j++)
			{
				fprintf(fctx->fileQvOut, "%.*s\n", rlen, qv[j]);
			}
		}
	}

	// cleanup
	free(data);

	return 1;
}

static void usage()
{
	printf("usage: [-ladX] [-bCxQgF <int>] [-ctqr <track>] [-f <patched.quiva>] [-T <file>] <db> <in.las> <patched.fasta>\n");
	printf("       -c ... convert track intervals (multiple -c possible)\n");
	printf("       -f ... patch quality streams\n");
	printf("       -x ... min length for fixed sequences (%d)\n", DEF_ARG_X);
	printf("       -Q ... segment quality threshold (%d)\n", DEF_ARG_Q);
	printf("       -g ... max gap length for patching (%d)\n", DEF_ARG_G);
	printf("       -t ... trim reads based on a track\n");
	printf("       -q ... quality track (default: %s)\n", TRACK_Q);
	printf("       -l ... low coverage mode\n");
	printf("EXPERIMENTAL OPTIONS\n");
	printf("       -X ... fix chimeric reads in repeat regions\n");
	printf("       -r ... repeat track\n");
	printf("       -d ... discard all chimeric reads, (default: 0, i.e. keep longest part of chimeric reads)\n");
	printf("       -b ... minimum border coverage to start a chimer detection (default: %d)\n", DEF_ARG_B);
	printf("       -C ... maximum chimer length (default: %d)\n", DEF_ARG_C);
	printf("       -T ... write trim interval after gap-, cross-, and chimer-detection into a file.\n");
	printf("       -F ... allow LAS chain to have a fuzzy begin and end overlap. LAS chains are used to find spanners over putative chimeric break intervals. (default: %d)\n", DEF_ARG_F);

}

int main(int argc, char* argv[])
{
	HITS_DB db;
	PassContext* pctx;
	FixContext fctx;
	FILE* fileOvlIn;

	bzero(&fctx, sizeof(FixContext));
	fctx.db = &db;
	fctx.minlen = DEF_ARG_X;
	fctx.lowq = DEF_ARG_Q;
	fctx.maxgap = DEF_ARG_G;
	fctx.trimName = NULL;
	fctx.repeatName = NULL;
	fctx.qName = TRACK_Q;
	fctx.discardChimericReads = 0;
	fctx.maxChimerLen = DEF_ARG_C;
	fctx.minChimerBorderCov = DEF_ARG_B;
	fctx.fuzzyChain = DEF_ARG_F;
	// process arguments

	char* pathQvOut = NULL;
	char* pathTrimOut = NULL;
	int c;
	int lowc = 0;
	opterr = 0;

	while ((c = getopt(argc, argv, "dlXx:f:c:Q:g:t:q:r:b:C:T:F:")) != -1)
	{
		switch (c)
		{
		case 'l':
			lowc = 1;
			break;

		case 'X':
			fctx.fixChimers = 1;
			break;

		case 'd':
			fctx.discardChimericReads = 1;
			break;

		case 'Q':
			fctx.lowq = atoi(optarg);
			break;

		case 'F':
 			fctx.fuzzyChain = atoi(optarg);
			break;

		case 'b':
			fctx.minChimerBorderCov = atoi(optarg);
			break;

		case 'C':
			fctx.maxChimerLen = atoi(optarg);
			break;

		case 'g':
			fctx.maxgap = atoi(optarg);
			break;

		case 'x':
			fctx.minlen = atoi(optarg);
			break;

		case 'f':
			pathQvOut = optarg;
			break;

		case 'T':
			pathTrimOut = optarg;
			break;

		case 'q':
			fctx.qName = optarg;
			break;

		case 't':
			fctx.trimName = optarg;
			break;

		case 'r':
			fctx.repeatName = optarg;
			break;

		case 'c':
			if (fctx.curctracks >= fctx.maxctracks)
			{
				fctx.maxctracks += 10;
				fctx.convertTracks = realloc(fctx.convertTracks, sizeof(HITS_TRACK*) * fctx.maxctracks);
			}

			// use the HITS_TRACK* array as temporary storage of the track names

			fctx.convertTracks[fctx.curctracks] = (HITS_TRACK*) optarg;
			fctx.curctracks++;

			break;

		default:
			usage();
			exit(1);
		}
	}

	if (opterr || argc - optind != 3)
	{
		usage();
		exit(1);
	}

	char* pcPathReadsIn = argv[optind++];
	char* pcPathOverlapsIn = argv[optind++];
	char* pcPathFastaOut = argv[optind++];

	if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL)
	{
		fprintf(stderr, "could not open '%s'\n", pcPathOverlapsIn);
		exit(1);
	}

	if ((fctx.fileFastaOut = fopen(pcPathFastaOut, "w")) == NULL)
	{
		fprintf(stderr, "could not open '%s'\n", pcPathFastaOut);
		exit(1);
	}

	if (pathQvOut)
	{
		if ((fctx.fileQvOut = fopen(pathQvOut, "w")) == NULL)
		{
			fprintf(stderr, "error: could not open '%s'\n", pathQvOut);
			exit(1);
		}
	}

	if (pathTrimOut)
	{
		if ((fctx.fileTrimOut = fopen(pathTrimOut, "w")) == NULL)
		{
			fprintf(stderr, "error: could not open '%s'\n", pathTrimOut);
			exit(1);
		}
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "could not open database '%s'\n", pcPathReadsIn);
		exit(1);
	}

	int i;
	for (i = 0; i < fctx.curctracks; i++)
	{
		char* track = (char*) fctx.convertTracks[i];
		fctx.convertTracks[i] = track_load(&db, track);

		if (fctx.convertTracks[i] == NULL)
		{
			fprintf(stderr, "could not open track '%s'\n", track);
			exit(1);
		}
	}

	if (lowc)
	{
		fctx.maxspanners = 3;
		fctx.minsupport = 2;
	}
	else
	{
		fctx.maxspanners = 7;       // 10
		fctx.minsupport = 4;        // 5
	}

	// pass

	if (fctx.fileQvOut)
	{
		if (Load_QVs(&db) != 0)
		{
			fprintf(stderr, "error: failed to load QVs\n");
			exit(1);
		}
	}

	pctx = pass_init(fileOvlIn, NULL);

	pctx->split_b = 0;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &fctx;

	fix_pre(pctx, &fctx);

	pass(pctx, fix_process);

	fix_post(pctx, &fctx);

	pass_free(pctx);

	// cleanup

	if (fctx.fileQvOut)
	{
		Close_QVs(&db);
		fclose(fctx.fileQvOut);
	}

	if (fctx.fileTrimOut)
		fclose(fctx.fileTrimOut);

	Close_DB(&db);

	fclose(fileOvlIn);
	fclose(fctx.fileFastaOut);

	if (fctx.curctracks)
		free(fctx.convertTracks);

	return 0;
}
