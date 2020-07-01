/*******************************************************************************************
 *
 *  Trims contigs based on valid contig chain overlaps 
 *
 *  Author :  DAmar Team
 *
 *  Date   :  May 2020
 *
 *******************************************************************************************/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>
#include <ctype.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/pass.h"
#include "lib/read_loader.h"
#include "lib/tracks.h"
#include "lib/trim.h"
#include "lib/utils.h"

#include "dalign/align.h"
#include "db/DB.h"
#include "CTtrim.h"

#define MIN_BIONANO_GAP_SIZE 13
#define TRIM_OFFSET 100
#define FUZZY_BASES 1500
#define FASTA_LINEWIDTH 80
#define MAX_TANDEMTRIM_PERC 70

#define DEBUG_MASKING
#undef DEBUG_MASKING2

void ensureLASchainBuffer(TrimEvidence *t, int numNewElements)
{
	assert(t != NULL);

	if (t->nLASchains + abs(numNewElements) >= t->maxLASchains)
	{
		int i = t->maxLASchains * 1.1 + MAX(numNewElements, 10);
		t->chains = (LASchain*) realloc(t->chains, sizeof(LASchain) * i);
		assert(t->chains != NULL);
		bzero(t->chains + t->maxLASchains, sizeof(LASchain) * (i - t->maxLASchains));
		t->maxLASchains = i;
	}
}

void ensureBionanoGapBuffer(TrimEvidence *t, int numNewElements)
{
	assert(t != NULL);

	if (t->nBioNanoGaps + abs(numNewElements) >= t->maxBionanoGaps)
	{
		int i = t->maxBionanoGaps * 1.1 + MAX(numNewElements, 10);
		t->gaps = (BionanoGap*) realloc(t->gaps, sizeof(BionanoGap) * i);
		assert(t->gaps != NULL);
		bzero(t->gaps + t->maxBionanoGaps, sizeof(BionanoGap) * (i - t->maxBionanoGaps));
		t->maxBionanoGaps = i;
	}
}

void addBionanoGAPInfoToTrimEvidence(TrimContext *ctx, int contigA, int aPartBeg, int aPartEnd, int contigB, int bPartBeg, int bPartEnd, int AdjustedGapLength)
{
	TrimEvidence *ta = find_TrimEvidence(ctx, contigA, contigB);
	TrimEvidence *tb = find_TrimEvidence(ctx, contigB, contigA);

	if (ta == NULL)
	{
		printf("[Warning] addBionanoGAPInfoToTrimEvidence: Could not find Bionano gap feature between contig %d and contig %d.\n", contigA, contigB);
		return;
	}

	if (tb == NULL)
	{
		printf("[Warning] addBionanoGAPInfoToTrimEvidence: Could not find Bionano gap feature between contig %d and contig %d.\n", contigB, contigA);
		return;
	}

	// check if the same gap feature is already present: it must be present;
	int i;
	BionanoGap *b;
	for (i = 0; i < ta->nBioNanoGaps; i++)
	{
		b = ta->gaps + i;

		if (((b->aBeg == aPartBeg && b->aEnd == aPartEnd) || (b->aBeg == aPartEnd && b->aEnd == aPartBeg)) && ((b->bBeg == bPartBeg && b->bEnd == bPartEnd) || (b->bBeg == bPartEnd && b->bEnd == bPartBeg)))
		{
			break;
		}
	}

	if (i == ta->nBioNanoGaps)
	{
		char *aName = getContigName(ctx, contigA);
		char *bName = getContigName(ctx, contigB);

		printf("[ERROR] - addBionanoGAPInfoToTrimEvidence 1: Cannot find bionano gap feature: Contig %d (%s) and Contig %d (%s): a[%d, %d] b[%d, %d] gapLen %d\n", contigA, aName, contigB, bName, aPartBeg, aPartEnd, bPartBeg, bPartEnd, AdjustedGapLength);

		for (i = 0; i < ta->nBioNanoGaps; i++)
		{
			b = ta->gaps + i;
			printBionanpGap(ctx, contigA, contigB, b);
		}

		return;
	}
	b->bionanoGapSize = AdjustedGapLength;

	for (i = 0; i < tb->nBioNanoGaps; i++)
	{
		b = tb->gaps + i;

		if (((b->aBeg == bPartEnd && b->aEnd == bPartBeg) || (b->aBeg == bPartBeg && b->aEnd == bPartEnd)) && ((b->bBeg == aPartBeg && b->bEnd == aPartEnd) || (b->bBeg == aPartEnd && b->bEnd == aPartBeg)))
		{
			break;
		}
	}

	if (i == tb->nBioNanoGaps)
	{
		char *aName = getContigName(ctx, contigA);
		char *bName = getContigName(ctx, contigB);

		printf("[ERROR] - addBionanoGAPInfoToTrimEvidence 2: Cannot find bionano gap feature: Contig %d (%s) and Contig %d (%s): a[%d, %d] b[%d, %d] gapLen %d\n", contigB, aName, contigA, bName, bPartEnd, bPartBeg, aPartEnd, aPartBeg, AdjustedGapLength);

		for (i = 0; i < tb->nBioNanoGaps; i++)
		{
			b = tb->gaps + i;
			printBionanpGap(ctx, contigA, contigB, b);
		}

		return;
	}
	b->bionanoGapSize = AdjustedGapLength;
}

char* getContigName(TrimContext *ctx, int id)
{

	assert(id >= 0);
	assert(id < DB_NREADS(ctx->db));

	int map = 0;
	while (id < ctx->findx[map - 1])
		map -= 1;
	while (id >= ctx->findx[map])
		map += 1;

	return ctx->flist[map];
}

void addBionanoAGPInfoToTrimEvidence(TrimContext *ctx, int contigA, int fromA, int toA, int contigB, int fromB, int toB, int gapLen)
{

	TrimEvidence *t;
	t = find_TrimEvidence(ctx, contigA, contigB);

	int sort = 0;
	if (t == NULL)
	{
		sort = 1;
		t = insert_TrimEvidence(ctx, contigA, contigB);
	}

	assert(t != NULL);

	// add contigA vs contigB
	ensureBionanoGapBuffer(t, 1);
	// check if the same gap feature is already present
	int i;
	BionanoGap *b;
	for (i = 0; i < t->nBioNanoGaps; i++)
	{
		b = t->gaps + i;

		if (intersect(b->aBeg, b->aEnd, fromA, toA) > 0 || intersect(b->bBeg, b->bEnd, fromB, toB) != 0)
		{
			printf("[ERROR] - 1: ambiguous Bioano gap for Contig %d and Contig %d: a[%d, %d] b[%d, %d] gapLen %d, "
					"collides witrh existing gap: Contig %d and Contig %d: a[%d, %d] b[%d, %d] gapLen %d\n", contigA, contigB, fromA, toA, fromB, toB, gapLen, contigA, contigB, b->aBeg, b->aEnd, b->bBeg, b->bEnd, b->agpGapSize);
			exit(1); // todo remove later, for now check if this occurs
		}
	}
	// add gap feature
	b = t->gaps + t->nBioNanoGaps;
	b->aBeg = fromA;
	b->aEnd = toA;
	b->bBeg = fromB;
	b->bEnd = toB;
	b->agpGapSize = gapLen;
	t->nBioNanoGaps++;

	t = find_TrimEvidence(ctx, contigB, contigA);
	if (t == NULL)
	{
		sort = 1;
		t = insert_TrimEvidence(ctx, contigB, contigA);
	}

	assert(t != NULL);

	// add contigB vs contigA
	ensureBionanoGapBuffer(t, 1);
	for (i = 0; i < t->nBioNanoGaps; i++)
	{
		b = t->gaps + i;

		if (intersect(b->aBeg, b->aEnd, toB, fromB) > 0 || intersect(b->bBeg, b->bEnd, toA, fromA) != 0)
		{
			printf("[ERROR] - 2: ambiguous Bioano gap for Contig %d and Contig %d: a[%d, %d] b[%d, %d] gapLen %d, "
					"collides witrh existing gap: Contig %d and Contig %d: a[%d, %d] b[%d, %d] gapLen %d\n", contigB, contigA, toB, fromB, toA, fromA, gapLen, contigB, contigA, b->aBeg, b->aEnd, b->bBeg, b->bEnd, b->agpGapSize);
			exit(1); // todo remove later, for now check if this occurs
		}
	}
	// add gap feature
	b = t->gaps + t->nBioNanoGaps;
	b->aBeg = toB;
	b->aEnd = fromB;
	b->bBeg = toA;
	b->bEnd = fromA;
	b->agpGapSize = gapLen;
	t->nBioNanoGaps++;

	// ensure sort order
	if (sort)
	{
		qsort(ctx->trimEvid, ctx->numTrimEvidence, sizeof(TrimEvidence), TrimEvidence_cmp);
	}
}

int addLASchainInfoToTrimEvidence(TrimContext *ctx, int aread, int bread, int alnLen, int unAlnLen, float erate, int cutPosInA)
{
	TrimEvidence *t = find_TrimEvidence(ctx, aread, bread);

	int result = 1;

	int sort = 0;
	if (t == NULL)
	{
		sort = 1;
		t = insert_TrimEvidence(ctx, aread, bread);
		result = 0; // create new TrimEvidence (i.e. no bionano evidence avaliable)
	}

	assert(t != NULL);

	// add contigA vs contigB
	ensureLASchainBuffer(t, 1);
	int i;
	LASchain *c;
	char *aName = getContigName(ctx, aread);
	char *bName = getContigName(ctx, bread);

	printf("addLASchainInfoToTrimEvidence %d (%s) vs %d (%s), aln %d, unAln: %d, err: %.3f, cut: %d\n", aread, aName, bread, bName, alnLen, unAlnLen, erate, cutPosInA);
	for (i = 0; i < t->nLASchains; i++)
	{
		c = t->chains + i;

		if ((c->trimPos < 0 && cutPosInA < 0) || (c->trimPos > 0 && cutPosInA > 0))
		{
			printf("[ERROR] addLASchainInfoToTrimEvidence: ambiguous contig %d (%s) vs contig %d (%s) overlap present!\n", aread, aName, bread, bName);
			printf("                                       new LASchain evidence: alnLen %d unAlnLen: %d, erate %f, cutPos: %d collides with: previously added LASchain evidence: alnLen %d unAlnLen: %d, erate %f, cutPos: %d\n", alnLen, unAlnLen, erate, cutPosInA, c->alnLen, c->unalignedBases, c->eRate,
					c->trimPos);
			return 2;
		}
	}

	c = t->chains + t->nLASchains;
	c->alnLen = alnLen;
	c->eRate = erate;
	c->unalignedBases = unAlnLen;
	c->trimPos = cutPosInA;

	t->nLASchains++;
	// ensure sort order
	if (sort)
	{
		qsort(ctx->trimEvid, ctx->numTrimEvidence, sizeof(TrimEvidence), TrimEvidence_cmp);
	}
	return result;
}

TrimEvidence*
find_TrimEvidence(TrimContext *ctx, const int contigA, const int contigB)
{
	TrimEvidence target;
	target.contigA = contigA;
	target.contigB = contigB;

	return bsearch(&target, ctx->trimEvid, ctx->numTrimEvidence, sizeof(TrimEvidence), TrimEvidence_cmp);
}

TrimEvidence*
insert_TrimEvidence(TrimContext *ctx, const int contigA, const int contigB)
{
	if (ctx->numTrimEvidence + 3 >= ctx->maxTrimEvidence)
	{
		int i = ctx->maxTrimEvidence * 1.2 + 10;
		ctx->trimEvid = (TrimEvidence*) realloc(ctx->trimEvid, sizeof(TrimEvidence) * i);
		assert(ctx->trimEvid != NULL);
		bzero(ctx->trimEvid + ctx->maxTrimEvidence, sizeof(TrimEvidence) * (i - ctx->maxTrimEvidence));
		ctx->maxTrimEvidence = i;
	}

	assert(ctx->trimEvid != NULL);

	TrimEvidence *result = ctx->trimEvid + ctx->numTrimEvidence;
	ctx->numTrimEvidence++;

	result->contigA = contigA;
	result->contigB = contigB;

	return result;
}

static void trim_pre(PassContext *pctx, TrimContext *tctx)
{
	if (tctx->verbose)
	{
		printf( ANSI_COLOR_GREEN "PASS contig trimming\n" ANSI_COLOR_RESET);

		printf( ANSI_COLOR_RED "OPTIONS\n" ANSI_COLOR_RESET);
		printf( ANSI_COLOR_RED "  verbose %d\n" ANSI_COLOR_RESET, tctx->verbose);
		printf( ANSI_COLOR_RED "  minBionanoGapLen %d\n" ANSI_COLOR_RESET, tctx->minBionanoGapLen);
		printf( ANSI_COLOR_RED "  maxTrimLength %d\n" ANSI_COLOR_RESET, tctx->maxTrimLength);
		printf( ANSI_COLOR_RED "  maxLowCompTrimPerc %d\n" ANSI_COLOR_RESET, tctx->maxLowCompTrimPerc);
		printf( ANSI_COLOR_RED "  trimOffset %d\n" ANSI_COLOR_RESET, tctx->trimOffset);
		printf( ANSI_COLOR_RED "  maxFuzzyBases %d\n" ANSI_COLOR_RESET, tctx->maxFuzzyBases);

		if (tctx->trackDust)
			printf( ANSI_COLOR_RED "  dust Track %s\n" ANSI_COLOR_RESET, tctx->trackDust->name);
		if (tctx->trackTan)
			printf( ANSI_COLOR_RED "  tandem Track %s\n" ANSI_COLOR_RESET, tctx->trackTan->name);
	}

	if (pctx != NULL)
		tctx->twidth = pctx->twidth;

	tctx->maxTrimEvidence = 100;
	tctx->numTrimEvidence = 0;
	tctx->trimEvid = (TrimEvidence*) malloc(sizeof(TrimEvidence) * tctx->maxTrimEvidence);
	assert(tctx->trimEvid != NULL);
	bzero(tctx->trimEvid, sizeof(TrimEvidence) * tctx->maxTrimEvidence);

	int i;
	tctx->trimCoord = (TrimCoordinates*) malloc(sizeof(TrimCoordinates) * DB_NREADS(tctx->db));
	bzero(tctx->trimCoord, sizeof(TrimCoordinates) * DB_NREADS(tctx->db));
	for (i = 0; i < DB_NREADS(tctx->db); i++)
	{
		TrimCoordinates *tc = tctx->trimCoord + i;
		tc->numCoordPairs = 1;
		tc->maxCoordPairs = 1;
		tc->coord = malloc(sizeof(int) * tc->numCoordPairs * 2);
		tc->coord[0] = 0;
		tc->coord[1] = DB_READ_LEN(tctx->db, i);
	}
}

static void trim_post(TrimContext *ctx)
{
	if (ctx->verbose)
	{
		if (ctx->statsTrimmedContigs > 0)
		{
			printf("#trimmed contigs %d\n", ctx->statsTrimmedContigs);
		}

		if (ctx->statsTrimmedBases > 0)
		{
			printf("#trimmed bases: %d\n", ctx->statsTrimmedBases);
		}
		if (ctx->statsNumValidLASchains > 0)
		{
			printf("#valid chains %d\n", ctx->statsNumValidLASchains);
		}
		if (ctx->statsNumInValidLASchains > 0)
		{
			printf("#invalid chains %d\n", ctx->statsNumInValidLASchains);
		}
		if (ctx->statsBionanoTrimmedContigs > 0)
		{
			printf("Bionano Gaps # contigs trimmed: %d; #trimmed bases: %d\n", ctx->statsBionanoTrimmedContigs, ctx->statsBionanoTrimmedBases);
		}
		if (ctx->statsBionanoGapsMissed > 0)
		{
			printf("#not trimmed bionano gaps <= %d: %d\n", ctx->minBionanoGapLen, ctx->statsBionanoGapsMissed);
		}

		if (ctx->statsBionanoGapsAll > 0)
		{
			printf("binano agp file: #ALL gaps: %d; #gaps (<= %d bases): %d, #gaps (<= %d bases at ContigBreaks): %d \n", ctx->statsBionanoGapsAll, ctx->minBionanoGapLen, ctx->statsBionanoGapsLtMinThresh, ctx->minBionanoGapLen, ctx->statsBionanoGapsLtMinThreshContigBreak);
		}
	}
}

int getTrimPositionsFromLAS(TrimContext *ctx, Overlap *ovl, int pointA, int *cutA, int *cutB)
{
	int abeg = ovl->path.abpos;
	int aend = ovl->path.aepos;

	int bbeg = ovl->path.bbpos;
	int bend = ovl->path.bepos;

	int twidth = ctx->twidth;

	if (pointA < abeg || pointA > aend)
		return 1;

	//printf("getTrimPositions %d x %d, a(%d, %d) %c b(%d, %d) pointA: %d\n", ovl->aread, ovl->bread, abeg, aend, (ovl->flags & OVL_COMP) ? 'C' : 'N', bbeg, bend, pointA);

	int dist = pointA - abeg;
	int apos, bpos;
	if (ovl->path.tlen)
	{
		ovl_trace *trace = ovl->path.trace;
		apos = abeg;
		bpos = ovl->path.bbpos;

		int j = 0;
		while (j < ovl->path.tlen)
		{
			//      printf("apos %6d, bpos %6d, oldDist %6d, newDist %6d\n", apos, bpos, dist, abs(pointA-((apos / twidth + 1) * twidth)));
			if (dist < abs(pointA - ((apos / twidth + 1) * twidth)))
				break;
			apos = (apos / twidth + 1) * twidth;
			bpos += trace[j + 1];
			//printf("apos %6d, bpos %6d\n", apos, bpos);
			dist = abs(pointA - apos);
			j += 2;
		}
	}
	else
	{
		apos = pointA;
		bpos = bbeg + dist;
	}

	//  printf("apos: %d, bpos: %d\n", apos, bpos);

	*cutA = apos;
	*cutB = bpos;

	if (*cutA < abeg || *cutA > aend)
		return 1;

	if (*cutB < bbeg || *cutB > bend)
		return 1;

	//printf("final range: %d, %d\n", *cutA, *cutB);
	return 0;
}

int analyzeContigOverlaps(TrimContext *ctx, Overlap *ovl, int novl)
{
	int i;

	int aLen = DB_READ_LEN(ctx->db, ovl->aread);
	int bLen = DB_READ_LEN(ctx->db, ovl->bread);

	char *aName = getContigName(ctx, ovl->aread);
	char *bName = getContigName(ctx, ovl->bread);

	// assumption: input overlaps must be chained with LAfilterChains !!!
	// one chain at the end and one chain at the beginning of a contig are possible! But BOT with the same contigA and contigB

	// sanity check
	Overlap *o1 = ovl;

	int cutA = -1;
	int cutB = -1;

	if (novl == 1)
	{
		// check if overlaps is valid
		// 1. check if overlap is a valid chain!!
		if ((o1->path.abpos > ctx->maxFuzzyBases && o1->path.aepos < aLen - ctx->maxFuzzyBases) || (o1->path.bbpos > ctx->maxFuzzyBases && o1->path.bepos < bLen - ctx->maxFuzzyBases))
		{
			if (ctx->verbose)
			{
				printf("[WARNGING] fuzzy base check failed! Ignore invalid chain [%d, %d] a[%d,%d] %c b[%d,%d]!\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o1->path.bepos);
			}
			ctx->statsNumInValidLASchains++;
			ctx->statsNumInValidLASchainOverlaps += novl;
			return 1;
		}
		// 2. check if one contig is contained within another one!!
		else if ((o1->path.abpos <= aLen / 2 && (aLen - o1->path.aepos) <= aLen / 2) || (o1->path.bbpos <= bLen / 2 && (bLen - o1->path.bepos) <= bLen / 2))
		{
			printf("[WARNGING] Containment found! Ignore invalid chain [%d, %d] a[%d,%d] %c b[%d,%d]!\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o1->path.bepos);
			ctx->statsNumInValidLASchains++;
			ctx->statsNumInValidLASchainOverlaps += novl;
			return 1;
		}

		int pointA = o1->path.abpos + (o1->path.aepos - o1->path.abpos) / 2;

		if (getTrimPositionsFromLAS(ctx, o1, pointA, &cutA, &cutB))
		{
			printf("Unable to get cutPosition for OVL [%d,%d] a[%d,%d] %c b[%d,%d] and pointA: %d\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o1->path.bepos, pointA);
			ctx->statsNumInValidLASchains++;
			ctx->statsNumInValidLASchainOverlaps += novl;
			return 1;
		}

		assert((cutA - ctx->trimOffset > 0) && (cutA + ctx->trimOffset < aLen));
		assert((cutB - ctx->trimOffset > 0) && (cutB + ctx->trimOffset < bLen));

		float erate = (200. * ovl->path.diffs) / ((ovl->path.aepos - ovl->path.abpos) + (ovl->path.bepos - ovl->path.bbpos));
		int resA = 2;

		// set cut position of contig_A
		if (o1->path.abpos < aLen - o1->path.aepos) // trim off contig at begin
		{
			resA = addLASchainInfoToTrimEvidence(ctx, ovl->aread, ovl->bread, ovl->path.aepos - ovl->path.abpos, ovl->path.abpos, erate, -(cutA + ctx->trimOffset));
		}
		else if (o1->path.abpos > aLen - o1->path.aepos) // trim off contig at end
		{
			resA = addLASchainInfoToTrimEvidence(ctx, ovl->aread, ovl->bread, ovl->path.aepos - ovl->path.abpos, aLen - ovl->path.aepos, erate, cutA - ctx->trimOffset);
		}
		else // containment
		{
			printf("[WARNGING] Containment found! Ignore invalid chain [%d (%s), %d (%s)]  a[%d,%d] %c b[%d,%d]!\n", o1->aread, aName, o1->bread, bName, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o1->path.bepos);
			exit(1);
		}

		int resB = 2;
		// set cut position of contig_B
		if (o1->path.bbpos < bLen - o1->path.bepos) // trim off contig at begin
		{
			if (o1->flags & OVL_COMP)
			{
				resB = addLASchainInfoToTrimEvidence(ctx, ovl->bread, ovl->aread, ovl->path.bepos - ovl->path.bbpos, ovl->path.bbpos, erate, bLen - (cutB + ctx->trimOffset));
			}
			else
			{
				resB = addLASchainInfoToTrimEvidence(ctx, ovl->bread, ovl->aread, ovl->path.bepos - ovl->path.bbpos, ovl->path.bbpos, erate, -(cutB + ctx->trimOffset));
			}
		}
		else if (o1->path.bbpos > bLen - o1->path.bepos) // trim off contig at end
		{
			if (o1->flags & OVL_COMP)
			{

				resB = addLASchainInfoToTrimEvidence(ctx, ovl->bread, ovl->aread, ovl->path.bepos - ovl->path.bbpos, bLen - ovl->path.bepos, erate, -(bLen - (cutB - ctx->trimOffset)));
			}
			else
			{
				resB = addLASchainInfoToTrimEvidence(ctx, ovl->bread, ovl->aread, ovl->path.bepos - ovl->path.bbpos, bLen - ovl->path.bepos, erate, cutB - ctx->trimOffset);
			}
		}

		assert(resA == resB);

		// check return value of res and update stats
		if (resA == 2)
		{
			ctx->statsNumValidLASchains++;
			ctx->statsNumValidLASchainOverlaps += novl;
			ctx->statsNumDuplicatedChains++;
			ctx->statsNumDuplicatedChains += novl;
		}
		else if (resA == 1)
		{
			ctx->statsNumValidLASchains++;
			ctx->statsNumValidLASchainOverlaps += novl;
			ctx->statsNumLASChainsWithBionanoSupport++;
		}
		else
		{
			ctx->statsNumValidLASchains++;
			ctx->statsNumValidLASchainOverlaps += novl;
			ctx->statsNumLASChainsWithoutBionanoSupport++;
		}
	}
	else
	{
		int validChain = 1;
		float avgErate = (200. * o1->path.diffs) / ((o1->path.aepos - o1->path.abpos) + (o1->path.bepos - o1->path.bbpos));
		int alignedBasesInA = o1->path.aepos - o1->path.abpos;
		int alignedBasesInB = o1->path.bepos - o1->path.bbpos;
		int unalignedBasesInA = MIN(o1->path.abpos, aLen - ovl[novl - 1].path.aepos);
		int unalignedBasesInB = MIN(o1->path.bbpos, bLen - ovl[novl - 1].path.bepos);

		int resA = 2;
		int resB = 2;
		Overlap *o2;
		// first: sanity check for LAS chain
		for (i = 1; i < novl; i++)
		{
			o2 = ovl + i;
			if (abs(o1->path.aepos - o2->path.abpos) > ctx->maxFuzzyBases || ((o1->flags & OVL_COMP) != (o2->flags & OVL_COMP)))
			{
				validChain = 0;
				break;
			}

			alignedBasesInA += (o2->path.aepos - o2->path.abpos);
			alignedBasesInB += (o2->path.bepos - o2->path.bbpos);

			if (o1->path.aepos > o2->path.abpos)
			{
				alignedBasesInA -= (o1->path.aepos - o2->path.abpos);
			}
			else
			{
				unalignedBasesInA += (o2->path.abpos - o1->path.aepos);
			}
			if (o1->path.bepos > o2->path.bbpos)
			{
				alignedBasesInB -= (o1->path.bepos - o2->path.bbpos);
			}
			else
			{
				unalignedBasesInB += (o2->path.bbpos - o1->path.bepos);
			}
			avgErate += (200. * o2->path.diffs) / ((o2->path.aepos - o2->path.abpos) + (o2->path.bepos - o2->path.bbpos));

			o1 = o2;
		}

		// reset overlap pointer o1 and o2 to first and last overlap of LASchain respectively
		o1 = ovl;
		o2 = ovl + (novl - 1);
		avgErate /= novl;

		// check for containment
		if (validChain)
		{
			if ((o1->path.abpos <= aLen / 2 && (aLen - o2->path.aepos) <= aLen / 2) || (o1->path.bbpos <= bLen / 2 && (bLen - o2->path.bepos) <= bLen / 2))
			{
				validChain = 0;
			}

			if (!validChain)
			{
				// special case for small contigs, with very good alignments allow up to 2/3 overlap and 1/3 overhang
				if (avgErate < 2.0 && unalignedBasesInA < ctx->maxFuzzyBases && unalignedBasesInB < ctx->maxFuzzyBases && (o1->path.abpos >= aLen / 3 || (aLen - o2->path.aepos) >= aLen / 3) && (o1->path.bbpos >= bLen / 3 || (bLen - o2->path.bepos) >= bLen / 3))
				{
					validChain = 1;
					printf("[WARNING] Containment found BUT IS KEPT %d in %d (e %.2f, aln: %d unaln: %d)! Ignore invalid chain [%d (%s), %d (%s)]  a[%d,%d] %c b[%d,%d]!\n", o1->aread, o1->bread, avgErate, alignedBasesInA, unalignedBasesInA, o1->aread, aName, o1->bread, bName, o1->path.abpos, o2->path.aepos,
							(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o2->path.bepos);
				}
				else
				{
					printf("[WARNING] Containment found %d in %d (e %.2f, aln: %d unaln: %d)! Ignore invalid chain [%d (%s), %d (%s)]  a[%d,%d] %c b[%d,%d]!\n", o1->aread, o1->bread, avgErate, alignedBasesInA, unalignedBasesInA, o1->aread, aName, o1->bread, bName, o1->path.abpos, o2->path.aepos,
							(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o2->path.bepos);
				}
			}
		}

		if (!validChain)
		{
			ctx->statsNumInValidLASchains++;
			ctx->statsNumInValidLASchainOverlaps += novl;

			if (ctx->verbose)
			{
				printf("INVALID chain: %d (%s) vs %d (%s)\n", ovl->aread, aName, ovl->bread, bName);

				for (i = 0; i < novl; i++)
				{
					printf("   a[%d,%d] %c b[%d,%d]\n", ovl[i].path.abpos, ovl[i].path.aepos, (ovl[i].flags & OVL_COMP) ? 'c' : 'n', ovl[i].path.bbpos, ovl[i].path.bepos);
				}
			}
			return 1;
		}

		// set cut position of contig_A
		if (o1->path.abpos < aLen - o2->path.aepos) // trim off contig at begin
		{
			if (o2->path.abpos >= o1->path.aepos)
			{
				cutA = o2->path.abpos + ctx->trimOffset;
			}
			else
			{
				cutA = o1->path.aepos + ctx->trimOffset;
			}

			resA = addLASchainInfoToTrimEvidence(ctx, o1->aread, o1->bread, alignedBasesInA, unalignedBasesInA, avgErate, -(cutA));
		}
		else if (o1->path.abpos > aLen - o2->path.aepos) // trim off contig at end
		{
			if (o2->path.abpos >= o1->path.aepos)
			{
				cutA = o1->path.aepos - ctx->trimOffset;
			}
			else
			{
				cutA = o2->path.abpos - ctx->trimOffset;
			}

			resA = addLASchainInfoToTrimEvidence(ctx, o1->aread, o1->bread, alignedBasesInA, unalignedBasesInA, avgErate, cutA);
		}
		else // containment
		{
			printf("Contained overlap (c1): [%d (%s),%d (%s)] a[%d,%d] %c b[%d,%d]\n", o1->aread, aName, o1->bread, bName, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o1->path.bepos);
			ctx->statsNumInValidLASchains++;
			ctx->statsNumInValidLASchainOverlaps += novl;
			return 1;
		}
		// set cut position of contig_B

		if (o1->path.bbpos < bLen - o2->path.bepos) // trim off contig at begin
		{
			if (o2->path.bbpos >= o1->path.bepos)
			{
				cutB = o2->path.bbpos + ctx->trimOffset;
				if (o1->flags & OVL_COMP)
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, bLen - cutB);
				}
				else
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, -(cutB));
				}
			}
			else
			{
				cutB = o1->path.bepos + ctx->trimOffset;
				if (o1->flags & OVL_COMP)
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, bLen - cutB);
				}
				else
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, -(cutB));
				}
			}
		}
		else if (o1->path.bbpos > bLen - o2->path.bepos) // trim off contig at end
		{
			if (o2->path.bbpos >= o1->path.bepos)
			{
				cutB = o1->path.bepos - ctx->trimOffset;
				if (o1->flags & OVL_COMP)
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, -(bLen - cutB));
				}
				else
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, cutB);
				}
			}
			else
			{
				cutB = o2->path.bbpos - ctx->trimOffset;
				if (o1->flags & OVL_COMP)
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, -(bLen - cutB));
				}
				else
				{
					resB = addLASchainInfoToTrimEvidence(ctx, o1->bread, o1->aread, alignedBasesInB, unalignedBasesInB, avgErate, cutB);
				}

			}
		}
		else // containment
		{
			printf("Contained overlap (c2): [%d (%s),%d (%s)] a[%d,%d] %c b[%d,%d]\n", o1->aread, aName, o1->bread, bName, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos, o1->path.bepos);
			ctx->statsNumInValidLASchains++;
			ctx->statsNumInValidLASchainOverlaps += novl;
			return 1;
		}

		assert(resA == resB);

		// check return value of res and update stats
		if (resA == 2)
		{
			ctx->statsNumValidLASchains++;
			ctx->statsNumValidLASchainOverlaps += novl;
			ctx->statsNumDuplicatedChains++;
			ctx->statsNumDuplicatedChains += novl;
		}
		else if (resA == 1)
		{
			ctx->statsNumValidLASchains++;
			ctx->statsNumValidLASchainOverlaps += novl;
			ctx->statsNumLASChainsWithBionanoSupport++;
		}
		else
		{
			ctx->statsNumValidLASchains++;
			ctx->statsNumValidLASchainOverlaps += novl;
			ctx->statsNumLASChainsWithoutBionanoSupport++;
		}

	}
	return 0;
}

static int trim_handler(void *_ctx, Overlap *ovl, int novl)
{
	TrimContext *ctx = (TrimContext*) _ctx;

	// analyze overlaps and find contig trim position
	analyzeContigOverlaps(ctx, ovl, novl);

	return 1;
}

int getMaskedBases(TrimContext *ctx, HITS_TRACK *t, int contigID, int beg, int end)
{
#ifdef DEBUG_MASKING
	printf("call getMaskedBases on track %s, contigID: %d, in: [%d, %d]\n", t->name, contigID, beg, end);
#endif
	if (t == NULL)
	{
#ifdef DEBUG_MASKING        
		printf(" --> masked bases 0 (track is Null)\n");
#endif
		return 0;
	}

	track_anno *mask_anno = t->anno;
	track_data *mask_data = t->data;

	if (contigID < 0 || contigID >= DB_NREADS(ctx->db))
	{
		fprintf(stderr, "[ERROR] - getMaskedBases contigID: %d out of bounds [0, %d]\n", contigID, DB_NREADS(ctx->db) - 1);
		fflush(stderr);
		exit(1);
	}

	track_anno rb, re;

	int maskBases = 0;
	int rBeg, rEnd;

	// repeat bases in a-read
	rb = mask_anno[contigID] / sizeof(track_data);
	re = mask_anno[contigID + 1] / sizeof(track_data);

	while (rb < re)
	{
		rBeg = mask_data[rb];
		rEnd = mask_data[rb + 1];

		maskBases += intersect(beg, end, rBeg, rEnd);

#ifdef DEBUG_MASKING2
        printf("     repInterval: [%d, %d] intersection with [%d, %d] is %d. cum sum %d\n", rBeg, rEnd, beg, end, intersect(beg, end, rBeg, rEnd), maskBases);
#endif

		rb += 2;
	}

#ifdef DEBUG_MASKING
	printf(" --> masked bases %d\n", maskBases);
#endif

	return maskBases;
}

char* trimwhitespace(char *str)
{
	//printf("trimwhitespace: %s\n", str);
	char *end;

	// Trim leading space
	while (isspace(*str))
		str++;

	if (*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while (end > str && isspace(*end))
		end--;

	// Write new null terminator
	*(end + 1) = '\0';

	return str;
}

int getDBcontigID(TrimContext *ctx, char *contigName, int *from, int *to)
{
//	printf("getDBcontigID(%s)\n", contigName);
	int i;
	for (i = 0; i < ctx->nfiles; i++)
	{
		if (strcmp(ctx->flist[i], contigName) == 0)
		{
			return i;
		}
	}

	int cNameLen;
	char *pchrf, *pchrl;

	pchrf = strstr(contigName, "_subseq_");

	if (pchrf == NULL)
		return -1;

	pchrl = strstr(pchrf, ":");

	if (pchrl == NULL)
		return -1;

	int agpCNameLen = pchrf - contigName;

	// printf("contigNameLen from %s is %d\n", contigName, agpCNameLen);

	for (i = 0; i < ctx->nfiles; i++)
	{
		cNameLen = strlen(ctx->flist[i]);

		if (agpCNameLen == cNameLen && strncmp(ctx->flist[i], contigName, cNameLen) == 0)
		{

			*pchrl = '\0';
			*from = strtol(pchrf + 8, NULL, 10);
			*pchrl = ':';
			*to = strtol(pchrl + 1, NULL, 10);

			//	printf("found: from: %d and to: %d\n", *from, *to);
			// sanity checks
			if (*from < 1 || *from > *to || *from > DB_READ_LEN(ctx->db, i) || *to > DB_READ_LEN(ctx->db, i))
			{
				*from = -1;
				*to = -1;
				return -1;
			}

			return i;
		}
	}

	return -1;
}

void parseBionanoGAPfile(TrimContext *ctx, char *pathInBionanoGAP)
{
	FILE *fileInBionanoGaps = NULL;

	if ((fileInBionanoGaps = fopen(pathInBionanoGAP, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] could not open %s\n", pathInBionanoGAP);
		exit(1);
	}

	char NGSId1[MAX_NAME];
	char NGSId2[MAX_NAME];
	int SuperScaffoldId;
	int XmapGapLength;
	int AdjustedGapLength;
	float NGSLength1;
	float NGSLength2;

	char *line = NULL;
	size_t maxline = 0;

	int nline = 0;
	int len;

	int r;
	int contigA = -1;
	int contigB = -1;
	int numInvalidLines = 0;
	printf("parseBionanoGapfile: %s\n", pathInBionanoGAP);
	while ((len = getline(&line, &maxline, fileInBionanoGaps)) > 0)
	{
		nline++;

		char *tline = trimwhitespace(line);

		if (tline[0] == '#')
			continue;

		r = sscanf(tline, "%s\t%s\t%d\t%d\t%d\t%f\t%f\n", NGSId1, NGSId2, &SuperScaffoldId, &XmapGapLength, &AdjustedGapLength, &NGSLength1, &NGSLength2);

		if (r != 7)
		{
			fprintf(stderr, "[ERROR] invalid bionano GAP file format %s. Expecting 7 columns, BUT parsed %d columns in line %d\n", pathInBionanoGAP, r, nline);
			exit(1);
		}

//		printf("line %d: %s\n", nline, tline);

		// try to match contig name with with DB contig ID

		int aPartBeg = -1;
		int aPartEnd = -1;

		int bPartBeg = -1;
		int bPartEnd = -1;

		contigA = getDBcontigID(ctx, NGSId1, &aPartBeg, &aPartEnd);
		contigB = getDBcontigID(ctx, NGSId2, &bPartBeg, &bPartEnd);

		if (contigA < 0 || contigB < 0)
		{
			printf("[WARNING] Could not match GAP contig names: %s and/or %s in current db! Ignore line %d in GAP file %s.\n", NGSId1, NGSId2, nline, pathInBionanoGAP);
			numInvalidLines++;
			continue;
		}

		// contig orientation is unknown !!!!
		if (aPartBeg < 0 || aPartEnd < 0)
		{
			aPartBeg = 1;
			aPartEnd = DB_READ_LEN(ctx->db, contigA);
		}

		if (bPartBeg < 0 || bPartEnd < 0)
		{
			bPartBeg = 1;
			bPartEnd = DB_READ_LEN(ctx->db, contigB);
		}

		assert(aPartBeg < aPartEnd);
		assert(bPartBeg < bPartEnd);

		// ignore lines were contigA equals contigB. Why do they exist?
		if (contigA == contigB)
		{
			continue;
		}

		addBionanoGAPInfoToTrimEvidence(ctx, contigA, aPartBeg, aPartEnd, contigB, bPartBeg, bPartEnd, AdjustedGapLength);
	}

	int negativeGaps = 0;
	int i, j;
	for (i = 0; i < ctx->numTrimEvidence; i++)
	{
		TrimEvidence *t = ctx->trimEvid + i;

		if (t->contigA > t->contigB)
			continue;

		for (j = 0; j < t->nBioNanoGaps; j++)
		{
			BionanoGap *b = t->gaps + j;
			if (b->bionanoGapSize < 0)
				negativeGaps++;
		}
	}
	printf("[INFO]  Number of invalid lines: %d (either format issues, or AGP contig names could not be matched to DB contig names.)\n", numInvalidLines);
	printf("[INFO]  #Bionano gaps < 0: %10d\n", negativeGaps);

	free(line);
	fclose(fileInBionanoGaps);
}

void printBionanpGap(TrimContext *ctx, int contigA, int contigB, BionanoGap *g)
{
	assert(g != NULL);

	char *aName = getContigName(ctx, contigA);
	char *bName = getContigName(ctx, contigB);

	printf("Bionano Gap: %d(%s)[%d,%d]-------GAP[%d, %d]------%d(%s)[%d,%d]\n", contigA, aName, g->aBeg, g->aEnd, g->agpGapSize, g->bionanoGapSize, contigB, bName, g->bBeg, g->bEnd);
}

void printLASchain(TrimContext *ctx, int contigA, int contigB, LASchain *c)
{
	assert(c != NULL);

	char *aName = getContigName(ctx, contigA);
	char *bName = getContigName(ctx, contigB);

	printf("LASchain: %3d (%s) vs %3d (%s): #alnBases: %6d #unalnBases: %6d eRate: %5.2f trimPosOfA: %6d\n", contigA, aName, contigB, bName, c->alnLen, c->unalignedBases, c->eRate, c->trimPos);
}

void parseBionanoAGPfile(TrimContext *ctx, char *pathInBionanoAGP)
{
	FILE *fileInBionanoGaps = NULL;

	if ((fileInBionanoGaps = fopen(pathInBionanoAGP, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] could not open %s\n", pathInBionanoAGP);
		exit(1);
	}

	char Prev_Obj_Name[MAX_NAME];
	char Obj_Name[MAX_NAME];
	int Obj_Start;
	int Obj_End;
	int PartNum;
	char Compnt_Type;
	char CompntId_GapLength[MAX_NAME];
	char CompntStart_GapType[MAX_NAME];
	char CompntEnd_Linkage[MAX_NAME];
	char Orientation_LinkageEvidence[MAX_NAME];

	char *line = NULL;
	size_t maxline = 0;

	int nline = 0;
	int len;

	int r;
	int contigA = -1;
	int oriA = 0;
	int fromA, toA;
	char contigNameA[MAX_NAME];
	int gapLen = -1;

	int contigB = -1;
	int oriB = 0;
	int fromB, toB;
	char contigNameB[MAX_NAME];

	printf("[INFO] parseBionanoAGPfile: %s\n", pathInBionanoAGP);

	int numInvalidLines = 0;

	int prevContigID = -1;
	int Prev_PartNum = -1;

	while ((len = getline(&line, &maxline, fileInBionanoGaps)) > 0)
	{
		nline++;

		char *tline = trimwhitespace(line);

		if (tline[0] == '#')
			continue;

		r = sscanf(tline, "%s\t%d\t%d\t%d\t%c\t%s\t%s\t%s\t%s\n", Obj_Name, &Obj_Start, &Obj_End, &PartNum, &Compnt_Type, CompntId_GapLength, CompntStart_GapType, CompntEnd_Linkage, Orientation_LinkageEvidence);

		if (r != 9)
		{
			fprintf(stderr, "[ERROR] invalid AGP file format %s. Expecting 9 columns, BUT parsed %d columns in line %d\n", pathInBionanoAGP, r, nline);
			exit(1);
		}

		if (ctx->verbose > 2)
			printf("line %d: %s\n", nline, tline);

		// try to match contig name with with DB contig ID
		if (Compnt_Type == 'W')
		{
			int from = -1;
			int to = -1;

			if (strcmp(Prev_Obj_Name, Obj_Name) != 0)
			{
				if(Prev_PartNum == 1)
				{
					printf(" SINGLETON ContigA %d prev: %s cur: %s\n", prevContigID, Prev_Obj_Name, Obj_Name);
				}

				strcpy(Prev_Obj_Name, Obj_Name);
				contigA = getDBcontigID(ctx, CompntId_GapLength, &from, &to);
				strcpy(contigNameA, CompntId_GapLength);
				if (contigA < 0)
				{
					printf("[WARNING] Could not match agp contig name: %s in current db! Ignore AGP file.\n", CompntId_GapLength);
					Prev_Obj_Name[0] = '\0'; // restart from scratch
					numInvalidLines++;
					continue;
				}

				if (strcmp(Orientation_LinkageEvidence, "+") == 0)
				{
					oriA = 1;
				}
				else if (strcmp(Orientation_LinkageEvidence, "-") == 0)
				{
					oriA = -1;
				}
				else
				{
					fprintf(stderr, "[ERROR] invalid AGP file format %s. Unknown orientation %s in line %d\n", pathInBionanoAGP, Orientation_LinkageEvidence, nline);
					Prev_Obj_Name[0] = '\0'; // restart from scratch
					numInvalidLines++;
					continue;
				}

				if (from == -1)
				{
					from = 1;
					to = DB_READ_LEN(ctx->db, contigA);
				}

				if (oriA < 0)
				{
					fromA = to;
					toA = from;
				}
				else
				{
					fromA = from;
					toA = to;
				}

				Prev_PartNum = PartNum;
				prevContigID = contigA;
			}
			else
			{
				Prev_PartNum = PartNum;
				contigB = getDBcontigID(ctx, CompntId_GapLength, &from, &to);
				strcpy(contigNameB, CompntId_GapLength);

				if (contigB < 0)
				{
					printf("[WARNING] Could not match agp contig name: %s in current db! Ignore AGP file.\n", CompntId_GapLength);
					Prev_Obj_Name[0] = '\0'; // restart from scratch
					numInvalidLines++;
					continue;
				}
				if (strcmp(Orientation_LinkageEvidence, "+") == 0)
				{
					oriB = 1;
				}
				else if (strcmp(Orientation_LinkageEvidence, "-") == 0)
				{
					oriB = -1;
				}
				else
				{
					fprintf(stderr, "[ERROR] invalid AGP file format %s. Unknown orientation %s in line %d\n", pathInBionanoAGP, Orientation_LinkageEvidence, nline);
					Prev_Obj_Name[0] = '\0'; // restart from scratch
					numInvalidLines++;
					continue;
				}

				if (from == -1)
				{
					from = 1;
					to = DB_READ_LEN(ctx->db, contigB);
				}

				if (oriB < 0)
				{
					fromB = to;
					toB = from;
				}
				else
				{
					fromB = from;
					toB = to;
				}

				printf("2: fromA: %d, toA: %d\n", fromA, toA);
				assert(gapLen > -1);

				// add trim evidence symmetrically: i.e. contigA-gap-contigB and contigB-gap-contigA
				addBionanoAGPInfoToTrimEvidence(ctx, contigA, fromA, toA, contigB, fromB, toB, gapLen);

				contigA = contigB;
				strcpy(contigNameA, contigNameB);
				oriA = oriB;
				fromA = fromB;
				toA = toB;
				// reset gap size
				gapLen = -1;

				//fprintf(stdout, "Assign B to A: ContigA[%d,%s,%d,%d,%d] - GAP [%d] - ContigB[%d,%s,%d,%d,%d]\n", contigA, contigNameA, oriA, fromA, toA, gapLen, contigB, contigNameB, oriB, fromB, toB);
			}

		}
		else if (Compnt_Type == 'N')
		{
			//printf("found gap\n");
			gapLen = strtol(CompntId_GapLength, NULL, 10);
			if (gapLen < 1)
			{
				//fprintf(stderr, "[ERROR] invalid AGP file format %s. Negative gap length: %s in line %d\n", pathInBionanoAGP, CompntId_GapLength, nline);
			}
		}
	}

	printf("[INFO]  Number of invalid lines: %d (either format issues, or AGP contig names could not be matched to DB contig names.)\n", numInvalidLines);

	int numGaps = 0;
	int numGapsSmallerThreshold = 0;
	int numContigBreaksPartOfAGap = 0;
	int numContigBreaksNotClosable = 0;
	int i, j;
	for (i = 0; i < ctx->numTrimEvidence; i++)
	{
		TrimEvidence *t = ctx->trimEvid + i;

		if (t->contigA > t->contigB)
			continue;

		int aLen = DB_READ_LEN(ctx->db, t->contigA);
		int bLen = DB_READ_LEN(ctx->db, t->contigB);

		for (j = 0; j < t->nBioNanoGaps; j++)
		{
			BionanoGap *b = t->gaps + j;

			if (ctx->verbose)
				printBionanpGap(ctx, t->contigA, t->contigB, b);

			numGaps++;
			if (b->agpGapSize <= ctx->minBionanoGapLen)
			{
				numGapsSmallerThreshold++;
			}
			if ((b->aBeg != 1 && b->aBeg != aLen) || (b->aEnd != 1 && b->aEnd != aLen))
			{
				numContigBreaksPartOfAGap++;
				if (ctx->verbose)
				{
					printf("BreaksPartOfAGap: if ((b->aBeg != 1 && b->aBeg != aLen) || (b->aEnd != 1 && b->aEnd != aLen))\n");
					printBionanpGap(ctx, t->contigA, t->contigB, b);
				}
			}
			else if ((b->bBeg != 1 && b->bBeg != bLen) || (b->bEnd != 1 && b->bEnd != bLen))
			{
				numContigBreaksPartOfAGap++;
				if (ctx->verbose)
				{
					printf("BreaksPartOfAGap: if ((b->bBeg != 1 && b->bBeg != bLen) || (b->bEnd != 1 && b->bEnd != bLen))\n");
					printBionanpGap(ctx, t->contigA, t->contigB, b);
				}
			}
			if ((b->aEnd != 1 && b->aEnd != aLen) || (b->bBeg != 1 && b->bBeg != bLen))
			{
				numContigBreaksNotClosable++;
				if (ctx->verbose)
				{
					printf("NotClosable: if((b->aEnd != 1 && b->aEnd != aLen) || (b->bBeg != 1 && b->bBeg != bLen))\n");
					printBionanpGap(ctx, t->contigA, t->contigB, b);
				}
			}
		}
	}

	ctx->statsBionanoGapsAll = numGaps;
	ctx->statsBionanoGapsLtMinThresh = numGapsSmallerThreshold;
	ctx->statsBionanoContigBreaksPartOfAGap = numContigBreaksPartOfAGap;
	ctx->statsBionanoContigBreaksNotClosable = numContigBreaksNotClosable;

	printf("[INFO]  #BionanoGaps: %15d\n", numGaps);
	printf("[INFO]  #BionanoGaps (<= %d): %7d\n", ctx->minBionanoGapLen, numGapsSmallerThreshold);
	printf("[INFO]  #ContigBreaksPartOfAGap: %4d\n", numContigBreaksPartOfAGap);
	printf("[INFO]  #ContigBreaksNotClosable: %3d\n", numContigBreaksNotClosable);

	free(line);
	fclose(fileInBionanoGaps);
}

void getDBFastaHeader(TrimContext *ctx, char *fullDBPath)
{
	char *pwd, *root;
	FILE *dstub;
	int i;

	root = Root(fullDBPath, ".db");
	pwd = PathTo(fullDBPath);
	if (ctx->db->part > 0)
	{
		fprintf(stderr, "[ERROR] - CTtrim can not work on blocks!");
		exit(1);
	}
	dstub = Fopen(Catenate(pwd, "/", root, ".db"), "r");
	if (dstub == NULL)
	{
		fprintf(stderr, "[ERROR] - Cannot open database file: %s\n", Catenate(pwd, "/", root, ".db"));
		exit(1);
	}
	free(pwd);
	free(root);

	if (fscanf(dstub, DB_NFILE, &(ctx->nfiles)) != 1)
	{
		fclose(dstub);
		fprintf(stderr, "[ERROR] - Cannot read files line '%s' in database file: %s\n",
		DB_NFILE, Catenate(pwd, "/", root, ".db"));
		exit(1);
	}

	ctx->flist = (char**) Malloc(sizeof(char*) * ctx->nfiles, "Allocating file list");
	ctx->hlist = (char**) Malloc(sizeof(char*) * ctx->nfiles, "Allocating header list");
	ctx->findx = (int*) Malloc(sizeof(int*) * (ctx->nfiles + 1), "Allocating file index");
	if (ctx->flist == NULL || ctx->findx == NULL || ctx->hlist == NULL)
	{
		fclose(dstub);
		fprintf(stderr, "[ERROR] - Cannot allocate file name and file index buffers!");
		exit(1);
	}

	ctx->findx += 1;
	ctx->findx[-1] = 0;

	for (i = 0; i < ctx->nfiles; i++)
	{
		char headername[MAX_NAME], filename[MAX_NAME];

		if (fscanf(dstub, DB_FDATA, ctx->findx + i, filename, headername) != 3)
		{
			fclose(dstub);
			fprintf(stderr, "[ERROR] - Cannot read %d-th fasta entry in database file %s\n", i + 1, Catenate(pwd, "/", root, ".db"));
			exit(1);
		}

		if ((ctx->flist[i] = Strdup(filename, "Adding to file list")) == NULL)
			exit(1);

		if ((ctx->hlist[i] = Strdup(headername, "Adding to file list")) == NULL)
			exit(1);
	}

	fclose(dstub);
}

//void trim_contigs(TrimContext *ctx)
/*{
 // open file handler
 FILE *trimmedContigsAll = NULL;
 FILE *purgedContigsAll = NULL;
 FILE *statsContigsAll = NULL;

 FILE *trimmedContigsNoTandem = NULL;
 FILE *purgedContigsNoTandem = NULL;
 FILE *statsContigsNoTandem = NULL;

 char *fout = malloc(strlen(ctx->fileOutPattern) + 50);
 assert(fout != NULL);

 sprintf(fout, "%s.trimmedContigs.fasta", ctx->fileOutPattern);
 printf("create file: %s\n", fout);
 if ((trimmedContigsAll = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }
 sprintf(fout, "%s.purgedContigs.fasta", ctx->fileOutPattern);
 if ((purgedContigsAll = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }
 sprintf(fout, "%s.trimmedContigs.stats", ctx->fileOutPattern);
 if ((statsContigsAll = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }

 sprintf(fout, "%s.ignoreTANtrimmedContigs.fasta", ctx->fileOutPattern);
 printf("create file: %s\n", fout);
 if ((trimmedContigsNoTandem = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }
 sprintf(fout, "%s.ignoreTANpurgedContigs.fasta", ctx->fileOutPattern);
 if ((purgedContigsNoTandem = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }
 sprintf(fout, "%s.ignoreTANtrimmedContigs.stats", ctx->fileOutPattern);
 if ((statsContigsNoTandem = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }

 fprintf(statsContigsAll, "#ContigID\tContigName\tnewContigLength\ttrimBegin\ttrimEnd\tcomments\n");
 fprintf(statsContigsNoTandem, "#ContigID\tContigName\tnewContigLength\ttrimBegin\ttrimEnd\tcomments\n");

 // debug report trim positions
 int nContigs = DB_NREADS(ctx->db);
 int i, j;
 char *read = New_Read_Buffer(ctx->db);
 for (i = 0; i < nContigs; i++)
 {
 int maxBeg = 0;
 int maxBegContigID = -1;
 int cLen = DB_READ_LEN(ctx->db, i);
 int minEnd = cLen;
 int minEndContigID = -1;
 for (j = 0; j < nContigs; j++)
 {
 int cutPos = ctx->LAStrimMatrix[i * nContigs + j];
 if (cutPos < 0 && abs(cutPos) > maxBeg)
 {
 maxBeg = abs(cutPos);
 maxBegContigID = j;
 }
 if (cutPos > 0 && cutPos < minEnd)
 {
 minEnd = cutPos;
 minEndContigID = j;
 }
 if (cutPos != 0)
 printf("FOUND CONTIG TRIM POSITION: CONTIG %d; TRIM: %d, TRIMLEN (%d) (OVL with: %d)\n", i, cutPos, (cutPos < 0) ? abs(cutPos) : cLen - cutPos, j);
 }
 float dustBegFract, dustEndFract, tanBegFract, tanEndFract;
 dustBegFract = dustEndFract = tanBegFract = tanEndFract = 0.0;
 if (maxBeg > 0 || minEnd != cLen)
 {
 if (maxBeg > 0)
 {
 dustBegFract = getMaskedBases(ctx, ctx->trackDust, i, 0, maxBeg) * 100.0 / maxBeg;
 tanBegFract = getMaskedBases(ctx, ctx->trackTan, i, 0, maxBeg) * 100.0 / maxBeg;
 }
 if (minEnd != cLen)
 {
 dustEndFract = getMaskedBases(ctx, ctx->trackDust, i, minEnd, cLen * 100.0 / cLen - minEnd);
 tanEndFract = getMaskedBases(ctx, ctx->trackTan, i, minEnd, cLen * 100.0 / cLen - minEnd);
 }

 printf(" --> final trim Interval: [%d, %d] -> trimmed [%d, %d] dustFract(in %%) [%.2f, %.2f] tanFract(in %%) [%.2f,%.2f]\n", maxBeg, minEnd, maxBeg, cLen - minEnd, dustBegFract, dustEndFract, tanBegFract, tanEndFract);

 }
 // int flags, qv;
 int map = 0;
 while (i < ctx->findx[map - 1])
 map -= 1;
 while (i >= ctx->findx[map])
 map += 1;

 Load_Read(ctx->db, i, read, 2);

 // write out trimmed contigs
 {
 fprintf(trimmedContigsAll, ">%s\n", ctx->flist[map]);
 for (j = maxBeg; j + ctx->lineWidth < minEnd; j += ctx->lineWidth)
 fprintf(trimmedContigsAll, "%.*s\n", ctx->lineWidth, read + j);
 if (j < minEnd)
 fprintf(trimmedContigsAll, "%.*s\n", minEnd - j, read + j);
 // write out purged sequence at begin of contig
 if (maxBeg > 0)
 {
 fprintf(purgedContigsAll, ">%s purged=%d,%d purgedLen=%d\n", ctx->flist[map], 0, maxBeg, maxBeg);
 for (j = 0; j + ctx->lineWidth < maxBeg; j += ctx->lineWidth)
 fprintf(purgedContigsAll, "%.*s\n", ctx->lineWidth, read + j);
 if (j < maxBeg)
 fprintf(purgedContigsAll, "%.*s\n", maxBeg - j, read + j);
 }
 // write out purged sequence at end of contig
 if (minEnd < cLen)
 {
 fprintf(purgedContigsAll, ">%s purged=%d,%d purgedLen=%d\n", ctx->flist[map], minEnd, cLen, cLen - minEnd);
 for (j = minEnd; j + ctx->lineWidth < cLen; j += ctx->lineWidth)
 fprintf(purgedContigsAll, "%.*s\n", ctx->lineWidth, read + j);
 if (j < cLen)
 fprintf(purgedContigsAll, "%.*s\n", cLen - j, read + j);
 }
 fprintf(statsContigsAll, "%d\t%s\t%d\t%d\t%d\ttrimBeg:LC=%.2f%%,TAN=%.2f%%;trimEnd=LC=%.2f%%,TAN=%.2f%%", i, ctx->flist[map], minEnd - maxBeg, maxBeg, cLen - minEnd, dustBegFract, tanBegFract, dustEndFract, tanEndFract);
 // contig support for trimBegin
 if (maxBeg != 0)
 {
 int bmap = 0;
 while (maxBegContigID < ctx->findx[bmap - 1])
 bmap -= 1;
 while (maxBegContigID >= ctx->findx[bmap])
 bmap += 1;
 fprintf(statsContigsAll, ";trimBegSupport:ID=%d,name=%s", maxBegContigID, ctx->flist[bmap]);
 }
 // contig support for trimEnd
 if (minEnd != cLen)
 {
 int bmap = 0;
 while (minEndContigID < ctx->findx[bmap - 1])
 bmap -= 1;
 while (minEndContigID >= ctx->findx[bmap])
 bmap += 1;
 fprintf(statsContigsAll, ";trimEndSupport:ID=%d,name=%s", minEndContigID, ctx->flist[bmap]);
 }
 fprintf(statsContigsAll, "\n");
 }
 // write out trimmed contigs but ignore tandem-induced overlaps
 {
 int tanMaxBeg = maxBeg;
 int tanMinEnd = minEnd;
 if (dustBegFract > ctx->maxLowCompTrimPerc || tanBegFract > ctx->maxLowCompTrimPerc)
 {
 tanMaxBeg = 0;
 }
 if (dustEndFract > ctx->maxLowCompTrimPerc || tanEndFract > ctx->maxLowCompTrimPerc)
 {
 tanMinEnd = cLen;
 }

 fprintf(trimmedContigsNoTandem, ">%s\n", ctx->flist[map]);
 for (j = tanMaxBeg; j + ctx->lineWidth < tanMinEnd; j += ctx->lineWidth)
 fprintf(trimmedContigsNoTandem, "%.*s\n", ctx->lineWidth, read + j);
 if (j < tanMinEnd)
 fprintf(trimmedContigsNoTandem, "%.*s\n", tanMinEnd - j, read + j);
 // write out purged sequence at begin of contig
 if (tanMaxBeg > 0)
 {
 fprintf(purgedContigsNoTandem, ">%s purged=%d,%d purgedLen=%d\n", ctx->flist[map], 0, tanMaxBeg, tanMaxBeg);
 for (j = 0; j + ctx->lineWidth < tanMaxBeg; j += ctx->lineWidth)
 fprintf(purgedContigsNoTandem, "%.*s\n", ctx->lineWidth, read + j);
 if (j < tanMaxBeg)
 fprintf(purgedContigsNoTandem, "%.*s\n", tanMaxBeg - j, read + j);
 }
 // write out purged sequence at end of contig
 if (tanMinEnd < cLen)
 {
 fprintf(purgedContigsNoTandem, ">%s purged=%d,%d purgedLen=%d\n", ctx->flist[map], tanMinEnd, cLen, cLen - tanMinEnd);
 for (j = tanMinEnd; j + ctx->lineWidth < cLen; j += ctx->lineWidth)
 fprintf(purgedContigsNoTandem, "%.*s\n", ctx->lineWidth, read + j);
 if (j < cLen)
 fprintf(purgedContigsNoTandem, "%.*s\n", cLen - j, read + j);
 }
 fprintf(statsContigsNoTandem, "%d\t%s\t%d\t%d\t%d\ttrimBeg:LC=%.2f%%,TAN=%.2f%%;trimEnd=LC=%.2f%%,TAN=%.2f%%", i, ctx->flist[map], tanMinEnd - tanMaxBeg, tanMaxBeg, cLen - tanMinEnd, dustBegFract, tanBegFract, dustEndFract, tanEndFract);
 // contig support for trimBegin
 if (tanMaxBeg != 0)
 {
 int bmap = 0;
 while (maxBegContigID < ctx->findx[bmap - 1])
 bmap -= 1;
 while (maxBegContigID >= ctx->findx[bmap])
 bmap += 1;
 fprintf(statsContigsNoTandem, ";trimBegSupport:ID=%d,name=%s", maxBegContigID, ctx->flist[bmap]);
 }
 // contig support for trimEnd
 if (tanMinEnd != cLen)
 {
 int bmap = 0;
 while (minEndContigID < ctx->findx[bmap - 1])
 bmap -= 1;
 while (minEndContigID >= ctx->findx[bmap])
 bmap += 1;
 fprintf(statsContigsNoTandem, ";trimEndSupport:ID=%d,name=%s", minEndContigID, ctx->flist[bmap]);
 }
 fprintf(statsContigsNoTandem, "\n");
 }
 }

 // write out bionano agp gaps
 if (ctx->BionanoAGPMatrix != NULL)
 {
 FILE *trimmedContigsBionano = NULL;
 FILE *purgedContigsBionano = NULL;
 FILE *statsContigsBionano = NULL;

 sprintf(fout, "%s.BionanoBasedTrimmedContigs.fasta", ctx->fileOutPattern);
 printf("create file: %s\n", fout);
 if ((trimmedContigsBionano = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }
 sprintf(fout, "%s.BionanoBasedPurgedContigs.fasta", ctx->fileOutPattern);
 if ((purgedContigsBionano = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }
 sprintf(fout, "%s.BionanoBasedContigs.stats", ctx->fileOutPattern);
 if ((statsContigsBionano = (FILE*) fopen(fout, "w")) == NULL)
 {
 fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
 exit(1);
 }

 fprintf(statsContigsBionano, "#ContigID\tContigName\tnewContigLength\ttrimBegin\ttrimEnd\tcomments\n");

 for (i = 0; i < nContigs; i++)
 {
 int maxBeg = 0;
 int maxBegContigID = -1;
 int cLen = DB_READ_LEN(ctx->db, i);
 int minEnd = cLen;
 int minEndContigID = -1;
 for (j = 0; j < nContigs; j++)
 {
 int cutPos = ctx->LAStrimMatrix[i * nContigs + j];
 int bionanoGap = ctx->BionanoAGPMatrix[i * nContigs + j];

 if (bionanoGap != 0 && cutPos != 0)
 {
 if (abs(bionanoGap) <= ctx->minBionanoGapLen)
 {
 printf("found matching bionano gap and contig ovl for contigs: %d vs %d, OVL: %d, GAP: %d\n", i, j, cutPos, bionanoGap);

 if (cutPos < 0 && abs(cutPos) > maxBeg)
 {
 maxBeg = abs(cutPos);
 maxBegContigID = j;
 }
 if (cutPos > 0 && cutPos < minEnd)
 {
 minEnd = cutPos;
 minEndContigID = j;
 }
 }
 else
 {
 printf("found matching bionano gap and contig ovl for contigs: %d vs %d, OVL: %d, GAP: %d - BUT GAP IS TO LARGE\n", i, j, cutPos, bionanoGap);
 }
 }
 else
 {
 if (bionanoGap != 0)
 {
 if (abs(bionanoGap) <= ctx->minBionanoGapLen)
 {
 printf("found bionano gap BUT NO contig ovl for contigs: %d vs %d, OVL: %d, GAP: %d\n", i, j, cutPos, bionanoGap);
 if (abs(bionanoGap) < 3000)
 printf("BIONANO GAP TO SMALL to find chains %d!\n", abs(bionanoGap));
 if (i < j)
 ctx->statsBionanoGapsMissed++;
 }
 }
 else if (cutPos != 0)
 {
 printf("found NO bionano gap but contig ovl for contigs: %d vs %d, OVL: %d, GAP: %d\n", i, j, cutPos, bionanoGap);
 }
 }
 }
 float dustBegFract, dustEndFract, tanBegFract, tanEndFract;
 dustBegFract = dustEndFract = tanBegFract = tanEndFract = 0.0;
 if (maxBeg > 0 || minEnd != cLen)	// contig i must be trimmed either at end or at begin, this corresponds with the Bionano AGP file
 {
 if (maxBeg > 0)
 {
 dustBegFract = getMaskedBases(ctx, ctx->trackDust, i, 0, maxBeg) * 100.0 / maxBeg;
 tanBegFract = getMaskedBases(ctx, ctx->trackTan, i, 0, maxBeg) * 100.0 / maxBeg;
 }
 if (minEnd != cLen)
 {
 dustEndFract = getMaskedBases(ctx, ctx->trackDust, i, minEnd, cLen * 100.0 / cLen - minEnd);
 tanEndFract = getMaskedBases(ctx, ctx->trackTan, i, minEnd, cLen * 100.0 / cLen - minEnd);
 }
 ctx->statsBionanoTrimmedContigs++;
 ctx->statsBionanoTrimmedBases += abs(maxBeg) + (cLen - minEnd);
 printf(" --> final trim Interval: [%d, %d] -> trimmed [%d, %d] dustFract(in %%) [%.2f, %.2f] tanFract(in %%) [%.2f,%.2f]\n", maxBeg, minEnd, maxBeg, cLen - minEnd, dustBegFract, dustEndFract, tanBegFract, tanEndFract);
 }

 int map = 0;
 while (i < ctx->findx[map - 1])
 map -= 1;
 while (i >= ctx->findx[map])
 map += 1;

 Load_Read(ctx->db, i, read, 2);

 // write out trimmed contigs
 {
 fprintf(trimmedContigsBionano, ">%s\n", ctx->flist[map]);
 for (j = maxBeg; j + ctx->lineWidth < minEnd; j += ctx->lineWidth)
 fprintf(trimmedContigsBionano, "%.*s\n", ctx->lineWidth, read + j);
 if (j < minEnd)
 fprintf(trimmedContigsBionano, "%.*s\n", minEnd - j, read + j);
 // write out purged sequence at begin of contig
 if (maxBeg > 0)
 {
 fprintf(purgedContigsBionano, ">%s purged=%d,%d purgedLen=%d\n", ctx->flist[map], 0, maxBeg, maxBeg);
 for (j = 0; j + ctx->lineWidth < maxBeg; j += ctx->lineWidth)
 fprintf(purgedContigsBionano, "%.*s\n", ctx->lineWidth, read + j);
 if (j < maxBeg)
 fprintf(purgedContigsBionano, "%.*s\n", maxBeg - j, read + j);
 }
 // write out purged sequence at end of contig
 if (minEnd < cLen)
 {
 fprintf(purgedContigsBionano, ">%s purged=%d,%d purgedLen=%d\n", ctx->flist[map], minEnd, cLen, cLen - minEnd);
 for (j = minEnd; j + ctx->lineWidth < cLen; j += ctx->lineWidth)
 fprintf(purgedContigsBionano, "%.*s\n", ctx->lineWidth, read + j);
 if (j < cLen)
 fprintf(purgedContigsBionano, "%.*s\n", cLen - j, read + j);
 }
 fprintf(statsContigsBionano, "%d\t%s\t%d\t%d\t%d\ttrimBeg:LC=%.2f%%,TAN=%.2f%%;trimEnd=LC=%.2f%%,TAN=%.2f%%", i, ctx->flist[map], minEnd - maxBeg, maxBeg, cLen - minEnd, dustBegFract, tanBegFract, dustEndFract, tanEndFract);
 // contig support for trimBegin
 if (maxBeg != 0)
 {
 int bmap = 0;
 while (maxBegContigID < ctx->findx[bmap - 1])
 bmap -= 1;
 while (maxBegContigID >= ctx->findx[bmap])
 bmap += 1;
 fprintf(statsContigsBionano, ";trimBegSupport:ID=%d,name=%s", maxBegContigID, ctx->flist[bmap]);
 }
 // contig support for trimEnd
 if (minEnd != cLen)
 {
 int bmap = 0;
 while (minEndContigID < ctx->findx[bmap - 1])
 bmap -= 1;
 while (minEndContigID >= ctx->findx[bmap])
 bmap += 1;
 fprintf(statsContigsBionano, ";trimEndSupport:ID=%d,name=%s", minEndContigID, ctx->flist[bmap]);
 }
 fprintf(statsContigsBionano, "\n");
 }
 }
 fclose(trimmedContigsBionano);
 fclose(purgedContigsBionano);
 fclose(statsContigsBionano);
 }

 fclose(trimmedContigsAll);
 fclose(purgedContigsAll);
 fclose(statsContigsAll);

 fclose(trimmedContigsNoTandem);
 fclose(purgedContigsNoTandem);
 fclose(statsContigsNoTandem);

 free(read - 1);
 free(fout);
 }*/

void trim_contigs(TrimContext *ctx)
{
	assert(ctx != NULL);
	int i, j, k, l;

	// split bionano gaps only
	if (ctx->purgeOpt == 0)
	{
		j = k = 0;

		while (j < ctx->numTrimEvidence)
		{
			while (k < ctx->numTrimEvidence - 1 && ctx->trimEvid[j].contigA == ctx->trimEvid[k + 1].contigA)
			{
				k++;
			}

			int n = k - j + 1;

			int aLen = DB_READ_LEN(ctx->db, ctx->trimEvid[j].contigA);
			int maxStart = 1;
			int minEnd = aLen;
			int tmp = 0;
			for (i = 0; i < n; i++)
			{
				TrimEvidence *te = ctx->trimEvid + j + i;
				for (l = 0; l < te->nBioNanoGaps; l++)
				{
					if (te->gaps[l].bionanoGapSize < ctx->minBionanoGapLen)
					{

						// todo: remember min and max cut positions
						//trimmedContigs++;

						printBionanpGap(ctx, te->contigA, te->contigB, te->gaps + l);
						if (te->gaps[l].aEnd == 1)
						{
							tmp = (1 + abs(te->gaps[l].bionanoGapSize) / 2 + ctx->trimOffset);
							if (tmp > maxStart)
							{
								maxStart = tmp;
							}
						}
						else if (te->gaps[l].aEnd == aLen)
						{
							tmp = aLen - (1 + abs(te->gaps[l].bionanoGapSize) / 2 + ctx->trimOffset);
							if (tmp < minEnd)
							{
								minEnd = tmp;
							}
						}

					}
				}
			}
			if (maxStart > 1 || minEnd < aLen)
			{
				printf("CUT POSITIONS (%d - %d, %d): %d, %d\n", ctx->trimEvid[j].contigA, 0, aLen, maxStart, minEnd);
				ctx->trimCoord[ctx->trimEvid[j].contigA].coord[0] =  maxStart;
				ctx->trimCoord[ctx->trimEvid[j].contigA].coord[1] =  minEnd;
			}
			k++;
			j = k;
		}
		printf("[INFO] num trim evidence: %d: \n", ctx->numTrimEvidence);
	}
	else if (ctx->purgeOpt == 1)
	{

		j = k = 0;

		while (j < ctx->numTrimEvidence)
		{
			while (k < ctx->numTrimEvidence - 1 && ctx->trimEvid[j].contigA == ctx->trimEvid[k + 1].contigA)
			{
				k++;
			}

			int n = k - j + 1;

			int aLen = DB_READ_LEN(ctx->db, ctx->trimEvid[j].contigA);
			int newStart = -1;
			int newEnd = -aLen;
			int tmp = 0;

			// ensure that buffer for cut positions is always big enough
			if(ctx->trimEvid[j].nBioNanoGaps > 1)
			{
				int tmp = (1+ctx->trimEvid[j].nBioNanoGaps);
				ctx->trimCoord[ctx->trimEvid[j].contigA].coord = (int*) realloc(ctx->trimCoord[ctx->trimEvid[j].contigA].coord, sizeof(int)*2*tmp);
				for (j=1; j<tmp; j++)
				{
					ctx->trimCoord[ctx->trimEvid[j].contigA].coord[2*j]=1;
					ctx->trimCoord[ctx->trimEvid[j].contigA].coord[2*j+1]=aLen;
				}
				ctx->trimCoord[ctx->trimEvid[j].contigA].maxCoordPairs = tmp;
			}

			for (i = 0; i < n; i++)
			{
				TrimEvidence *te = ctx->trimEvid + j + i;

				for (l = 0; l < te->nBioNanoGaps; l++)
				{

					// normal gap, i.e. bionano did not cut the contigs before scaffolding them
					if(((te->gaps[l].aEnd == 1) || (te->gaps[l].aEnd == aLen)) && ((te->gaps[l].bBeg == 1) || (te->gaps[l].aBeg == aLen)))
					{
						// check for negative gap
						if (te->gaps[l].bionanoGapSize < ctx->minBionanoGapLen)
						{
							if (te->gaps[l].aEnd == 1)
							{
								tmp = (1 + abs(te->gaps[l].bionanoGapSize) / 2 + ctx->trimOffset);
								// newStart
								int m;
								for (m=0; m < ctx->trimCoord[ctx->trimEvid[j].contigA].maxCoordPairs; m++)
								{
									if(intersect(
											abs(ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2]),
											abs(ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2+1]),
											te->gaps[l].aEnd, te->gaps[l].aBeg)
									)
									{
										 break;
									}
								}
								assert(m < ctx->trimCoord[ctx->trimEvid[j].contigA].maxCoordPairs);
								ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2] = tmp;
								ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2+1] = te->gaps[l].aBeg;
							}
							else if (te->gaps[l].aEnd == aLen)
							{
								tmp = aLen - (1 + abs(te->gaps[l].bionanoGapSize) / 2 + ctx->trimOffset);
								// newEnd
								int m;
								for (m=0; m < ctx->trimCoord[ctx->trimEvid[j].contigA].maxCoordPairs; m++)
								{
									if(intersect(
											abs(ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2]),
											abs(ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2+1]),
											te->gaps[l].aBeg, te->gaps[l].aEnd)
									)
									{
										 break;
									}
								}
								assert(m < ctx->trimCoord[ctx->trimEvid[j].contigA].maxCoordPairs);
								ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2] = te->gaps[l].aBeg;
								ctx->trimCoord[ctx->trimEvid[j].contigA].coord[m*2+1] = tmp;
							}
							else // should not occur
							{
								assert(0);
							}
						}
					}
					else // at least one contig at current gap was split before the scaffolding
					{


					}
				}
			}
			k++;
			j = k;
		}
		printf("[INFO] num trim evidence: %d: \n", ctx->numTrimEvidence);

	}
	else if (ctx->purgeOpt == 2)
	{

	}
	else if (ctx->purgeOpt == 3)
	{

	}
	else if (ctx->purgeOpt == 4)
	{

	}
	else if (ctx->purgeOpt == 5)
	{

	}
	else
	{
		/// UNKNOWN OPTION
	}

	// CUT CONTIGS!!!
	// open file handler
	FILE *trimmedContigs = NULL;
	FILE *removedContigParts = NULL;
	FILE *statsContigs = NULL;

	char *fout = malloc(strlen(ctx->fileOutPattern) + 50);
	assert(fout != NULL);

	sprintf(fout, "%s.trimmedContigs.fasta", ctx->fileOutPattern);
	printf("create file: %s\n", fout);
	if ((trimmedContigs = (FILE*) fopen(fout, "w")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
		exit(1);
	}

	sprintf(fout, "%s.removedContigParts.fasta", ctx->fileOutPattern);
	if ((removedContigParts = (FILE*) fopen(fout, "w")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
		exit(1);
	}
	sprintf(fout, "%s.trim.stats", ctx->fileOutPattern);
	if ((statsContigs = (FILE*) fopen(fout, "w")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
		exit(1);
	}
	free(fout);
	fprintf(statsContigs, "#ContigID\tContigName\tnewContigLength\ttrimBegin\ttrimEnd\tcomments\n");
	char *read = New_Read_Buffer(ctx->db);

	for (i = 0; i < DB_NREADS(ctx->db); i++)
	{
		// debug report trim position
		Load_Read(ctx->db, i, read, 2);
		TrimCoordinates *tc = ctx->trimCoord + i;
		int aLen = DB_READ_LEN(ctx->db, i);

		int amap = 0;
		while (i < ctx->findx[amap - 1])
			amap -= 1;
		while (i >= ctx->findx[amap])
			amap += 1;

		// write out removed part a the beginning !!!
		if (tc->coord[0] > 1)
		{
			ctx->statsRemovedContigPartBases += tc->coord[0];
			ctx->statsRemovedContigParts++;
			fprintf(removedContigParts, ">%s part=%d,%d\n", ctx->flist[amap], 1, tc->coord[0]);
			for (k = 0; k + ctx->lineWidth < tc->coord[0]; k += ctx->lineWidth)
				fprintf(removedContigParts, "%.*s\n", ctx->lineWidth, read + k);
			if (k < tc->coord[0])
				fprintf(removedContigParts, "%.*s\n", tc->coord[0] - k, read + k);
		}

		for (j = 0; j < tc->numCoordPairs; j++)
		{
			int index = j*2;
			ctx->statsTrimmedBases += tc->coord[index+1] - tc->coord[index];
			ctx->statsTrimmedContigs++;
			fprintf(trimmedContigs, ">%s part=%d,%d\n", ctx->flist[amap], tc->coord[index], tc->coord[index+1]);
			for (k = tc->coord[index]; k + ctx->lineWidth < tc->coord[index+1]; k += ctx->lineWidth)
				fprintf(trimmedContigs, "%.*s\n", ctx->lineWidth, read + k);
			if (k < tc->coord[index+1])
				fprintf(trimmedContigs, "%.*s\n", tc->coord[index+1] - k, read + k);

			if (j+1 < tc->numCoordPairs)
			{
				assert(tc->coord[index+2] > tc->coord[index+1]);
				ctx->statsRemovedContigPartBases += tc->coord[index+2]-tc->coord[index+1];
				ctx->statsRemovedContigParts++;
				fprintf(removedContigParts, ">%s part=%d,%d\n", ctx->flist[amap], tc->coord[index+1], tc->coord[index]);
				for (k = tc->coord[index+1]; k + ctx->lineWidth < tc->coord[index+2]; k += ctx->lineWidth)
					fprintf(removedContigParts, "%.*s\n", ctx->lineWidth, read + k);
				if (k < tc->coord[index+2])
					fprintf(removedContigParts, "%.*s\n", tc->coord[index+2] - k, read + k);
			}
		}

		// write out removed part at the end!!!
		if (tc->coord[tc->numCoordPairs * 2 - 1] < aLen)
		{
			ctx->statsRemovedContigPartBases += aLen - tc->coord[tc->numCoordPairs * 2 - 1];
			ctx->statsRemovedContigParts++;
			fprintf(removedContigParts, ">%s part=%d,%d\n", ctx->flist[amap], tc->coord[tc->numCoordPairs * 2 - 1], aLen);
			for (k = tc->coord[tc->numCoordPairs * 2 - 1]; k + ctx->lineWidth < aLen; k += ctx->lineWidth)
				fprintf(removedContigParts, "%.*s\n", ctx->lineWidth, read + k);
			if (k < aLen)
				fprintf(removedContigParts, "%.*s\n", aLen - k, read + k);

		}
	}
	fclose(trimmedContigs);
	fclose(removedContigParts);
	fclose(statsContigs);
	free(read -1);
}

void usage()
{
fprintf(stderr, "[-v] [-GTLOFwt <int>] [-ago <file>] [-dt <track>] <db> <contigs_out_prefix>\n");

fprintf(stderr, "options: -v        verbose\n");
fprintf(stderr, "         -d <trc>  low complexity track (e.g. dust)\n");
fprintf(stderr, "         -t <trc>  tandem repeat track  (e,f, tan)\n");
fprintf(stderr, "         -a <file> bionano agp file");
fprintf(stderr, "                   If a bionano-agp file is given, then only gaps up the minimum gaps size of (default: %d) and a valid overlap chain are trimmed\n", MIN_BIONANO_GAP_SIZE);
fprintf(stderr, "         -g <file> bionano gap file");
fprintf(stderr, "         -o <file> overlap file that was filtered with LAfilterChain");
fprintf(stderr, "         -G <int>  min Bionano gap size (default: %d)\n", MIN_BIONANO_GAP_SIZE);
fprintf(stderr, "         -T <int>  maximum trim length (default: -1)\n");
fprintf(stderr, "         -L <int>  maximum tandem repeat overlap fraction (in %%) (default: %d, valid range: [0,100])\n", MAX_TANDEMTRIM_PERC);
fprintf(stderr, "         -O <int>  trim offset in bases (default %d), i.e. in best case (if we have single overlap between 2 contigs) a gap of size 2xtrim_offset is created )\n", TRIM_OFFSET);
fprintf(stderr, "                   in case a valid alignment chain consisting of multiple alignments is present (representing heterozygous variations). The first last and the last alignment are used, (- trimOffset and + trimOffset, accordingly) \n");
fprintf(stderr, "                   (- trimOffset and + trimOffset, accordingly) creates a larger gap size, but heopefully removes the heterozygous difference.\n");
fprintf(stderr, "         -F <int>  number of fuzzy bases for chaining. For sanity check only, must correspond to -f of LAfilterChain. (default: %d)\n", FUZZY_BASES);
fprintf(stderr, "         -w <int>  specify number of characters per fasta line (default: %d)\n", FASTA_LINEWIDTH);
fprintf(stderr, "         -p <int>  trim Options, as follows:\n");
fprintf(stderr, "         					0: only trim negative bionano gaps (requires bionano AGP and GAP files)\n");
fprintf(stderr, "         					1: -t 0 AND split and trim contigs if bionano AGP file reports a contig break (requires bionano AGP and GAP files)\n");
fprintf(stderr, "         					2: trimming based on LAS chains only (requires a chain filtered overlap file)\n");
fprintf(stderr, "         					3: -t 2 BUT tandem induced contig overlap are not trimmed (requires chain filtered overlap file and tandem repeat track) - \n");
fprintf(stderr, "         					4: intersection of -t 0 and -t 2, i.e. LAS chains based trimming only if negative Bionano gaps support them + small negative Bionano gaps that are below the minimum daligner alignment length cutoff\n");
fprintf(stderr, "         					5: -t 4 AND LAS chain based trimming of small contigs that are below the bionano minimum length threshold (<150K)\n");
}

int main(int argc, char *argv[])
{
HITS_DB db;
TrimContext tctx;
PassContext *pctx = NULL;

FILE *fileOvlIn = NULL;

bzero(&tctx, sizeof(TrimContext));

tctx.db = &db;

// args

char *pcTrackDust = NULL;
char *pcTrackTan = NULL;

char *pathInBionanoAGP = NULL;
char *pathInBionanoGAP = NULL;

char *pcPathReadsIn = NULL;
char *pcPathOverlapsIn = NULL;

int c, tmp;

tctx.minBionanoGapLen = MIN_BIONANO_GAP_SIZE;
tctx.maxTrimLength = -1;
tctx.maxLowCompTrimPerc = MAX_TANDEMTRIM_PERC;
tctx.trimOffset = TRIM_OFFSET;
tctx.maxFuzzyBases = FUZZY_BASES;
tctx.lineWidth = FASTA_LINEWIDTH;

opterr = 0;
while ((c = getopt(argc, argv, "vd:t:a:g:o:G:T:L:O:F:w:p:")) != -1)
{
switch (c)
{
case 'v':
	tctx.verbose++;
	break;
case 'd':
	pcTrackDust = optarg;
	break;
case 't':
	pcTrackTan = optarg;
	break;
case 'a':
	pathInBionanoAGP = optarg;
	break;
case 'g':
	pathInBionanoGAP = optarg;
	break;
case 'o':
	pcPathOverlapsIn = optarg;
	break;
case 'G':
	tctx.minBionanoGapLen = atoi(optarg);
	break;
case 'T':
	tctx.maxTrimLength = atoi(optarg);
	break;
case 'L':
	tmp = atoi(optarg);
	if (tmp < 0 || tmp > 100)
	{
		fprintf(stderr, "[ERROR] Invalid range for tandem repeat fraction %d. Must be in [0,100]\n", tmp);
		exit(1);
	}
	tctx.maxLowCompTrimPerc = tmp;
	break;
case 'O':
	tctx.trimOffset = atoi(optarg);
	break;
case 'F':
	tctx.maxFuzzyBases = atoi(optarg);
	break;
case 'w':
	tctx.lineWidth = atoi(optarg);
	break;
case 'p':
	tctx.purgeOpt = atoi(optarg);
	break;

default:
	fprintf(stderr, "unknown option %c\n", optopt);
	usage();
	exit(1);
}
}

if (argc - optind != 2)
{
usage();
exit(1);
}

pcPathReadsIn = argv[optind++];
tctx.fileOutPattern = argv[optind++];

if (pcPathOverlapsIn)
{
if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL)
{
fprintf(stderr, "[ERROR] - could not open %s\n", pcPathOverlapsIn);
exit(1);
}
}

if (Open_DB(pcPathReadsIn, &db))
{
fprintf(stderr, "[ERROR] - could not open %s\n", pcPathReadsIn);
exit(1);
}

if (pcTrackDust)
{
tctx.trackDust = track_load(&db, pcTrackDust);
if (!tctx.trackDust)
{
fprintf(stderr, "[ERROR] - could not load track %s\n", pcTrackDust);
exit(1);
}
}

if (pcTrackTan)
{
tctx.trackTan = track_load(&db, pcTrackTan);
if (!tctx.trackTan)
{
fprintf(stderr, "[ERROR] - could not load track %s\n", pcTrackTan);
exit(1);
}
}

	// check if purge options are valid
if (tctx.purgeOpt == 0 || tctx.purgeOpt == 1)
{
if (pathInBionanoAGP == NULL || pathInBionanoGAP == NULL)
{
fprintf(stderr, "[ERROR] - trim option -p 0 and -p 1: requires a bionano.agp and bionano.gap file\n");
exit(1);
}
}
else if (tctx.purgeOpt == 2 || tctx.purgeOpt == 3)
{
if (fileOvlIn == NULL)
{
fprintf(stderr, "[ERROR] - trim option -p 2 and -p 3: requires a chain filtered LAS file\n");
exit(1);
}
if (tctx.purgeOpt == 3 && tctx.trackTan == NULL)
{
fprintf(stderr, "[ERROR] - trim option -p 3: requires a chain filtered LAS file and a tandem repeat file!\n");
exit(1);
}
}
else if (tctx.purgeOpt == 4 || tctx.purgeOpt == 5)
{
if (fileOvlIn == NULL || pathInBionanoAGP == NULL || pathInBionanoGAP == NULL)
{
fprintf(stderr, "[ERROR] - trim option -p 4 and -p 5: requires a chain filtered LAS file and a bionano.agp and a bionano.gap file\n");
exit(1);
}
}
else
{
fprintf(stderr, "[ERROR] - unknown trim option -p %d\n", tctx.purgeOpt);
usage();
exit(1);
}

if (fileOvlIn)
{
	// passes

pctx = pass_init(fileOvlIn, NULL);

pctx->split_b = 1;
pctx->load_trace = 1;
pctx->unpack_trace = 1;
pctx->data = &tctx;
pctx->write_overlaps = 0;
pctx->purge_discarded = 0;
}

trim_pre(pctx, &tctx);

getDBFastaHeader(&tctx, pcPathReadsIn);

if (pathInBionanoAGP)
{
parseBionanoAGPfile(&tctx, pathInBionanoAGP);
}
if (pathInBionanoGAP)
{
parseBionanoGAPfile(&tctx, pathInBionanoGAP);
}

if (fileOvlIn)
{
pass(pctx, trim_handler);
}

trim_contigs(&tctx);

trim_post(&tctx);

if (fileOvlIn)
{
pass_free(pctx);
fclose(fileOvlIn);
}

Close_DB(&db);

int i;
for (i = 0; i < tctx.nfiles; i++)
{
free(tctx.flist[i]);
free(tctx.hlist[i]);
}

for (i = 0; i < tctx.numTrimEvidence; i++)
{
if (tctx.trimEvid->nLASchains)
free(tctx.trimEvid[i].chains);
if (tctx.trimEvid->nBioNanoGaps)
free(tctx.trimEvid[i].gaps);
}
free(tctx.trimEvid);

free(tctx.flist);
free(tctx.hlist);
free(tctx.findx - 1);

for (i = 0; i < DB_NREADS(tctx.db); i++)
{
free(tctx.trimCoord[i].coord);
}
free(tctx.trimCoord);
return 0;
}
