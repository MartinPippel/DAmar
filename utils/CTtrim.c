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

#define MIN_BIONANO_GAP_SIZE 13
#define TRIM_OFFSET 100
#define FUZZY_BASES 1500
#define FASTA_LINEWIDTH 80

#define DEBUG_MASKING
#undef DEBUG_MASKING2

typedef struct {
	// stats counters
	int statsNumInvalidChains;
	int statsNumValidChains;
	int statsTrimmedContigs;
	int statsTrimmedBases;

	// db and I/O files
	HITS_DB *db;
	HITS_TRACK *trackDust;
	HITS_TRACK *trackTan;

	char *fileOutPattern;

	ovl_header_twidth twidth;

	// vector to store overlapping contigs: length: #contigs x #contigs
	int *LAStrimMatrix;

	// other options
	int verbose;
	int minBionanoGapLen;
	int maxTrimLength;
	int maxLowCompTrimPerc;
	int trimOffset;
	int maxFuzzyBases;
	int lineWidth;

	// fasta header
	int nfiles;
	char **flist;
	char **hlist;
	int *findx;
} TrimContext;

static void trim_pre(PassContext *pctx, TrimContext *tctx) {
	if (tctx->verbose) {
		printf( ANSI_COLOR_GREEN "PASS contig trimming\n" ANSI_COLOR_RESET);

		printf( ANSI_COLOR_RED "OPTIONS\n" ANSI_COLOR_RESET);
		printf( ANSI_COLOR_RED "  verbose %d\n" ANSI_COLOR_RESET,
				tctx->verbose);
		printf( ANSI_COLOR_RED "  minBionanoGapLen %d\n" ANSI_COLOR_RESET,
				tctx->minBionanoGapLen);
		printf( ANSI_COLOR_RED "  maxTrimLength %d\n" ANSI_COLOR_RESET,
				tctx->maxTrimLength);
		printf( ANSI_COLOR_RED "  maxLowCompTrimPerc %d\n" ANSI_COLOR_RESET,
				tctx->maxLowCompTrimPerc);
		printf( ANSI_COLOR_RED "  trimOffset %d\n" ANSI_COLOR_RESET,
				tctx->trimOffset);
		printf( ANSI_COLOR_RED "  maxFuzzyBases %d\n" ANSI_COLOR_RESET,
				tctx->maxFuzzyBases);

		if (tctx->trackDust)
			printf( ANSI_COLOR_RED "  dust Track %s\n" ANSI_COLOR_RESET,
					tctx->trackDust->name);
		if (tctx->trackTan)
			printf( ANSI_COLOR_RED "  tandem Track %s\n" ANSI_COLOR_RESET,
					tctx->trackTan->name);
	}

	tctx->twidth = pctx->twidth;
	tctx->LAStrimMatrix = (int*) malloc(
	DB_NREADS(tctx->db) * sizeof(int) * DB_NREADS(tctx->db));
	assert(tctx->LAStrimMatrix != NULL);
	bzero(tctx->LAStrimMatrix,
	DB_NREADS(tctx->db) * sizeof(int) * DB_NREADS(tctx->db));
}

static void trim_post(TrimContext *ctx) {
	if (ctx->verbose) {
		if (ctx->statsTrimmedContigs > 0) {
			printf("#trimmed contigs %d\n", ctx->statsTrimmedContigs);
		}

		if (ctx->statsTrimmedBases > 0) {
			printf("#trimmed bases: %d\n", ctx->statsTrimmedBases);
		}
		if (ctx->statsNumValidChains > 0) {
			printf("#valid chains %d\n", ctx->statsNumValidChains);
		}
		if (ctx->statsNumInvalidChains > 0) {
			printf("#invalid chains %d\n", ctx->statsNumInvalidChains);
		}

	}
	free(ctx->LAStrimMatrix);
}

static int getTrimPositions(TrimContext *ctx, Overlap *ovl, int pointA,
		int *cutA, int *cutB) {
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
	if (ovl->path.tlen) {
		ovl_trace *trace = ovl->path.trace;
		apos = abeg;
		bpos = ovl->path.bbpos;

		int j = 0;
		while (j < ovl->path.tlen) {
			//      printf("apos %6d, bpos %6d, oldDist %6d, newDist %6d\n", apos, bpos, dist, abs(pointA-((apos / twidth + 1) * twidth)));
			if (dist < abs(pointA - ((apos / twidth + 1) * twidth)))
				break;
			apos = (apos / twidth + 1) * twidth;
			bpos += trace[j + 1];
			//printf("apos %6d, bpos %6d\n", apos, bpos);
			dist = abs(pointA - apos);
			j += 2;
		}
	} else {
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

static int analyzeContigOverlaps(TrimContext *ctx, Overlap *ovl, int novl) {
	int i;

	int aLen = DB_READ_LEN(ctx->db, ovl->aread);
	int bLen = DB_READ_LEN(ctx->db, ovl->bread);

	// assumption: input overlaps must be chained with LAfilterChains !!!
	// one chain at the end and one chain at the beginning of a contig are possible!!!

	// sanity check
	Overlap *o1 = ovl;
	int cutA = -1;
	int cutB = -1;

	if (novl == 1) {
		// check if overlaps is valid

		if ((o1->path.abpos > ctx->maxFuzzyBases
				&& o1->path.aepos < aLen - ctx->maxFuzzyBases)
				|| (o1->path.bbpos > ctx->maxFuzzyBases
						&& o1->path.bepos < bLen - ctx->maxFuzzyBases)) {
			ctx->statsNumInvalidChains++;
			if (ctx->verbose) {
				printf(
						"[WARNGING] fuzzy base check failed! Ignore invalid chain [%d, %d] a[%d,%d] %c b[%d,%d]!\n",
						o1->aread, o1->bread, o1->path.abpos, o1->path.aepos,
						(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos,
						o1->path.bepos);
			}
			return 1;
		}

		ctx->statsNumValidChains++;
		int pointA = o1->path.abpos + (o1->path.aepos - o1->path.abpos) / 2;

		if (getTrimPositions(ctx, o1, pointA, &cutA, &cutB)) {
			printf(
					"Unable to get cutPosition for OVL [%d,%d] a[%d,%d] %c b[%d,%d] and pointA: %d\n",
					o1->aread, o1->bread, o1->path.abpos, o1->path.aepos,
					(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos,
					o1->path.bepos, pointA);
			return 1;
		}

		assert((cutA - ctx->trimOffset > 0) && (cutA + ctx->trimOffset < aLen));
		assert((cutB - ctx->trimOffset > 0) && (cutB + ctx->trimOffset < bLen));

		// set cut position of contig_A
		if (o1->path.abpos < aLen - o1->path.aepos) // trim off contig at begin
				{
			ctx->LAStrimMatrix[o1->aread * DB_NREADS(ctx->db) + o1->bread] =
					-(cutA + ctx->trimOffset);
		} else if (o1->path.abpos > aLen - o1->path.aepos) // trim off contig at end
				{
			ctx->LAStrimMatrix[o1->aread * DB_NREADS(ctx->db) + o1->bread] =
					cutA - ctx->trimOffset;
		} else // containment
		{
			printf(
					"Contained overlap: [%d,%d] a[%d,%d] %c b[%d,%d] and pointA: %d\n",
					o1->aread, o1->bread, o1->path.abpos, o1->path.aepos,
					(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos,
					o1->path.bepos, pointA);
			return 1;
		}

		// set cut position of contig_B
		if (o1->path.bbpos < bLen - o1->path.bepos) // trim off contig at begin
				{
			if (o1->flags & OVL_COMP) {
				ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db) + o1->aread] =
						bLen - (cutB + ctx->trimOffset);
			} else {
				ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db) + o1->aread] =
						-(cutB + ctx->trimOffset);
			}
		} else if (o1->path.bbpos > bLen - o1->path.bepos) // trim off contig at end
				{
			if (o1->flags & OVL_COMP) {
				ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db) + o1->aread] =
						-(bLen - (cutB - ctx->trimOffset));
			} else {
				ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db) + o1->aread] =
						cutB - ctx->trimOffset;
			}
		}
	} else {
		int validChain = 1;
		Overlap *o2;
		// first: sanity check for LAS chain
		for (i = 1; i < novl; i++) {
			o2 = ovl + i;
			if (abs(o1->path.aepos - o2->path.abpos) > ctx->maxFuzzyBases
					|| ((o1->flags & OVL_COMP) != (o2->flags & OVL_COMP))) {
				validChain = 0;
				break;
			}
			o1 = o2;
		}

		if (!validChain) {
			ctx->statsNumInvalidChains++;

			int mapA = 0;
			while (ovl->aread < ctx->findx[mapA - 1])
				mapA -= 1;
			while (ovl->aread >= ctx->findx[mapA])
				mapA += 1;

			int mapB = 0;
			while (ovl->bread < ctx->findx[mapB - 1])
				mapB -= 1;
			while (ovl->bread >= ctx->findx[mapB])
				mapB += 1;

			printf("INVALID chain: %d (%s) vs %d (%s)\n", ovl->aread,
					ctx->flist[mapA], ovl->bread, ctx->flist[mapB]);

			for (i = 0; i < novl; i++) {
				printf("   a[%d,%d] %c b[%d,%d]\n", ovl[i].path.abpos,
						ovl[i].path.aepos,
						(ovl[i].flags & OVL_COMP) ? 'c' : 'n',
						ovl[i].path.bbpos, ovl[i].path.bepos);
			}
			return 1;
		}
		ctx->statsNumValidChains++;

		o1 = ovl;
		o2 = ovl + (novl - 1);

		// set cut position of contig_A
		if (o1->path.abpos < aLen - o2->path.aepos) // trim off contig at begin
				{
			if (o2->path.abpos >= o1->path.aepos) {
				cutA = o2->path.abpos + ctx->trimOffset;
			} else {
				cutA = o1->path.aepos + ctx->trimOffset;
			}
			ctx->LAStrimMatrix[o1->aread * DB_NREADS(ctx->db) + o1->bread] =
					-(cutA);
		} else if (o1->path.abpos > aLen - o1->path.aepos) // trim off contig at end
				{
			if (o2->path.abpos >= o1->path.aepos) {
				cutA = o1->path.aepos - ctx->trimOffset;
			} else {
				cutA = o2->path.abpos - ctx->trimOffset;
			}
			ctx->LAStrimMatrix[o1->aread * DB_NREADS(ctx->db) + o1->bread] =
					cutA;
		} else // containment
		{
			printf("Contained overlap: [%d,%d] a[%d,%d] %c b[%d,%d]\n",
					o1->aread, o1->bread, o1->path.abpos, o1->path.aepos,
					(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos,
					o1->path.bepos);
			return 1;
		}
		// set cut position of contig_B

		if (o1->path.bbpos < bLen - o2->path.bepos) // trim off contig at begin
				{
			if (o2->path.bbpos >= o1->path.bepos) {
				cutB = o2->path.bbpos + ctx->trimOffset;
				if (o1->flags & OVL_COMP) {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = bLen - cutB;
				} else {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = -(cutB);
				}
			} else {
				cutB = o1->path.bepos + ctx->trimOffset;
				if (o1->flags & OVL_COMP) {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = bLen - cutB;
				} else {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = -(cutB);
				}
			}
		} else if (o1->path.bbpos > bLen - o2->path.bepos) // trim off contig at end
				{
			if (o2->path.bbpos >= o1->path.bepos) {
				cutB = o1->path.bepos - ctx->trimOffset;
				if (o1->flags & OVL_COMP) {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = -(bLen - cutB);
				} else {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = cutB;
				}
			} else {
				cutB = o2->path.bbpos - ctx->trimOffset;
				if (o1->flags & OVL_COMP) {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = -(bLen - cutB);
				} else {
					ctx->LAStrimMatrix[o1->bread * DB_NREADS(ctx->db)
							+ o1->aread] = cutB;
				}

			}
		} else // containment
		{
			printf("Contained overlap: [%d,%d] a[%d,%d] %c b[%d,%d]\n",
					o1->aread, o1->bread, o1->path.abpos, o1->path.aepos,
					(o1->flags & OVL_COMP) ? 'c' : 'n', o1->path.bbpos,
					o1->path.bepos);
			return 1;
		}
	}
	return 0;
}

static int trim_handler(void *_ctx, Overlap *ovl, int novl) {
	TrimContext *ctx = (TrimContext*) _ctx;

	// analyze overlaps and find contig trim position
	analyzeContigOverlaps(ctx, ovl, novl);

	return 1;
}

static int getMaskedBases(TrimContext *ctx, HITS_TRACK *t, int contigID,
		int beg, int end) {
#ifdef DEBUG_MASKING
	printf("call getMaskedBases on track %s, contigID: %d, in: [%d, %d]\n",
			t->name, contigID, beg, end);
#endif
	if (t == NULL) {
#ifdef DEBUG_MASKING        
		printf(" --> masked bases 0 (track is Null)\n");
#endif
		return 0;
	}

	track_anno *mask_anno = t->anno;
	track_data *mask_data = t->data;

	if (contigID < 0 || contigID >= DB_NREADS(ctx->db)) {
		fprintf(stderr,
				"[ERROR] - getMaskedBases contigID: %d out of bounds [0, %d]\n",
				contigID, DB_NREADS(ctx->db) - 1);
		fflush(stderr);
		exit(1);
	}

	track_anno rb, re;

	int maskBases = 0;
	int rBeg, rEnd;

	// repeat bases in a-read
	rb = mask_anno[contigID] / sizeof(track_data);
	re = mask_anno[contigID + 1] / sizeof(track_data);

	while (rb < re) {
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

static char *trimwhitespace(char *str)
{
	printf("trimwhitespace: %s\n",str);
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

static int getDBcontigID(TrimContext *ctx, char* contigName)
{
	printf("getDBcontigID(%s)\n", contigName);
	int i;
	for (i = 0; i < ctx->nfiles; i++)
	{
		if(strcmp(ctx->flist[i], contigName) == 0)
		{
			return i;
		}
	}
	return -1;
}

static void parseBionanoAGPfile(TrimContext *ctx, char *pathInBionanoAGP) {
	FILE *fileInBionanoGaps = NULL;

	if ((fileInBionanoGaps = fopen(pathInBionanoAGP, "r")) == NULL) {
		fprintf(stderr, "[ERROR] could not open %s\n", pathInBionanoAGP);
		exit(1);
	}

	char Obj_Name[MAX_NAME];
	int Obj_Start;
	int Obj_End;
	int PartNum;
	char Compnt_Type;
	char CompntId_GapLength[MAX_NAME];
	int CompntStart_GapType;
	int CompntEnd_Linkage;
	char Orientation_LinkageEvidence[MAX_NAME];

    char* line = NULL;
    size_t maxline = 0;

    int nline = 0;
    int len;
    char *pchrf, *pchrl;

    int r;
    int contigA;

    printf("parseBionanoAGPfile: %s\n", pathInBionanoAGP);
    while ((len = getline(&line, &maxline, fileInBionanoGaps)) > 0)
    {
        nline++;

        printf("line %d: %s\n", nline, line);
        char *tline = trimwhitespace(line);
        printf("line %d: %s\n", nline, tline);

        if (tline[0] == '#')
        	continue;

        printf("run scanf");
        r = scanf(tline, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", Obj_Name, &Obj_Start, &Obj_End, &PartNum, &Compnt_Type, CompntId_GapLength, &CompntStart_GapType, &CompntEnd_Linkage, Orientation_LinkageEvidence);
        printf("finish scanf: return %d\n",r);

        if( r != 9)
        {
        	fprintf(stderr, "[ERROR] invalid AGP file format %s. Expecting 9 columns, BUT parsed %d columns in line %d\n", pathInBionanoAGP, r, nline);
        	exit(1);
        }
        printf("line %d: %s\n", nline, tline);

        // try to match contig name with with DB contig ID
        contigA = getDBcontigID(ctx, CompntId_GapLength);
        if(contigA < 0)
        {
        	printf("Could not match agp contig name: %s in current db! Ignore AGP file.\n", CompntId_GapLength);
        	return;
        }
        else
        {
        	printf("found db contig id %d for agp contig name %s", contigA, CompntId_GapLength);
        }

	}

}

static void getDBFastaHeader(TrimContext *ctx, char *fullDBPath) {
	char *pwd, *root;
	FILE *dstub;
	int i;

	root = Root(fullDBPath, ".db");
	pwd = PathTo(fullDBPath);
	if (ctx->db->part > 0) {
		fprintf(stderr, "[ERROR] - CTtrim can not work on blocks!");
		exit(1);
	}
	dstub = Fopen(Catenate(pwd, "/", root, ".db"), "r");
	if (dstub == NULL) {
		fprintf(stderr, "[ERROR] - Cannot open database file: %s\n",
				Catenate(pwd, "/", root, ".db"));
		exit(1);
	}
	free(pwd);
	free(root);

	if (fscanf(dstub, DB_NFILE, &(ctx->nfiles)) != 1) {
		fclose(dstub);
		fprintf(stderr,
				"[ERROR] - Cannot read files line '%s' in database file: %s\n",
				DB_NFILE, Catenate(pwd, "/", root, ".db"));
		exit(1);
	}

	ctx->flist = (char**) Malloc(sizeof(char*) * ctx->nfiles,
			"Allocating file list");
	ctx->hlist = (char**) Malloc(sizeof(char*) * ctx->nfiles,
			"Allocating header list");
	ctx->findx = (int*) Malloc(sizeof(int*) * (ctx->nfiles + 1),
			"Allocating file index");
	if (ctx->flist == NULL || ctx->findx == NULL || ctx->hlist == NULL) {
		fclose(dstub);
		fprintf(stderr,
				"[ERROR] - Cannot allocate file name and file index buffers!");
		exit(1);
	}

	ctx->findx += 1;
	ctx->findx[-1] = 0;

	for (i = 0; i < ctx->nfiles; i++) {
		char headername[MAX_NAME], filename[MAX_NAME];

		if (fscanf(dstub, DB_FDATA, ctx->findx + i, filename, headername)
				!= 3) {
			fclose(dstub);
			fprintf(stderr,
					"[ERROR] - Cannot read %d-th fasta entry in database file %s\n",
					i + 1, Catenate(pwd, "/", root, ".db"));
			exit(1);
		}

		if ((ctx->flist[i] = Strdup(filename, "Adding to file list")) == NULL)
			exit(1);

		if ((ctx->flist[i] = Strdup(headername, "Adding to file list")) == NULL)
			exit(1);
	}

	fclose(dstub);
}

static void trim_contigs(TrimContext *ctx) {
	// open file handler
	FILE *trimmedContigsAll = NULL;
	FILE *purgedContigsAll = NULL;
	FILE *statsContigsAll = NULL;

	FILE *trimmedContigsNoTandem = NULL;
	FILE *purgedContigsNoTandem = NULL;
	FILE *statsContigsNoTandem = NULL;

	FILE *trimmedContigsBionano = NULL;
	FILE *purgedContigsBionano = NULL;
	FILE *statsContigsBionano = NULL;

	char *fout = malloc(strlen(ctx->fileOutPattern) + 30);
	assert(fout != NULL);

	sprintf(fout, "%s.trimmedContigs.fasta", ctx->fileOutPattern);
	printf("create file: %s\n", fout);
	if ((trimmedContigsAll = (FILE*) fopen(fout, "w")) == NULL) {
		fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
		exit(1);
	}
	sprintf(fout, "%s.purgedContigs.fasta", ctx->fileOutPattern);
	if ((purgedContigsAll = (FILE*) fopen(fout, "w")) == NULL) {
		fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
		exit(1);
	}
	sprintf(fout, "%s.trimmedContigs.stats", ctx->fileOutPattern);
	if ((statsContigsAll = (FILE*) fopen(fout, "w")) == NULL) {
		fprintf(stderr, "[ERROR] - could not open file %s\n", fout);
		exit(1);
	}

	fprintf(statsContigsAll,
			"#ContigID\tContigName\tnewContigLength\ttrimBegin\ttrimEnd\tcomments\n");

	// debug report trim positions
	int nContigs = DB_NREADS(ctx->db);
	int i, j;
	char *read = New_Read_Buffer(ctx->db);
	for (i = 0; i < nContigs; i++) {
		int maxBeg = 0;
		int maxBegContigID = -1;
		int cLen = DB_READ_LEN(ctx->db, i);
		int minEnd = cLen;
		int minEndContigID = -1;
		for (j = 0; j < nContigs; j++) {
			int cutPos = ctx->LAStrimMatrix[i * nContigs + j];
			if (cutPos < 0 && abs(cutPos) > maxBeg) {
				maxBeg = abs(cutPos);
				maxBegContigID = j;
			}
			if (cutPos > 0 && cutPos < minEnd) {
				minEnd = cutPos;
				minEndContigID = j;
			}
			if (cutPos != 0)
				printf(
						"FOUND CONTIG TRIM POSITION: CONTIG %d; TRIM: %d, TRIMLEN (%d) (OVL with: %d)\n",
						i, cutPos, (cutPos < 0) ? abs(cutPos) : cLen - cutPos,
						j);
		}
		float dustBegFract, dustEndFract, tanBegFract, tanEndFract;
		dustBegFract = dustEndFract = tanBegFract = tanEndFract = 0.0;
		if (maxBeg > 0 || minEnd != cLen) {
			if (maxBeg > 0) {
				dustBegFract = getMaskedBases(ctx, ctx->trackDust, i, 0, maxBeg)
						* 100.0 / maxBeg;
				tanBegFract = getMaskedBases(ctx, ctx->trackTan, i, 0, maxBeg)
						* 100.0 / maxBeg;
			}
			if (minEnd != cLen) {
				dustEndFract = getMaskedBases(ctx, ctx->trackDust, i, minEnd,
						cLen * 100.0 / cLen - minEnd);
				tanEndFract = getMaskedBases(ctx, ctx->trackTan, i, minEnd,
						cLen * 100.0 / cLen - minEnd);
			}

			printf(
					" --> final trim Interval: [%d, %d] -> trimmed [%d, %d] dustFract(in %%) [%.2f, %.2f] tanFract(in %%) [%.2f,%.2f]\n",
					maxBeg, minEnd, maxBeg, cLen - minEnd, dustBegFract,
					dustEndFract, tanBegFract, tanEndFract);

		}
		// int flags, qv;
		int map = 0;
		HITS_READ *r = ctx->db->reads + i;
		while (i < ctx->findx[map - 1])
			map -= 1;
		while (i >= ctx->findx[map])
			map += 1;

		Load_Read(ctx->db, i, read, 2);

		// write out trimmed contigs
		fprintf(trimmedContigsAll, ">%s\n", ctx->flist[map]);
		for (j = maxBeg; j + ctx->lineWidth < minEnd; j += ctx->lineWidth)
			fprintf(trimmedContigsAll, "%.*s\n", ctx->lineWidth, read + j);
		if (j < minEnd)
			fprintf(trimmedContigsAll, "%.*s\n", minEnd - j, read + j);
		// write out purged sequence at begin of contig
		if (maxBeg > 0) {
			fprintf(purgedContigsAll, ">%s purged=%d,%d purgedLen=%d\n",
					ctx->flist[map], 0, maxBeg, maxBeg);
			for (j = 0; j + ctx->lineWidth < maxBeg; j += ctx->lineWidth)
				fprintf(purgedContigsAll, "%.*s\n", ctx->lineWidth, read + j);
			if (j < maxBeg)
				fprintf(purgedContigsAll, "%.*s\n", maxBeg - j, read + j);
		}
		// write out purged sequence at end of contig
		if (minEnd < cLen) {
			fprintf(purgedContigsAll, ">%s purged=%d,%d purgedLen=%d\n",
					ctx->flist[map], minEnd, cLen, cLen - minEnd);
			for (j = minEnd; j + ctx->lineWidth < cLen; j += ctx->lineWidth)
				fprintf(purgedContigsAll, "%.*s\n", ctx->lineWidth, read + j);
			if (j < cLen)
				fprintf(purgedContigsAll, "%.*s\n", cLen - j, read + j);
		}
		fprintf(statsContigsAll,
				"%d\t%s\t%d\t%d\t%d\ttrimBeg:LC=%.2f%%,TAN=%.2f%%;trimEnd=LC=%.2f%%,TAN=%.2f%%",
				i, ctx->flist[map], minEnd - maxBeg, maxBeg, cLen - minEnd,
				dustBegFract, tanBegFract, dustEndFract, tanEndFract);
		// contig support for trimBegin
		if (maxBeg != 0) {
			int bmap = 0;
			while (maxBegContigID < ctx->findx[bmap - 1])
				bmap -= 1;
			while (maxBegContigID >= ctx->findx[bmap])
				bmap += 1;
			fprintf(statsContigsAll, ";trimBegSupport:ID=%d,name=%s",
					maxBegContigID, ctx->flist[bmap]);
		}
		// contig support for trimEnd
		if (minEnd != cLen) {
			int bmap = 0;
			while (minEndContigID < ctx->findx[bmap - 1])
				bmap -= 1;
			while (minEndContigID >= ctx->findx[bmap])
				bmap += 1;
			fprintf(statsContigsAll, ";trimEndSupport:ID=%d,name=%s",
					minEndContigID, ctx->flist[bmap]);
		}
		fprintf(statsContigsAll, "\n");
	}

	fclose(trimmedContigsAll);
	fclose(purgedContigsAll);
	fclose(statsContigsAll);
}

static void usage() {
	fprintf(stderr,
			"[-v] [-GTLOFw <int>] [-b <file>] [-dt <track>] <db> <overlaps_in> <contigs_out_prefix>\n");

	fprintf(stderr, "options: -v        verbose\n");
	fprintf(stderr, "         -d <trc>  low complexity track (e.g. dust)\n");
	fprintf(stderr, "         -t <trc>  tandem repeat track  (e,f, tan)\n");
	fprintf(stderr,
			"         -b <file> gap csv file i.e. based on bionano agp. Must be in following format: +|-<ContigID1> <gapLen> +|-<contigID2>\n");
	fprintf(stderr,
			"                   where: ContigID1, and ContigID2 must refer to the DAmar seqeuence database (i.e. 0-based contig IDs)\n");
	fprintf(stderr,
			"                   the mandatory +|- prefix of the ContigID describes the orientation of the contig\n");
	fprintf(stderr,
			"                   If a bionano-gap file is given, then only gaps up the minimum gaps size of (default: %d) are trimmed\n",
			MIN_BIONANO_GAP_SIZE);
	fprintf(stderr,
			"                   --> idea behind this: If Bionano inserts a 13bp gap, then it's most probable that the adjacent contigs overlap with each other\n");
	fprintf(stderr, "         -G <int>  min Bionano gap size (default: %d)\n",
	MIN_BIONANO_GAP_SIZE);
	fprintf(stderr, "         -T <int>  maximum trim length (default: -1)\n");
	fprintf(stderr,
			"         -L <int>  do not trim contigs if trim length contains more then -S bases (in %%) of tandem repeats (default: -1, valid range: [0,100])\n");
	fprintf(stderr,
			"         -O <int>  trim offset in bases (default %d), i.e. in best case (if we have single overlap between 2 contigs) a gap of size 2xtrim_offset is created )\n",
			TRIM_OFFSET);
	fprintf(stderr,
			"                   in case a valid alignment chain consisting of multiple alignments is present (representing heterozygous variations). The first last and the last alignment are used, (- trimOffset and + trimOffset, accordingly) \n");
	fprintf(stderr,
			"                   (- trimOffset and + trimOffset, accordingly) creates a larger gap size, but heopefully removes the heterozygous difference.\n");
	fprintf(stderr, "         -F <int>  number of fuzzy bases (default: %d)\n",
	FUZZY_BASES);
	fprintf(stderr,
			"         -w <int>  specify number of characters per fasta line (default: %d)\n",
			FASTA_LINEWIDTH);
}

int main(int argc, char *argv[]) {
	HITS_DB db;
	TrimContext tctx;
	PassContext *pctx;

	FILE *fileOvlIn = NULL;

	bzero(&tctx, sizeof(TrimContext));

	tctx.db = &db;

// args

	char *pcTrackDust = NULL;
	char *pcTrackTan = NULL;

	char *pathInBionanoAGP = NULL;

	char *pcPathReadsIn = NULL;
	char *pcPathOverlapsIn = NULL;

	int c, tmp;

	tctx.minBionanoGapLen = MIN_BIONANO_GAP_SIZE;
	tctx.maxTrimLength = -1;
	tctx.maxLowCompTrimPerc = -1;
	tctx.trimOffset = TRIM_OFFSET;
	tctx.maxFuzzyBases = FUZZY_BASES;
	tctx.lineWidth = FASTA_LINEWIDTH;

	opterr = 0;
	while ((c = getopt(argc, argv, "vd:t:b:G:T:L:O:F:w:")) != -1) {
		switch (c) {
		case 'v':
			tctx.verbose = 1;
			break;
		case 'd':
			pcTrackDust = optarg;
			break;
		case 't':
			pcTrackTan = optarg;
			break;
		case 'b':
			pathInBionanoAGP = optarg;
			break;
		case 'G':
			tctx.minBionanoGapLen = atoi(optarg);
			break;
		case 'T':
			tctx.maxTrimLength = atoi(optarg);
			break;
		case 'L':
			tmp = atoi(optarg);
			if (tmp < 0 || tmp > 100) {
				fprintf(stderr,
						"[ERROR] Invalid range for tandem repeat fraction %d. Must be in [0,100]\n",
						tmp);
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

		default:
			fprintf(stderr, "unknown option %c\n", optopt);
			usage();
			exit(1);
		}
	}

	if (argc - optind != 3) {
		usage();
		exit(1);
	}

	pcPathReadsIn = argv[optind++];
	pcPathOverlapsIn = argv[optind++];
	tctx.fileOutPattern = argv[optind++];

	if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL) {
		fprintf(stderr, "[ERROR] - could not open %s\n", pcPathOverlapsIn);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db)) {
		fprintf(stderr, "[ERROR] - could not open %s\n", pcPathReadsIn);
		exit(1);
	}

	if (pcTrackDust) {
		tctx.trackDust = track_load(&db, pcTrackDust);
		if (!tctx.trackDust) {
			fprintf(stderr, "[ERROR] - could not load track %s\n", pcTrackDust);
			exit(1);
		}
	}

	if (pcTrackTan) {
		tctx.trackTan = track_load(&db, pcTrackTan);
		if (!tctx.trackTan) {
			fprintf(stderr, "[ERROR] - could not load track %s\n", pcTrackTan);
			exit(1);
		}
	}

	getDBFastaHeader(&tctx, pcPathReadsIn);

	if (pathInBionanoAGP) {
		parseBionanoAGPfile(&tctx, pathInBionanoAGP);
	}

// passes

	pctx = pass_init(fileOvlIn, NULL);

	pctx->split_b = 1;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &tctx;
	pctx->write_overlaps = 0;
	pctx->purge_discarded = 0;

	trim_pre(pctx, &tctx);

	pass(pctx, trim_handler);

	// todo
	// 1. do the trimming + create the stats + different outpu files + do not trim if tandem repeat dfraction is above a certain threshold
	// 2. read in bionano agp file + do the trimming based on short gaps
	// 3. cleanup + testing

	trim_contigs(&tctx);

	trim_post(&tctx);

	pass_free(pctx);

	//todo cleanup

	Close_DB(&db);
	fclose(fileOvlIn);

	return 0;
}
