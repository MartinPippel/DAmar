#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <sys/stat.h>

#include "db/DB.h"
#include "lib/pass.h"
#include "dalign/align.h"
#include "dalign/filter.h"

typedef struct
{
	int alnLen;
	float eRate;
	int unalignedBases;

	int trimPos;	// if < 0 trim at contig begin, else trim at contig end
} LASchain;

typedef struct
{
	int flag;
	// info from AGP file
	int aBeg;	// most times this is the full contig length,
	int aEnd; // only in case Bionano split contigs those are relevant

	int bBeg;
	int bEnd;

	int agpGapSize; 		// this value must be positive: [13,N]
	int bionanoGapSize; // this value can also be negative in case of overlapping contigs
} BionanoGap;

typedef struct
{
	int contigA;
	int contigB;

	// evidence from LAS files
	int nLASchains;
	int maxLASchains;
	LASchain *chains;

	// evidence from Bionano AGP and GAP
	int nBioNanoGaps;
	int maxBionanoGaps;
	BionanoGap *gaps;
} TrimEvidence;

int TrimEvidence_cmp(const void *x, const void *y)
{
	TrimEvidence *te1 = (TrimEvidence*) x;
	TrimEvidence *te2 = (TrimEvidence*) y;

	if (te1->contigA != te2->contigA)
		return (te1->contigA - te2->contigA);

	return (te1->contigB - te2->contigB);
}

typedef struct
{
	// stats counters
	int statsNumInvalidChains;
	int statsNumValidChains;

	int statsTrimmedContigs;
	int statsTrimmedBases;
	int statsBionanoTrimmedContigs;
	int statsBionanoGapsMissed;
	int statsBionanoTrimmedBases;

	int statsBionanoGapsLtMinThresh;
	int statsBionanoGapsLtMinThreshContigBreak;
	int statsBionanoGapsAll;

	int numTrimEvidence;
	int maxTrimEvidence;
	TrimEvidence *trimEvid;

	// db and I/O files
	HITS_DB *db;
	HITS_TRACK *trackDust;
	HITS_TRACK *trackTan;

	char *fileOutPattern;

	ovl_header_twidth twidth;

	// vector to store overlapping contigs: length: #contigs x #contigs
	//int *LAStrimMatrix;
	// bionano AGP vector!

	//int *BionanoAGPMatrix;
	// bionano Gap length information
	//int *BionanoGapMatrix;

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

/*
 * increase bionano gap buffer if necessary
 */
void ensureBionanoGapBuffer(TrimEvidence *t, int numNewElements);

/*
 * increase LAS chain buffer if necessary
 */
void ensureLASchainBuffer(TrimEvidence *t, int numNewElements);

/*
 * add Bionano gap feature as new TrimEvidence
 * bionano feature: Contig_A - GAP - Contig_B (plus oriention, i.e. coordinates of the contigs)
 */
void addBionanoAGPInfoToTrimEvidence(TrimContext *ctx, int contigA, int fromA, int toA, int contigB, int fromB, int toB, int gapLen);

/*
 * add Bionano gap feature to an existing TrimEvidence.
 * assumption: bionano AGP must be parsed before the Bionano GAP file, because Bionano gap features must be present before adding the actual gap information,
 * this is mandatory as the contig orientation is unknown in the gap files and will become ambiguous if Bionano splitted the input contigs
 */
void addBionanoGAPInfoToTrimEvidence(TrimContext *ctx, int contigA, int aPartBeg, int aPartEnd, int contigB, int bPartBeg, int bPartEnd, int AdjustedGapLength);

/*
 * add/create Contig-Vs-Contig overlap to an existing/a new TrimEvidence in an aymmetric way
 */
void addLASchainInfoToTrimEvidence(TrimContext *ctx, int aread, int bread, int alnLen, int unAlnLen, float erate, int cutPosInA);

/*
 * find if a TrimEvidence entity between contigA and contigB is present.
 * if yes: return  TrimEvidence pointer
 * if no: return NULL
 */
TrimEvidence* find_TrimEvidence(TrimContext *ctx, const int contigA, const int contigB);

/*
 * Insert a new trimEvidence entity for contigA and contigB
 * returns a newly created TrimEvidanec pointer
 */
TrimEvidence* insert_TrimEvidence(TrimContext *ctx, const int contigA, const int contigB);

/*
 * assign trim positions for contigA and contigB for a given overlap and a given trimPoint
 * assumption: contigA and contigB have a single proper overlap (i.e. not a chain with multiple overlaps)
 * if trace points are present, then find closest trace point to pointA, and assign cutA to this tracepoint and cutB to the corresonding coordinate in B
 *   otherwise: assign pointA to cutA, and set cutB to ovl->path.bbeg + (pointA-ovl->path.abeg)
 *
 *  return 0: if everything was ok,
 *  return 1: if something went wrong
 */
int getTrimPositionsFromLAS(TrimContext *ctx, Overlap *ovl, int pointA, int *cutA, int *cutB);

/*
 * Analyze alignment chains, and add potential trim positions to the TrimEvidence
 */
int analyzeContigOverlaps(TrimContext *ctx, Overlap *ovl, int novl);

/*
 * get masked bases for a given mask track, contig and coordinates
 */
int getMaskedBases(TrimContext *ctx, HITS_TRACK *t, int contigID, int beg, int end);

/*
 * trim of leading and tailing white spaces from a given string
 * used when parsing Bionano AGP and GAP files
 */
char* trimwhitespace(char *str);

/*
 * try to find a given contigContigName with a contig ID in the current database
 */
int getDBcontigID(TrimContext *ctx, char *contigName, int *from, int *to);

/*
 * parse Bionano Gap file, which contains real gap estimates (even negative sizes, when contigs overlap)
 * unfortunately: the contig orientations are missing, therefore we need the AGP file as well
 */
void parseBionanoGAPfile(TrimContext *ctx, char *pathInBionanoGAP);

/*
 * parse Bionano AGP file
 */
void parseBionanoAGPfile(TrimContext *ctx, char *pathInBionanoAGP);

/*
 * auxiliary function: get sequence names from the database file
 * its used to match the DB-contig names with the bionano-contig names
 */
void getDBFastaHeader(TrimContext *ctx, char *fullDBPath);

/*
 * This function trims the contigs based on the available TrimEvidence
 */
void trim_contigs(TrimContext *ctx);

/*
 * print out how to use CTtrim tool
 */
void usage();

void printBionanpGap(TrimContext *ctx, int contigA, int contigB, BionanoGap *g);
