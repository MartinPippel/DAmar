/*
 * analyzeContigs.h
 *
 *  Created on: 3 Aug 2017
 *      Author: pippel
 */

#ifndef UTILS_LAANALYZE_H_
#define UTILS_LAANALYZE_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "db/DB.h"
#include "dalign/align.h"

// if abs(beg) > abs(end), i.e. overlap is in complement orientation
// if beg < 0 || end < 0, i.e. abs(positions) are not exact, only derived from nearest trace point
typedef struct
{
	int patchedID;   // fixed read ID that overlaps with a fixed read X that was incorporated in a contig
	int beg;  // corresponding beg, end positions
	int end;  // within read ID

	int cBeg; // corresponding begin, end positions
	int cEnd; // within fixed read X --> is then converted into corresponding contig position

} OvlRead;


#define CONTIG_READ_VALID     (1 << 0)  // 	1 initially all reads are valid, because they are part of a contig
#define CONTIG_READ_CHECKED	 (1 << 1)  // 	2
#define CONTIG_READ_DISCARD   (1 << 2) // 	4 contig read should be completely cut out from contig
#define CONTIG_READ_NEWEND    (1 << 3)  // 	8c	ontig ends at current read at position contigPosEnd
#define CONTIG_READ_NEWSTART  (1 << 4)  // 	16	contig starts at current read at position: contigPosBeg
#define CONTIG_READ_BREAK_UNKNOWN       	(1 << 5)		//	32
#define CONTIG_READ_BREAK_MULTIREAD 		(1 << 6)		//	64
#define CONTIG_READ_BREAK_STARTEND 		(1 << 7)		//	128
#define CONTIG_READ_BREAK_FALSEJOINLEFT 	(1 << 8)		//	256
#define CONTIG_READ_BREAK_FALSEJOINRIGHT	(1 << 9)		//	512
#define CONTIG_READ_BREAK_HETEROZYGLEFT	(1 << 10)	//	1024
#define CONTIG_READ_BREAK_HETEROZYGRIGHT	(1 << 11)	//	2048
#define CONTIG_READ_BREAK_DEADENDLEFT	(1 << 12)	//	4096
#define CONTIG_READ_BREAK_DEADENDRIGHT	(1 << 13)	//	8192

typedef struct
{

	int type; // see above types, only contigPosBeg and contigPosEnd modified if necessary, (original position are still available ovlReads[0])

	// information according to corrected reads
	int correctedID;
	int correctedContigPosBeg;		// absolute coordinates in
	int correctedContigPosEnd;		// corrected contig
	int correctedReadPosBeg;			// absolute coordinates in
	int correctedReadPosEnd;			// corrected contig
	// information according to patched reads and patched LAS
	int patchedID;
	int patchedContigPosBeg;			// absolute coordinates in
	int patchedContigPosEnd;		  // uncorrected (patched only) contig

	int nOutReads;
	int nInReads;

	int nComReadsWithNextContigRead; // out

	float avgCov;
	int lowCov;    // yes or no, indicates a problem e.g. chimeric read

	OvlRead *ovlReads;
	int numOvlReads;
	int maxOvlReads;

} ContigRead;

//// overlap group flags
#define OVLGRP_COMP          (1 << 0)        // if overlaps are in complement direction
#define OVLGRP_DISCARD       (1 << 1)        // overlap group tagged for removal (if a contig is tagged invalid, or the ovlgrp is covers only a few bases, etc)
// all defines are in aread - bread direction
#define OVLGRP_AREAD_IS_CONTAINED  (1 << 2)  // a-read is contained in b read
#define OVLGRP_AREAD_HAS_CONTAINED (1 << 3)  // b-read is contained in a read
#define OVLGRP_TOO_SHORT     (1 << 4)        // overlap group is shorter than MIN_CONTIG_LENGTH


typedef struct
{
	Overlap **ovls;
	int novl;
	int maxOvl;
} Chain;


typedef struct
{
	int aread;
	int bread;

	int coveredBasesInA;
	int coveredBasesInB;
	int gapBasesInA;
	int gapBasesInB;
	int repeatBasesInA;
	int repeatBasesInB;

	int first_abpos, first_aepos;
	int first_bbpos, first_bepos;  // corresponding to bread

	int last_abpos, last_aepos;
	int last_bbpos, last_bepos;  // corresponding to bread

	int flag;
} OverlapGroup;

//// contig flags
#define CONTIG_UNIQUE        (1 << 0)        //	1	contigs has no valid overlap group
#define CONTIG_DISCARD       (1 << 1)        //	2	 contig tagged for removal
#define CONTIG_IS_CONTAINED  (1 << 2)        // 	4	contig is contained in another contig
#define CONTIG_HAS_CONTAINED (1 << 3)        // 	8	contig includes contained
#define CONTIG_UNCLASSIFIED  (1 << 4)        // 	16
#define CONTIG_NORELATION    (1 << 5)        // 	32
#define CONTIG_AL50PCTCOVERED (1 << 6)    	//  64

typedef struct _ContigFlag2Label ContigFlag2Label;
struct _ContigFlag2Label
{
	int mask;

	char* label;
	char indicator;
};

typedef struct
{
	int flag;

	int correspID;     /// contained in Contig id, || contains read id || has left exit with read id etc. ....
	int numCoveredIntervals;
	int *coveredIntervals; // begPos, endPos, avgCov --> positions according to the contig  where this ContigGraphClassification belongs to
	int nJointReads;
} ContigGraphClassification;

///// contig classification flags
#define CONTIG_CLASS_UNKNOWN			(1 << 0)			//	1
#define CONTIG_CLASS_HAPLOID			(1 << 1)     	// 	2	haploid contig
#define CONTIG_CLASS_ALT   			(1 << 2)     	// 	4	alternative contig (heterozygous difference)
#define CONTIG_CLASS_ALT_SPUR		(1 << 3)    		// 	8	alternative contig is a Spur
#define CONTIG_CLASS_ALT_BUBBLE		(1 << 4)    		// 	16	alternative contig is a Bubble
#define CONTIG_CLASS_ALT_HUGEDIFF	(1 << 5)  		// 	32	alternative contig with hiuge difference
#define CONTIG_CLASS_REPEAT			(1 << 6)     	// 	64	repeat
#define CONTIG_CLASS_WEIRD 			(1 << 7)     	// 	128	weird
#define CONTIG_VISITED 				(1 << 8)     	// 	256	visited flag

typedef struct
{
	int duplicRead;
	int contigPos1;
	int contigPos2;
	int type;

	int leftProperRead;
	int leftProperReadPos1;
	int leftProperReadPos2;

	int rightProperRead;
	int rightProperReadPos1;
	int rightProperReadPos2;
} SplitEvent;

typedef struct
{
	int corContigIdx;
	int *abpos;
	int *aepos;
	int numPos;
} ContigChain;

typedef struct
{
	int flag;              // (above) CONTIG flags
	int len;               // contig length
	int idx;               // initial index from Contig-DB

	ContigRead *cReads;
	int numcReads;
	float avgCov;

	// valid overlap chains with other contigs, i.e. all repeat induced overlaps are filtered out !!!!
	Chain *contigChains;
	int curContigChains;
	int maxContigChains;

	ContigChain *cChains;
	int numCChains;
	int maxCChains;

	ContigGraphClassification *gClassific;
	int gClassificFlag;
	int numGClassific;
	int maxGClassific;

	int tmp;
	int class;
} Contig;

typedef struct
{
	// corrected contigs
	char *corContigDBName;
	char *corContigLASName;

	HITS_DB* corContigDB;
	HITS_TRACK *corContigPatchReadsAndPos_track;
	HITS_TRACK *corContigCorrectReadsAndPos_track;
	HITS_TRACK *corContigRawReads_track;
	HITS_TRACK *corContigRepeats_track;

	// patched reads DB and LAS
	char *patchedReadDBName;
	char *patchedReadLASName;

	HITS_DB* patchedReadDB;
	HITS_TRACK *patchedReadSource_track;
	HITS_TRACK *patchedReadTrim_track;

	// corrected reads DB
	char *corReadDBName;
	HITS_DB* corReadDB;

	// some filter for patched reads LAS
	int min_rlen;
	int min_olen;
	int exp_cov;

	// save output files in sub dir
	char *outDir;
	// output stuff
	char *readSeq;

	int numContigs;
	int **vreadMask; // len: number of reads of fixedDB times number of times a read can occur in different contigs
	int maxVReadMask; // number of times the same read can occur in different contigs (should not be greater than 2-3??)

	Contig* contigs;
	Contig** contigs_sorted;

	int twidth;             // trace point spacing
	int VERBOSE;
	int SPUR_LEN;
	int TIP_LEN;

	// TODO change this, use file indexes and not names!!!
	int *contigFileNameOffsets;
	char **ContigFileNames;
	int numContigFileNames;

	// detect and store temporarily contig chains
	Chain *ovlChains;
	int curChains;
	int maxChains;
} AnalyzeContext;

void initAnalyzeContext(AnalyzeContext *actx);
void finalContigValidation(AnalyzeContext *actx);

//void printSplitEvent(Contig *contig, SplitEvent *split);
//void printFinalAltContigClassification(FILE *out, FinalContigAltClass facc);
//void printFinalContigClassification(FILE *out, FinalContigClass fcc);
int getFullReadPositionInContig(AnalyzeContext *actx, Contig* c, int readID, int *cbeg, int *cend);
ContigRead* getContigRead(Contig* contig, int readID);

char* getFastaFileNameFromDB(AnalyzeContext *actx, int contigID);
int getNumberOfSequencesFromFastaFileName(AnalyzeContext *actx, int contigID);
int getFirstDBIdxFromFastaFileName(AnalyzeContext *actx, int contigID);
int createOutDir(char *out);

//////// analyze contigs based on fixed-read overlaps and intersection of read ids
int processFixedReadOverlap_handler1(void* _ctx, Overlap* ovls, int novl);
int processFixedReadOverlap_handler2(void* _ctx, Overlap* ovls, int novl);
void getCorrespondingPositionInBRead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2);
void getCorrespondingPositionInARead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2);
/// do coverage analysis for contig-read-overlaps, set absolute contig positions
void analyzeContigOverlapGraph(AnalyzeContext *actx);
int cmpOVLreadById(const void *a, const void *b);
int parseDBFileNames(char *dbName, int **fileOffsets, char ***fileNames, char ***readNames);
void classifyContigsByBReadsAndPath(AnalyzeContext *actx);
typedef struct
{
	int idx;
	int len;

} ContigIDAndLen;
int cmpContigLength(const void *a, const void *b);
int getPathID(AnalyzeContext *actx, int contigID);
void getContigsEndIDs(AnalyzeContext *actx, int contigID, int* beg, int *end);

//////// analyze contigs on their alignments with each other
void chainContigOverlaps(AnalyzeContext* ctx, Overlap* ovls, int n);
int analyzeChains(AnalyzeContext *ctx);
int processContigOverlap_handler(void* _ctx, Overlap* ovls, int novl);
int compareInt(const void * a, const void * b);
int isOverlapGroupValid(AnalyzeContext *actx, OverlapGroup *ovlgrp);
void addOverlapGroupToContig(AnalyzeContext *actx, OverlapGroup *ovlgrp, int symmetricAdd);
int contained(int ab, int ae, int bb, int be);
void updateOverlapGroup(AnalyzeContext *actx, OverlapGroup *ovlgrp, Overlap *pOvls, int *ovlIdx, int numOvlIdx);
void classifyContigsByOverlaps(AnalyzeContext *actx);
OverlapGroup* copySymmetricOverlapGroup(OverlapGroup *ovlgrp);

/// use all info from classifyContigsByOverlaps and classifyContigsByBReadsAndPath and tracks
/// to make final decision about contig class
void classify(AnalyzeContext *actx);
void writeFasta(AnalyzeContext *actx, Contig* contig, FILE* contigFile, int splitIdx, int cBegPos, int cEndPos);

void analyzeCoveredBases(AnalyzeContext* actx, Overlap* pOvls, int n);
//char getClassification(FinalContigClass fcc);
//char* getSplitType(SplitEvent* se);
void analyzeJunction(AnalyzeContext* actx, int read, int numContigs);
void createSplitEvent(AnalyzeContext * actx, Contig *contig, int read, int sType);

int checkJunctionCoverage(AnalyzeContext *actx, Contig *contig, int read); // returns true if everything is ok, false otherwise

int getRepeatBasesFromInterval(AnalyzeContext* actx, int readID, int beg, int end);
#endif /* UTILS_LAANALYZE_H_ */
