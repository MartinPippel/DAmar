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

	int repeatBases;
	// overlap information from final patched-reads overlap graph, that went into touring
	int nOutReads;
	int nInReads;

	int nComReadsWithNextContigRead; // out

	float avgCov;
	int lowCov;    // yes or no, indicates a problem e.g. chimeric read

	OvlRead *ovlReads;
	int numOvlReads;
	int maxOvlReads;

} ContigRead;

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
} ReadRelation;

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
	int flag;
	int corContigIdx;
	int *abpos;
	int *aepos;
	int *bbpos;
	int *bepos;
	int numPos;
} ContigRelation;

typedef struct
{
	int contigID1;
	int contigID2;	// only necessary if a contig bridges 2 other contigs
	int flag;
} TourRelation;

// CONTIG RELATION FLAGS
/// TourRelationFlags - based on touring evaluates ends arguments
#define REL_TOUR_PRIMARY						(1 << 0)
#define REL_TOUR_IS_ALT							(1 << 1) // contig is bubble or spur
#define REL_TOUR_HAS_ALT						(1 << 2) // contig is bubble or spur
#define REL_TOUR_UNIQUE							(1 << 3) // contig has no other tour relationship
#define REL_TOUR_IS_BRIDGE					(1 << 4) // contig bridges two other contigs
#define REL_TOUR_HAS_BRIDGE					(1 << 5) // contig has another contig that is a bridge contig
/// ContigRelationFlags - based on Contig vs Contig overlaps

#define REL_CONTIG_IS_ALT						(1 << 6)
#define REL_CONTIG_IS_REPEAT_ALT		(1 << 7)
#define REL_CONTIG_HAS_ALT 					(1 << 8)
#define REL_CONTIG_UNIQUE 					(1 << 9)
#define REL_CONTIG_BRIDGE 					(1 << 10)
/// ReadRelationFlags - based on Contig vs patched read overlaps

#define REL_READ_IS_ALT							(1 << 12)
#define REL_READ_HAS_ALT 						(1 << 13)
#define REL_READ_BRIDGE 						(1 << 14)
#define REL_READ_UNIQE	 						(1 << 15)    // 32768


// CONTIG CLASSIFICATION FLAGS
#define CLASS_CONTIG_CLASSIFIED					(1 << 0)
#define CLASS_CONTIG_PRIMARY						(1 << 1)
#define CLASS_CONTIG_ALT								(1 << 2)
#define CLASS_CONTIG_DISCARD						(1 << 3)
#define CLASS_CONTIG_DISCARD_REPEAT			(1 << 4)
#define CLASS_CONTIG_DISCARD_LEN				(1 << 5)
#define CLASS_CONTIG_DISCARD_CREADS			(1 << 6)
#define CLASS_CONTIG_MODIFIED_FALSEJOIN	(1 << 7)		// i.e. cut out corresponding Cread and keep classification, blacklist C-read!!!!
#define CLASS_CONTIG_MODIFIED_TRIMLEFT	(1 << 8)		// i.e. trim position: find trim_ab of corresponding Cread -> left contig part becomes ALT right part PRIM
#define CLASS_CONTIG_MODIFIED_TRIMRIGHT	(1 << 9)		// i.e. trim position: trim_ae of corresponding Cread -> left contig part stays PRIM right part becomes ALT
/// READ JUNCTION FLAGS // TODO in development
#define CLASS_READ_VALID     						(1 << 10)  	// 	1024  initially all reads are valid, because they are part of a contig
#define CLASS_READ_CHECKED	 						(1 << 11)  	// 	2048
#define CLASS_READ_DISCARD   						(1 << 12) 	// 	4096  contig read should be completely cut out from contig
#define CLASS_READ_NEWEND    						(1 << 13)  	// 	8192	contig ends at current read at position contigPosEnd
#define CLASS_READ_NEWSTART  						(1 << 14)  	// 	16384	contig starts at current read at position: contigPosBeg

typedef struct
{
	int len;
	int contigID;
	int fileID; 	  // contigs are usually split into connected compounds and stored different files
	int pathID;     // pathID represents an individual contig ID within a connected compound

	int cflag;			// contig classification flags
	int rflag; 			// contig relation flags

	int repBasesFromContigLAS;
	int repBasesFromReadLAS;

	int tmp;
} ContigProperties;

typedef struct
{
	TourRelation 		*correspTourRelation;
	ContigRelation 	*correspContigRelation;
	ReadRelation 		*correspReadRelation;
} Classification;


typedef struct
{
  ContigProperties property;

	ContigRead *cReads;
	int numcReads;
	float avgCov;

	ContigRelation *contigRelations;
	int numContigRelations;
	int maxContigRelations;

	ReadRelation *readRelations;
	int numReadRelations;
	int maxReadRelations;

	TourRelation *tourRelations;
	int numTourRelations;

	Classification *classif;
} Contig;

typedef struct
{
	int numFileNames;

	char **fileNames;
	char **fastaNames;
	int  *fromDbIdx;
	int  *toDbIdx;		// inclusive
} FileNamesAndOffsets;

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
	HITS_TRACK *patchedReadRepeat_track;

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
	int nFuzzBases;
	int contByReads_CommonReadFraction; 	// how many raw reads (in percent) must a contained contig have in common with a primary contig
	int contByReads_CoveredLenPerc;			  // length of a contained contig (in percent) that is covered from a primary contig
	int contByContigs_CoveredLenPerc;

	int minPrimContigLen;
	int minPrimContigReads;
	int maxPrimContigRepeatPerc;

	FileNamesAndOffsets *contigFileNamesAndOffsets;
	FileNamesAndOffsets *rawReadFileNamesAndOffsets;
} AnalyzeContext;

void initAnalyzeContext(AnalyzeContext *actx);
void finalContigValidation(AnalyzeContext *actx);

//void printSplitEvent(Contig *contig, SplitEvent *split);
//void printFinalAltContigClassification(FILE *out, FinalContigAltClass facc);
//void printFinalContigClassification(FILE *out, FinalContigClass fcc);
int getFullReadPositionInContig(AnalyzeContext *actx, Contig* c, int readID, int *cbeg, int *cend);
ContigRead* getContigRead(Contig* contig, int readID);

char* getContigFastaFileName(AnalyzeContext *actx, Contig *contig);
int createOutDir(char *out);

//////// analyze contigs based on fixed-read overlaps and intersection of read ids
int processFixedReadOverlap_handler1(void* _ctx, Overlap* ovls, int novl);
int processFixedReadOverlap_handler2(void* _ctx, Overlap* ovls, int novl);
void getCorrespondingPositionInBRead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2);
void getCorrespondingPositionInARead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2);
/// do coverage analysis for contig-read-overlaps, set absolute contig positions
void analyzeContigOverlapGraph(AnalyzeContext *actx);
int cmpOVLreadById(const void *a, const void *b);
void parseDBFileNames(char *dbName, FileNamesAndOffsets *fileNamesAndOffsets);
void classifyContigsByBReadsAndPath(AnalyzeContext *actx);
typedef struct
{
	int idx;
	int len;

} ContigIDAndLen;
int cmpContigLength(const void *a, const void *b);
int getPathID(AnalyzeContext *actx, int contigID);
int getFileID(FileNamesAndOffsets *fileAndNameOffset, int contigID);
int getNumberOfSubgraphContigs(FileNamesAndOffsets *fileAndNameOffset, int contigID);

void getContigsEndIDs(AnalyzeContext *actx, int contigID, int* beg, int *end);

//////// analyze contigs on their alignments with each other
int processContigOverlap_handler(void* _ctx, Overlap* ovls, int novl);
int compareInt(const void * a, const void * b);
int contained(int ab, int ae, int bb, int be);
void classifyContigsByOverlaps(AnalyzeContext *actx);

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

int getRepeatBasesFromInterval(AnalyzeContext* actx, int contigDB, int readID, int beg, int end);
#endif /* UTILS_LAANALYZE_H_ */
