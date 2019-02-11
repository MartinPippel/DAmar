/*
 * analyzeContigs.c
 *
 *  assign contigs into following sets: haploid, alternative, repeat
 *
 *  use:
 *  		self alignment information from contigs
 *  		union of fixed reads that map to the contigs (derived from fixedReadOverlaps)
 *
 *
 *  Created on: 3 Aug 2017
 *      Author: pippel
 */

#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/pass.h"
#include "lib/tracks.h"
#include "lib/utils.h"

#include "lib.ext/types.h"
#include "lib.ext/bitarr.h"
#include "utils/CTanalyze.h"

#undef DEBUG_STEP1A
#undef DEBUG_STEP1B
#define DEBUG_STEP1C

#define DEBUG_STEP2A

#define DEBUG
#define DEBUG_READ_INTERSECTION
#define DEBUG_CONTIG_COVERAGE
#define DEBUG_CONTIGOVERLAP

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define DEF_REPEATS_TRACK "repeats"

#define MIN_OVL_LOOKAHEAD 5000
#define STRIDE_OVL_LOOKAHEAD 5000
#define MAX_OVL_LOOKAHEAD 30000    // todo actually should be adapted according to underlying repeat length

#define	DEF_SPUR_LEN 50000
#define DEF_TIP_LEN 100000
#define DEF_ARG_F 10000
#define DEF_ARG_CRF 50
#define DEF_ARG_CL 50

#define DEF_ARG_L	100000 		// 	... minimum contig length, to be considered as a primary contig
#define DEF_ARG_N	5					//	... minimum number of contig reads,  to be considered as a primary contig
#define DEF_ARG_P	100				// 	... maximum repeat percentage, to be considered as a primary contig, not considered

/// defines for contig vs contig alignments
#define ConVsConMinAlign 5000
#define ConVsConAnchorBases 2000
#define ConVsConContainmentThreshold 0.7 // contig is contained if 70% of the bases are covered by another larger contig
#define ConVsConShortContigLen 100000
#define ConVsConShortContigContainmentThreshold 0.5

ContigFlag2Label cflag2label[] =
{
{ CLASS_CONTIG_DISCARD, "DIS", 'D' },
{
 CLASS_CONTIG_ALT, "ALT", 'A' },
{ CLASS_CONTIG_CLASSIFIED, "CLAS", 'C' },
{ CLASS_CONTIG_MODIFIED, "MOD", 'M' },
{
CLASS_CONTIG_PRIMARY, "PRI", 'P' },
{ CLASS_CONTIG_VISITED, "VIS", 'V' },
{ 0, NULL, '\0' } };

void cflags2str(char* pc, int flags)
{
	int i;
	int last = -1;
	for (i = 0; cflag2label[i].mask; i++)
	{
		if (flags & cflag2label[i].mask)
		{
			pc[i] = cflag2label[i].indicator;
			last = i;
		}
		else
		{
			pc[i] = ' ';
		}
	}

	pc[last + 1] = '\0';
}

static void printContigClassification(FILE* out, int flags, char sep)
{
	int i;
	for (i = 0; cflag2label[i].mask; i++)
	{
		if (flags & cflag2label[i].mask)
		{
			fprintf(out, "%s%c", cflag2label[i].label, sep);
		}
	}
}

//compare contigs by tmp flag, smallest first
static int cmp_contig_bytmp(const void* a, const void* b)
{
	Contig* c1 = *(Contig**) a;
	Contig* c2 = *(Contig**) b;

	return (c1->property.tmp - c2->property.tmp);
}

//compare contigs by length, longest first
static int cmp_contig_length(const void* a, const void* b)
{
	Contig* c1 = *(Contig**) a;
	Contig* c2 = *(Contig**) b;

	int cmp = c2->property.len - c1->property.len;

	if (!cmp)
	{
		cmp = (c2->avgCov - c1->avgCov);
	}

	return cmp;
}

void initAnalyzeContext(AnalyzeContext *actx)
{
	int i, j, k, l;
	int numOfContigs = DB_NREADS(actx->corContigDB);

	if (actx->VERBOSE)
		printf("numOfContigs: %d\n", numOfContigs);

	printf("allocate contig buffer");
	Contig *contigs = (Contig*) malloc(sizeof(Contig) * numOfContigs);
	actx->numContigs = numOfContigs;
	bzero(contigs, sizeof(Contig) * numOfContigs);

	// patched read IDs and Beg,End position from contig
	track_anno* patchedRead_anno = actx->corContigPatchReadsAndPos_track->anno;
	track_data* patchedRead_data = actx->corContigPatchReadsAndPos_track->data;

	// corrected read IDs and Beg,End position from corrected contig
	track_anno* correctedRead_anno = actx->corContigCorrectReadsAndPos_track->anno;
	track_data* correctedRead_data = actx->corContigCorrectReadsAndPos_track->data;

	printf(" - done\n");

	printf("allocate VreadMask");
	fflush(stdout);
	actx->maxVReadMask = 100;
	actx->vreadMask = (int**) malloc(sizeof(int*) * actx->maxVReadMask);
	assert(actx->vreadMask != NULL);
	for (i = 0; i < actx->maxVReadMask; i++)
	{
		actx->vreadMask[i] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
		bzero(actx->vreadMask[i], sizeof(int) * DB_NREADS(actx->patchedReadDB));
	}
	printf(" - done\n");

	actx->readSeq = New_Read_Buffer(actx->corContigDB);
	bzero(actx->readSeq, DB_READ_MAXLEN(actx->corContigDB));

	// parse Contig File Names and Offset, must be done before getFileID is called the first time
	actx->contigFileNamesAndOffsets = malloc(sizeof(FileNamesAndOffsets));
	bzero(actx->contigFileNamesAndOffsets, sizeof(FileNamesAndOffsets));

	parseDBFileNames(actx->corContigDBName, actx->contigFileNamesAndOffsets);

	printf(" init contigs - START\n");
	for (i = 0; i < numOfContigs; i++)
	{
		// init contigs with reads
		Contig *contig = contigs + i;

		contig->property.len = DB_READ_LEN(actx->corContigDB, i);
		contig->property.contigID = i;
		contig->property.repBasesFromContigLAS = getRepeatBasesFromInterval(actx, 1, i, 0, contig->property.len);
		contig->property.pathID = getPathID(actx, i);
		contig->property.fileID = getFileID(actx->contigFileNamesAndOffsets, i);

		track_anno patched_rb, patched_re;
		track_anno corrected_rb, corrected_re;

		patched_rb = patchedRead_anno[i] / sizeof(track_data);
		patched_re = patchedRead_anno[i + 1] / sizeof(track_data);
		assert(patched_rb < patched_re);

		corrected_rb = correctedRead_anno[i] / sizeof(track_data);
		corrected_re = correctedRead_anno[i + 1] / sizeof(track_data);
		assert(corrected_rb < corrected_re);

		// get number of contig reads
		contig->numcReads = (int) ((patched_re - patched_rb) / 3);
		assert(contig->numcReads == (int) ((corrected_re - corrected_rb) / 3));

		if (contig->numcReads == 0)
		{
			fprintf(stderr, "!strange! - no reads found for contig: %d. DISCARDED\n", i);
			contig->property.cflag |= CLASS_CONTIG_DISCARD;
			continue;
		}

		contig->numReadRelations = 0;
		contig->maxReadRelations = 1;
		contig->readRelations = (ReadRelation*) malloc(sizeof(ReadRelation) * contig->maxReadRelations);

		contig->numContigRelations = 0;
		contig->maxContigRelations = 0;
		contig->contigRelations = NULL;

		// init contig read buffer for overlapping reads
		contig->cReads = (ContigRead*) malloc(sizeof(ContigRead) * contig->numcReads);

		// init position buffer
		if (contig->cReads == NULL)
		{
			fprintf(stderr, "[ERROR] - ananlyzeContigs: Unable to allocate read buffers for Contig %d!\n", i);
			exit(1);
		}
		bzero(contig->cReads, sizeof(ContigRead) * contig->numcReads);

		int cov = MAX((int )(1.5 * actx->exp_cov), 30);
		for (j = 0, k = patched_rb; j < contig->numcReads; j++, k += 3)
		{
			ContigRead *cread = contig->cReads + j;

			cread->type |= CLASS_READ_VALID;
			cread->patchedID = patchedRead_data[k];
			cread->correctedID = correctedRead_data[k];
			if (j)
			{
				cread->patchedContigPosBeg = contig->cReads[j - 1].patchedContigPosEnd;
				cread->patchedContigPosEnd = cread->patchedContigPosBeg + abs(patchedRead_data[k + 2] - patchedRead_data[k + 1]);

				cread->correctedContigPosBeg = contig->cReads[j - 1].correctedContigPosEnd;
				cread->correctedContigPosEnd = cread->correctedContigPosBeg + abs(correctedRead_data[k + 2] - correctedRead_data[k + 1]);
			}
			else
			{
				cread->patchedContigPosBeg = 0;
				cread->patchedContigPosEnd = abs(patchedRead_data[k + 2] - patchedRead_data[k + 1]);

				cread->correctedContigPosBeg = 0;
				cread->correctedContigPosEnd = abs(correctedRead_data[k + 2] - correctedRead_data[k + 1]);
			}
			cread->correctedReadPosBeg = correctedRead_data[k + 1];
			cread->correctedReadPosEnd = correctedRead_data[k + 2];

			cread->ovlReads = (OvlRead*) malloc(cov * sizeof(OvlRead));
			if (cread->ovlReads == NULL)
			{
				fprintf(stderr, "[ERROR] - ananlyzeContigs: Unable to allocate OvlRead buffer for Contig %d at patchedRead %d!\n", i, patchedRead_data[k]);
				exit(1);
			}

			cread->maxOvlReads = cov;

			// add repeat bases from patched-reads repeat track to cread
			if (patchedRead_data[k + 1] < patchedRead_data[k + 2])
				cread->repeatBases = getRepeatBasesFromInterval(actx, 0, cread->patchedID, patchedRead_data[k + 1], patchedRead_data[k + 2]);
			else
				cread->repeatBases = getRepeatBasesFromInterval(actx, 0, cread->patchedID, patchedRead_data[k + 2], patchedRead_data[k + 1]);

			contig->property.repBasesFromReadLAS += cread->repeatBases;

			// insert cRead into ovlReads, cBeg and cEnd will be later translated into absolute contig positions!
			cread->ovlReads->patchedID = -patchedRead_data[k];
			cread->ovlReads->beg = patchedRead_data[k + 1];
			cread->ovlReads->end = patchedRead_data[k + 2];
			cread->ovlReads->cBeg = cread->patchedContigPosBeg;
			cread->ovlReads->cEnd = cread->patchedContigPosEnd;

			cread->numOvlReads = 1;


			// set vmaskread
			l = 0;
			while (l < actx->maxVReadMask && actx->vreadMask[l][patchedRead_data[k]] != 0)
				l++;

			if (l == actx->maxVReadMask)
			{
				printf("TERROR read %d is part of %d contigs!!!!", patchedRead_data[k], l);
				actx->maxVReadMask++;

				actx->vreadMask = (int**) realloc(actx->vreadMask, sizeof(int*) * actx->maxVReadMask);
				actx->vreadMask[l] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
				bzero(actx->vreadMask[l], sizeof(int) * DB_NREADS(actx->patchedReadDB));
			}
			// store negative 1-based contig id

			actx->vreadMask[l][patchedRead_data[k]] = -(i + 1);
		}

		if (actx->VERBOSE)
		{
			printf(
					"Contig_%d, reads (%d): [ unCorRead bpos epos len cBpos cEpos | corRead bpos epos len cBpos cEpos| corCumLen unCorCumLen | NReadRep PReadRep CumPReadRep ]\n",
					i, contig->numcReads);
			int corCumLen = 0;
			int unCorCumLen = 0;
			int unCorCumRepBases = 0;
			for (j = 0; j < contig->numcReads; j++)
			{
				unCorCumLen += contig->cReads[j].patchedContigPosEnd - contig->cReads[j].patchedContigPosBeg;
				corCumLen += contig->cReads[j].correctedContigPosEnd - contig->cReads[j].correctedContigPosBeg;
				unCorCumRepBases += contig->cReads[j].repeatBases;
				// uncorrected (patched) coordinates
				printf("%8d", contig->cReads[j].patchedID);
				printf(" %8d", contig->cReads[j].ovlReads[0].beg);
				printf(" %8d", contig->cReads[j].ovlReads[0].end);
				printf(" %5d", abs(contig->cReads[j].ovlReads[0].end - contig->cReads[j].ovlReads[0].beg));
				printf(" %8d", contig->cReads[j].patchedContigPosBeg);
				printf(" %8d", contig->cReads[j].patchedContigPosEnd);
				//  correx=cted coordinates
				printf(" |	 %8d", contig->cReads[j].correctedID);
				printf(" %8d", contig->cReads[j].correctedReadPosBeg);
				printf(" %8d", contig->cReads[j].correctedReadPosEnd);
				printf(" %5d", abs(contig->cReads[j].correctedReadPosEnd - contig->cReads[j].correctedReadPosBeg));
				printf(" %8d", contig->cReads[j].correctedContigPosBeg);
				printf(" %8d", contig->cReads[j].correctedContigPosEnd);

				printf(" | %8d", corCumLen);
				printf(" %8d", unCorCumLen);
				printf(" %5.3f%%", unCorCumLen * 100.0 / corCumLen);

				printf(" | %5d %2.f%% %.2f%%\n", contig->cReads[j].repeatBases,
						contig->cReads[j].repeatBases * 100.0 / (contig->cReads[j].patchedContigPosEnd - contig->cReads[j].patchedContigPosBeg),
						unCorCumRepBases * 100.0 / unCorCumLen);
			}
			printf(" --> len: %10d | fileID %4d | path %3d \n", contig->property.len, contig->property.fileID, contig->property.pathID);
		}

		//todo reduce this to actual number of relations
		int numContigsOfCurSubgraph = getNumberOfSubgraphContigs(actx->contigFileNamesAndOffsets, i);
		printf("numContigsOfCurSubgraph: %d\n", numContigsOfCurSubgraph);

		contig->tourRelations = (TourRelation*) malloc(sizeof(TourRelation) * numContigsOfCurSubgraph);
		bzero(contig->tourRelations, sizeof(TourRelation) * numContigsOfCurSubgraph);
	}

	printf("create touring relations\n");
	// create Touring relations
	for (i = 0; i < numOfContigs; i++)
	{
		Contig *contig = contigs + i;

		printf("contig id %d, len %d\n", contig->property.contigID, contig->property.len);
		int contigEndPathIdx1, contigEndPathIdx2;
		getContigsEndIDs(actx, i, &contigEndPathIdx1, &contigEndPathIdx2);
		printf("contigEndPathIdx1: %d, contigEndPathIdx2 %d \n", contigEndPathIdx1, contigEndPathIdx2);

		printf("contig->property.fileID: %d\n", contig->property.fileID);
		printf("actx->contigFileNamesAndOffsets->toDbIdx[contig->property.fileID] %d\n", actx->contigFileNamesAndOffsets->toDbIdx[contig->property.fileID]);
		printf("actx->contigFileNamesAndOffsets->fromDbIdx[contig->property.fileID]: %d\n", actx->contigFileNamesAndOffsets->fromDbIdx[contig->property.fileID]);

		fflush(stdout);
		int numContigsOfCurSubgraph = actx->contigFileNamesAndOffsets->toDbIdx[contig->property.fileID]
				- actx->contigFileNamesAndOffsets->fromDbIdx[contig->property.fileID] + 1;
		printf("numContigsOfCurSubgraph: %d\n", numContigsOfCurSubgraph);
		assert(MAX(contigEndPathIdx1, contigEndPathIdx2) <= numContigsOfCurSubgraph - 1);

		// 1. bubble
		if (contigEndPathIdx1 != -1 && contigEndPathIdx1 == contigEndPathIdx2)
		{
			int curPathID = contig->property.pathID;

			Contig *relatedContig = contigs + i - (curPathID - contigEndPathIdx1);

			assert(relatedContig != NULL);

			assert(contigEndPathIdx1 == relatedContig->property.pathID);

			contig->tourRelations[contig->numTourRelations].contigID1 = relatedContig->property.contigID;
			contig->tourRelations[contig->numTourRelations].flag |= REL_TOUR_IS_ALT;
			contig->numTourRelations++;

			relatedContig->tourRelations[relatedContig->numTourRelations].contigID1 = contig->property.contigID;
			relatedContig->tourRelations[relatedContig->numTourRelations].flag |= REL_TOUR_HAS_ALT;
			relatedContig->numTourRelations++;
		}

		// 2. spur
		else if ((contigEndPathIdx1 == -1 && contigEndPathIdx2 != -1) || (contigEndPathIdx1 != -1 && contigEndPathIdx2 == -1))
		{
			int curPathID = contig->property.pathID;
			int curPathRelID = MAX(contigEndPathIdx1, contigEndPathIdx2);
			Contig *relatedContig = contigs + i - (curPathID - curPathRelID);
			assert(relatedContig != NULL);
			assert(curPathRelID == relatedContig->property.pathID);

			contig->tourRelations[contig->numTourRelations].contigID1 = relatedContig->property.contigID;
			contig->tourRelations[contig->numTourRelations].flag |= REL_TOUR_IS_ALT;
			contig->numTourRelations++;

			relatedContig->tourRelations[relatedContig->numTourRelations].contigID1 = contig->property.contigID;
			relatedContig->tourRelations[relatedContig->numTourRelations].flag |= REL_TOUR_HAS_ALT;
			relatedContig->numTourRelations++;
		}

		// 3. link between two different contigs
		else if (contigEndPathIdx1 != -1 && contigEndPathIdx2 != -1)
		{
			int curPathID = contig->property.pathID;

			Contig *relatedContig1 = contigs + i - (curPathID - contigEndPathIdx1);
			assert(relatedContig1 != NULL);
			assert(contigEndPathIdx1 == relatedContig1->property.pathID);

			Contig *relatedContig2 = contigs + i - (curPathID - contigEndPathIdx2);
			assert(relatedContig2 != NULL);
			assert(contigEndPathIdx2 == relatedContig2->property.pathID);

			contig->tourRelations[contig->numTourRelations].contigID1 = relatedContig1->property.contigID;
			contig->tourRelations[contig->numTourRelations].contigID2 = relatedContig2->property.contigID;
			contig->tourRelations[contig->numTourRelations].flag |= REL_TOUR_IS_BRIDGE;
			contig->numTourRelations++;

			relatedContig1->tourRelations[relatedContig1->numTourRelations].contigID1 = contig->property.contigID;
			relatedContig1->tourRelations[relatedContig1->numTourRelations].contigID2 = relatedContig2->property.contigID;
			relatedContig1->tourRelations[relatedContig1->numTourRelations].flag |= REL_TOUR_HAS_BRIDGE;
			relatedContig1->numTourRelations++;

			relatedContig2->tourRelations[relatedContig2->numTourRelations].contigID1 = contig->property.contigID;
			relatedContig2->tourRelations[relatedContig2->numTourRelations].contigID2 = relatedContig1->property.contigID;
			relatedContig2->tourRelations[relatedContig2->numTourRelations].flag |= REL_TOUR_HAS_BRIDGE;
			relatedContig2->numTourRelations++;
		}
		else if(contigEndPathIdx1 == -1 && contigEndPathIdx2 == -1)
		{
			contig->tourRelations[contig->numTourRelations].contigID1 = -1;
			contig->tourRelations[contig->numTourRelations].flag |= REL_TOUR_PRIMARY;

			if(getNumberOfSubgraphContigs(actx->contigFileNamesAndOffsets, contig->property.contigID) == 1)
			{
				contig->tourRelations[contig->numTourRelations].flag |= REL_TOUR_UNIQUE;
			}
			contig->numTourRelations++;
		}
	}

	printf(" init contigs - DONE\n");

	actx->contigs = contigs;
	printf("sort contigs by length - START\n");
	// sort contig pointer according length (longest first)
	actx->contigs_sorted = malloc(sizeof(Contig*) * actx->numContigs);
	assert(actx->contigs_sorted != NULL);
	for (i = 0; i < actx->numContigs; i++)
	{
		actx->contigs_sorted[i] = actx->contigs + i;
	}
	qsort(actx->contigs_sorted, actx->numContigs, sizeof(Contig*), cmp_contig_length);
	printf("sort contigs by length - DONE\n");
}

char* getContigFastaFileName(AnalyzeContext *actx, Contig *contig)
{
	int i;

	for (i = 0; i < actx->contigFileNamesAndOffsets->numFileNames; i++)
	{
		if (contig->property.contigID >= actx->contigFileNamesAndOffsets->fromDbIdx[i] && contig->property.contigID < actx->contigFileNamesAndOffsets->fromDbIdx[i])
		{
			return actx->contigFileNamesAndOffsets->fileNames[i];
		}
	}

	fprintf(stderr, "[ERROR] - Cannot find contig file name of contig %d\n", contig->property.contigID);
	exit(1);
}

void parseDBFileNames(char *dbName, FileNamesAndOffsets *fileNamesAndOffsets)
{
	// parse old database file and store bax names as well as their read offsets

	assert(fileNamesAndOffsets != NULL);

	FILE *f;
	int numFiles = 0;
	char name[MAX_NAME];

	f = fopen(dbName, "r");
	if (f == NULL) // try to add .db extension
	{
		sprintf(name, "%s.db", dbName);
		f = fopen(name, "r");
	}
	if (f == NULL)
	{
		fprintf(stderr, "[ERROR] parseDBFileNames: Cannot open file %s for reading\n", dbName);
		exit(1);
	}

	if (fscanf(f, DB_NFILE, &numFiles) != 1)
	{
		fprintf(stderr, "[ERROR] parseDBFileNames: Cannot parse \"%s\" from file %s.\n",
		DB_NFILE, dbName);
		fclose(f);
		exit(1);
	}

	fileNamesAndOffsets->numFileNames = numFiles;

	fileNamesAndOffsets->fromDbIdx = (int *) malloc(sizeof(int) * (numFiles));
	memset(fileNamesAndOffsets->fromDbIdx, 0, numFiles);
	fileNamesAndOffsets->toDbIdx = (int *) malloc(sizeof(int) * (numFiles));
	memset(fileNamesAndOffsets->toDbIdx, 0, numFiles);

	fileNamesAndOffsets->fileNames = (char **) malloc(sizeof(char*) * (numFiles));
	fileNamesAndOffsets->fastaNames = (char **) malloc(sizeof(char*) * (numFiles));

	int i;
	for (i = 0; i < numFiles; i++)
	{
		fileNamesAndOffsets->fileNames[i] = (char*) malloc(sizeof(char) * MAX_NAME);
		fileNamesAndOffsets->fastaNames[i] = (char*) malloc(sizeof(char) * MAX_NAME);
	}

	int prevOffset, curOffset;
	prevOffset = 0;
	for (i = 0; i < numFiles; i++)
	{
		if (fscanf(f, DB_FDATA, &curOffset, fileNamesAndOffsets->fileNames[i], fileNamesAndOffsets->fastaNames[i]) != 3)
		{
			fprintf(stderr, "[ERROR] parseDBFileNames: Cannot parse \"%s\" in line %d from file %s.\n",
			DB_FDATA, i, dbName);
			fclose(f);
			exit(1);
		}

		assert(curOffset > prevOffset);
		fileNamesAndOffsets->fromDbIdx[i] = prevOffset;
		fileNamesAndOffsets->toDbIdx[i] = curOffset - 1;
		printf("fileNamesAndOffsets[%d]: %d %d\n", i, fileNamesAndOffsets->fromDbIdx[i], fileNamesAndOffsets->toDbIdx[i]);
		prevOffset = curOffset;
	}
}

static void pre(PassContext* pctx, AnalyzeContext* actx)
{
	actx->twidth = pctx->twidth;
}

ContigRead* getContigRead(Contig* contig, int readID)
{
	int i;

	for (i = 0; i < contig->numcReads; i++)
	{
		if (contig->cReads[i].patchedID == readID)
		{
			return contig->cReads + i;
		}
	}
	return NULL;
}

int getFullReadPositionInContig(AnalyzeContext *actx, Contig* contig, int readID, int *cbeg, int *cend)
{
	int rbeg;
	int tbeg, tend;
	get_trim(actx->patchedReadDB, actx->patchedReadTrim_track, readID, &tbeg, &tend);

	int i, j;
	int readPos = -1;

	for (i = 0; i < contig->numcReads; i++)
	{
		if (contig->cReads[i].patchedID == readID)
		{
			readPos = i;
			break;
		}
	}

	assert(readPos >= 0 && readPos < contig->numcReads);

	ContigRead *cread = contig->cReads + readPos;
//	printf("reads %d %d %d %d %d cIdx: %d\n", reads->id, reads->cBeg, reads->cEnd, reads->beg, reads->end, readPos);

	assert(abs(cread->patchedID) == readID);

	if (cread->ovlReads->beg < cread->ovlReads->end) // forward orientation
	{
		*cend = cread->ovlReads->cEnd;
		assert(cread->ovlReads->end == tend);
		*cbeg = cread->ovlReads->cBeg;
		if (cread->ovlReads->beg != tbeg) // try to find overlaps for readID with previous cReads
		{
			i = 1;
			while (readPos - i >= 0)
			{
//				printf("try to find aread %d in overlaps with previous Cread %d at position %d\n", readID, neighborCread, readPos - i);
				cread = contig->cReads + (readPos - i);
//				printf("cur cread %d %d %d %d %d\n", reads->id, reads->cBeg, reads->cEnd, reads->beg, reads->end);
				int found = 0;
				for (j = 1; j < cread->numOvlReads; j++)
				{
//					printf("check %d c(%d,%d) r(%d, %d)\n", reads[j].id, reads[j].cBeg, reads[j].cEnd, reads[j].beg, reads[j].end);
					if (cread->ovlReads[j].patchedID == readID)
					{
//						printf("found\n");
						assert(*cbeg > cread->ovlReads[j].cBeg);
						*cbeg = cread->ovlReads[j].cBeg;
						found = 1;

						// i.e. previous cread is in forward orientation
						if (((cread->ovlReads->beg < cread->ovlReads->end)))
						{
							rbeg = cread->ovlReads[j].beg;
							if (cread->ovlReads[j].beg == tbeg)
								found = 2;
						}
						if (cread->ovlReads->beg > cread->ovlReads->end)
						{
							rbeg = cread->ovlReads[j].end;
							if (cread->ovlReads[j].end == tbeg)
								found = 2;
						}
						break; // found start of current readID in contig
					}
				}
				if (found == 2)
				{
//					printf("iteration %d: found proper start of aread. i.e. proper cut position\n", i);
					break;
				}
				if (found == 0)
				{
					if (abs(rbeg) - tbeg < 1000) // i.e. there is no overlap, due to the minimum overlap length of daligner
					{
						*cbeg = (*cbeg) - (abs(rbeg) - tbeg);
						assert((*cbeg) >= 0 && (*cbeg) < (*cend) && (*cend) <= contig->property.len);
					}
					else
					{
						fprintf(stderr, "TERROR: iteration %d: UNABLE to find proper start of aread. i.e. unclear cut position. Remaining bases: %d\n", i,
								abs(rbeg) - tbeg);
					}
					break;
				}
				i++;
			}
		}
	}
	else // reverse complement orientation
	{
		*cend = cread->ovlReads->cEnd;
		assert(cread->ovlReads->end == tbeg);
		*cbeg = cread->ovlReads->cBeg;
		if (cread->ovlReads->beg != tend) // try to find overlaps for readID with previous cReads
		{
			i = 1;
			while (readPos - i >= 0)
			{
//				int neighborCread = contig->contigReadIds[readPos - i];
//				printf("try to find aread %d in overlaps with previous Cread %d at position %d\n", readID, neighborCread, readPos - i);
				cread = contig->cReads + (readPos - i);
//				printf("cur cread %d %d %d %d %d\n", reads->id, reads->cBeg, reads->cEnd, reads->beg, reads->end);
				int found = 0;
				for (j = 1; j < cread->numOvlReads; j++)
				{
//					printf("check %d c(%d,%d) r(%d, %d)\n", reads[j].id, reads[j].cBeg, reads[j].cEnd, reads[j].beg, reads[j].end);
					if (cread->ovlReads[j].patchedID == readID)
					{
//						printf("found\n");
						assert(*cbeg > cread->ovlReads[j].cBeg);
						*cbeg = cread->ovlReads[j].cBeg;

						found = 1;

						// i.e. previous cread is in forward orientation
						if (cread->ovlReads->beg < cread->ovlReads->end)
						{
							rbeg = cread->ovlReads[j].beg;
							if (cread->ovlReads[j].beg == tend)
								found = 2;
						}

						if (cread->ovlReads->beg > cread->ovlReads->end)
						{
							rbeg = cread->ovlReads[j].end;
							if (cread->ovlReads[j].end == tend)
								found = 2;
						}

						break; // found start of current readID in contig
					}
				}
				if (found == 2)
				{
//					printf("iteration %d: found proper start of aread. i.e. proper cut position\n", i);
					break;
				}
				if (found == 0)
				{
					if (tend - abs(rbeg) < 1000) // i.e. there is no overlap, due to the minimum overlap length of daligner
					{
						*cbeg = (*cbeg) - (tend - abs(rbeg));
						assert((*cbeg) >= 0 && (*cbeg) < (*cend) && (*cend) <= contig->property.len);
					}
					else
					{
						fprintf(stderr, "TERROR: iteration %d: UNABLE to find proper start of aread. i.e. unclear cut position. Remaining bases: %d\n", i,
								tend - abs(rbeg));
					}
					break;
				}
				i++;
			}
		}
	}

	return 0;
}

//void printSplitEvent(Contig *contig, SplitEvent *split)
//{
//	printf("CONTIG %d len %d new SplitEvent [%s] duplRead %d absPos (%d, %d, %d) leftRead (%d, %d,%d) rightRead (%d, %d,%d)\n", contig->property.contigID, contig->property.len,
//			getSplitType(split), split->duplicRead, split->type, split->contigPos1, split->contigPos2, split->leftProperRead, split->leftProperReadPos1,
//			split->leftProperReadPos2, split->rightProperRead, split->rightProperReadPos1, split->rightProperReadPos2);
//}

void analyzeJunction(AnalyzeContext* actx, int read, int numContigs)
{
	printf("analyzeJunction(AnalyzeContext* actx, int read %d , int numContigs %d)\n", read, numContigs);
	if (numContigs < 2)
		return;

	int i, j;

	// first do the easy cases:
	if (numContigs == 2)
	{
		// 1. check classification of both contigs
		// 1a) haploid vs contained todo check if we hav

		Contig *contig_i, *contig_j;
		contig_i = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		contig_j = actx->contigs + abs(actx->vreadMask[1][read]) - 1;

		// check if we have a normal primary and secondary contig connection
		if (contig_i->property.len < contig_j->property.len)
		{
			contig_i = actx->contigs + abs(actx->vreadMask[1][read]) - 1;
			contig_j = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		}

		// check putative heterozygous bubbles and spurs and their coverage (switch contents if necessary)
		// i.e.
		//						A1     								H1
		//    i.e.			 _______  						     ________
		//   		H1	____/	    \________ H1    OR    H1 ____/        \_____ H1
		//   				\_______/							   |____	/
		//                      H1        							A1
		// check haploid coverage
		int path_j_beg, path_j_end;
		getContigsEndIDs(actx, contig_j->property.contigID, &path_j_beg, &path_j_end);

		int containedByContigOverlaps = 0;
//		for (i = 0; i < contig_j->numOvlGrps; i++)
//		{
//			OverlapGroup *olg = contig_j->ovlGrps[i];
//			if ((olg->bread == contig_i->property.contigID) && (olg->flag & OVLGRP_AREAD_IS_CONTAINED))
//			{
//				containedByContigOverlaps = 1;
//				break;
//			}
//		}

		int containedByFixedReadCoverage = 0;
		for (i = 0; i < contig_j->numReadRelations; i++)
		{
			ReadRelation *readRel = contig_j->readRelations + i;
			if ((readRel->correspID == contig_i->property.contigID) && (readRel->flag & REL_CONTIG_IS_ALT))
			{
				containedByFixedReadCoverage = 1;
				break;
			}
		}

		// check coverage along junction
		checkJunctionCoverage(actx, contig_i, read);
		checkJunctionCoverage(actx, contig_j, read);

	}
	return;

	int cutOut = 0;
	//////////////////////////////////////////////////////////////////////// MULTIJOIN - EXACTLY 3 CONTIGS ////////////////////////////////////////////////////////////////////////
	if (numContigs == 3)
	{
		//ignore split if, two consecutive alt contigs collapse at same read in hap contig
		//						A1      H1 									H1         H2
		//    i.e.			 _______  ______								 ________  _________
		//   		H1	____/	    \/      \________ H1    OR    H1 ____/        \/         \_____ H2
		//   				\_______/\______/							   |____	/ \________/
		//                      H1       A2 									H2        A1

		Contig *contig_i, *contig_j, *contig_k;
		contig_i = contig_j = contig_k = NULL;

		contig_i = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		contig_j = actx->contigs + abs(actx->vreadMask[1][read]) - 1;
		contig_k = actx->contigs + abs(actx->vreadMask[2][read]) - 1;

		// all 3 contigs are haploid | or none of them is haploid
		if (((contig_i->property.cflag & CLASS_CONTIG_PRIMARY) && (contig_j->property.cflag & CLASS_CONTIG_PRIMARY)
				&& (contig_k->property.cflag & CLASS_CONTIG_PRIMARY))
				|| (!(contig_i->property.cflag & CLASS_CONTIG_PRIMARY) && !(contig_j->property.cflag & CLASS_CONTIG_PRIMARY)
						&& !(contig_k->property.cflag & CLASS_CONTIG_PRIMARY)))
		{
			cutOut = 1;
		}
		else
		{
			Contig *contig_hap1 = NULL;
			Contig *contig_hap2 = NULL;
			Contig *contig_alt1 = NULL;
			Contig *contig_alt2 = NULL;

			if (contig_i->property.cflag & CLASS_CONTIG_PRIMARY)
			{
				contig_hap1 = contig_i;

				if (contig_j->property.cflag & CLASS_CONTIG_PRIMARY)
				{
					contig_hap2 = contig_j;
					contig_alt1 = contig_k;
				}
				else if (contig_k->property.cflag & CLASS_CONTIG_PRIMARY)
				{
					contig_hap2 = contig_k;
					contig_alt1 = contig_j;
				}
				else
				{
					contig_alt1 = contig_j;
					contig_alt2 = contig_k;
				}
			}
			else if (contig_j->property.cflag & CLASS_CONTIG_PRIMARY)
			{
				contig_hap1 = contig_j;
				contig_alt1 = contig_i;
				if (contig_k->property.cflag & CLASS_CONTIG_PRIMARY)
				{
					contig_hap2 = contig_k;
				}
				else
				{
					contig_alt2 = contig_k;
				}
			}
			else if (contig_k->property.cflag & CLASS_CONTIG_PRIMARY)
			{
				contig_hap1 = contig_k;
				contig_alt1 = contig_i;
				contig_alt2 = contig_j;
			}
			else
			{
				fprintf(stderr, "SHIT happens\n");
				exit(1);
			}

			// check for first pattern
			if (contig_alt1 != NULL && contig_alt2 != NULL && contig_hap1 != NULL)
			{
				int a1End1, a1End2, a2End1, a2End2;
				getContigsEndIDs(actx, contig_alt1->property.contigID, &a1End1, &a1End2);
				getContigsEndIDs(actx, contig_alt2->property.contigID, &a2End1, &a2End2);

				if (
						(
							(((contig_alt1->property.cflag & CLASS_CONTIG_ALT)) && (a1End1 == contig_hap1->property.contigID || a1End2 == contig_hap1->property.contigID))
						)
						  && ((((contig_alt2->property.cflag & CLASS_CONTIG_ALT)) && (a1End1 == contig_hap1->property.contigID || a1End2 == contig_hap1->property.contigID)
						)

				   )
						//&& checkJunctionCoverage(actx, read, numContigs) // returns true if everything is ok, false otherwise
					)
				{
					// assign to both alt contigs a heterozygous split
					if (contig_alt1->cReads[0].patchedID == read && contig_alt2->cReads[contig_alt2->numcReads - 1].patchedID == read)
					{
//						createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
//						createSplitEvent(actx, contig_alt2, read, SPLIT_HETEROZYGRIGHT);
//						createSplitEvent(actx, contig_hap1, read, (SPLIT_HETEROZYGLEFT | SPLIT_HETEROZYGRIGHT));
					}
					else if (contig_alt2->cReads[0].patchedID == read && contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
					{
//						createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
//						createSplitEvent(actx, contig_alt2, read, SPLIT_HETEROZYGLEFT);
//						createSplitEvent(actx, contig_hap1, read, (SPLIT_HETEROZYGLEFT | SPLIT_HETEROZYGRIGHT));
					}
					else // always fall back to split behavior, i.e. cut out read from all contigs
					{
						cutOut = 1;
					}

				}
				else
				{
					cutOut = 1;
				}
			}
			// check for second pattern
			else if (contig_alt1 != NULL && contig_hap1 != NULL && contig_hap2 != NULL)
			{
				int end1, end2;
				getContigsEndIDs(actx, contig_alt1->property.contigID, &end1, &end2);

				// check if haploid contigs have dead end
				int deadEnd = 0;
				int contigCutPos1, contigCutPos2;
				contigCutPos1 = contigCutPos2 = -1;
				if (contig_hap1->cReads[0].patchedID == read || contig_hap1->cReads[contig_hap1->numcReads - 1].patchedID == read)
				{
					getFullReadPositionInContig(actx, contig_hap1, read, &contigCutPos1, &contigCutPos2);

					if (contigCutPos1 <= actx->TIP_LEN)
					{
						deadEnd = 1;
					}
					else if (contig_hap1->property.len - contigCutPos2 <= actx->TIP_LEN)
					{
						deadEnd = 2;
					}
				}
				else if (contig_hap2->cReads[0].patchedID == read || contig_hap2->cReads[contig_hap2->numcReads - 1].patchedID == read)
				{
					getFullReadPositionInContig(actx, contig_hap2, read, &contigCutPos1, &contigCutPos2);

					if (contigCutPos1 <= actx->TIP_LEN)
					{
						deadEnd = 1;
					}
					else if (contig_hap1->property.len - contigCutPos2 <= actx->TIP_LEN)
					{
						deadEnd = 2;
					}
				}

				if (deadEnd
						&& ((((contig_alt1->property.cflag & CLASS_CONTIG_ALT))
								&& (end1 == contig_hap1->property.contigID || end2 == contig_hap1->property.contigID || end1 == contig_hap2->property.contigID
										|| end2 == contig_hap2->property.contigID))
									)
						//&& checkJunctionCoverage(actx, read, numContigs) // returns true if everything is ok, false otherwise
						)
				{
					if (end1 == contig_hap1->property.contigID && end2 == contig_hap1->property.contigID) // ALT contig is a bubble of contig_hap1
					{
						if (contig_alt1->cReads[0].patchedID == read)
						{
//							createSplitEvent(actx, contig_hap1, read, SPLIT_HETEROZYGLEFT);
//							createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
						}
						else if (contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
						{
//							createsplitevent(actx, contig_hap1, read, split_heterozygright);
//							createsplitevent(actx, contig_alt1, read, split_heterozygright);
						}
						else
						{
							cutOut = 1;
						}
					}
					else if (end1 == contig_hap2->property.contigID && end2 == contig_hap2->property.contigID) // ALT contig is a bubble of contig_hap2
					{
						if (contig_alt1->cReads[0].patchedID == read)
						{
//							createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGLEFT);
//							createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
						}
						else if (contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
						{
//							createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGRIGHT);
//							createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
						}
						else
						{
							cutOut = 1;
						}
					}
					else if (contig_alt1->property.cflag == CLASS_CONTIG_ALT) // TODO --- this should be an ALT only  !!!!
					{
						if (contig_alt1->cReads[0].patchedID == read)
						{
							if (end1 == contig_hap1->property.contigID)
							{
//								createSplitEvent(actx, contig_hap1, read, SPLIT_HETEROZYGLEFT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
							}
							else if (end1 == contig_hap2->property.contigID)
							{
//								createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGLEFT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
							}
							else
							{
								cutOut = 1;
							}
						}
						else if (contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
						{
							if (end2 == contig_hap1->property.contigID)
							{
//								createSplitEvent(actx, contig_hap1, read, SPLIT_HETEROZYGRIGHT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
							}
							else if (end2 == contig_hap2->property.contigID)
							{
//								createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGRIGHT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
							}
							else
							{
								cutOut = 1;
							}
						}
						else
							cutOut = 1;
					}
					else
						cutOut = 1;

					if (!cutOut)
					{
						if (deadEnd == 1)
						{
//							createSplitEvent(actx, contig_hap1, read, SPLIT_DEADENDLEFT);
//							createSplitEvent(actx, contig_hap2, read, SPLIT_DEADENDLEFT);
						}
						else
						{
//							createSplitEvent(actx, contig_hap1, read, SPLIT_DEADENDRIGHT);
//							createSplitEvent(actx, contig_hap2, read, SPLIT_DEADENDRIGHT);
						}
					}

				}
				else
				{
					cutOut = 1;
				}

			}
			// split otherwise
			else
			{
				cutOut = 1;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////// MULTIJOIN - MORE THAN 3 CONTIGS ////////////////////////////////////////////////////////////////////////
	if (numContigs > 3 || (numContigs == 3 && cutOut == 1)) // thats really strange and should never happen
	{
		fprintf(stderr, "[WARNING] - Read %d is part of multiple contigs: ", read);

		for (i = 0; i < numContigs; i++)
		{
			printf(" %d", actx->vreadMask[i][read]);
		}
		printf("\n");
		assert(i == numContigs);

		j = 0;
		while (j < numContigs)
		{
			Contig *contig = actx->contigs + abs(actx->vreadMask[j][read]) - 1;
//			createSplitEvent(actx, contig, read, SPLIT_MULTIREAD);

		}
	}

	//////////////////////////////////////////////////////////////////////// READ IS PART OF EXACTLY 2 CONTIGS ////////////////////////////////////////////////////////////////////////
	///
	/// can be of type
	///
	/// 1. START_END      ------contig_i------X-----contig_j----, where X is read that is part of both contigs ---> always split those things
	///
	/// 2. BUBBLE/SPUR
	///
	///  			A)	ALT BUBBLE						b) ALT SPUR							C) DEADEND (SPUR or HAP and HAP)                  D) FALSE JOIN
	///
	///                   H1									H1										H1/H2
	///                _________							 __________								 ___________
	///      H1 ______/         \______ H1       H1 _____/          \________ H1			H1  ____/    		\______ H1/H2		     CONTIG_I_________X____________CONTIG_I
	///               \_________/						\_____|									\______|								    			      |
	///                   A1								   A1									   A1/H1											      |________ CONTIG_J
	else
	{

		// check START_END
		Contig *contig_j = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		Contig *contig_k = actx->contigs + abs(actx->vreadMask[1][read]) - 1;

		printf("else: %d vs %d\n", contig_j->property.contigID, contig_k->property.contigID);
		/////
		//// check the weird stuff first
		///
		// contig_k is circular
		if (contig_k->cReads[0].patchedID == read && contig_k->cReads[contig_k->numcReads - 1].patchedID == read)
		{
			fprintf(stderr, "-----------------------> circle in contig_k %d\n", contig_k->property.contigID);
			exit(1);
		}
		// contig_j is circular
		else if (contig_j->cReads[0].patchedID == read && contig_j->cReads[contig_j->numcReads - 1].patchedID == read)
		{
			fprintf(stderr, "-----------------------> circle in contig_j %d\n", contig_j->property.contigID);
			exit(1);
		}
		// contig_j stops at read X and contig_k starts at read X
		else if ((contig_k->cReads[0].patchedID == read || contig_k->cReads[contig_k->numcReads - 1].patchedID == read)
				&& (contig_j->cReads[contig_j->numcReads - 1].patchedID == read || contig_j->cReads[0].patchedID == read))
		{
			printf("// contig_j stops at read X and contig_k starst at read X\n");
//			createSplitEvent(actx, contig_j, read, SPLIT_STARTEND);
//			createSplitEvent(actx, contig_k, read, SPLIT_STARTEND);
		}
		else
		{
			//// first, get position of read X in both contigs
			int contigJCutPos1, contigJCutPos2;
			int contigKCutPos1, contigKCutPos2;

			contigJCutPos1 = contigJCutPos2 = contigKCutPos1 = contigKCutPos2 = -1;

			Contig *contig_j = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
			Contig *contig_k = actx->contigs + abs(actx->vreadMask[1][read]) - 1;

			getFullReadPositionInContig(actx, contig_j, read, &contigJCutPos1, &contigJCutPos2);
			getFullReadPositionInContig(actx, contig_k, read, &contigKCutPos1, &contigKCutPos2);

			assert(contigJCutPos1 != -1 && contigJCutPos2 != -1);
			assert(contigKCutPos1 != -1 && contigKCutPos2 != -1);

			printf("contig %d read %d at [%d, %d] contig %d read %d at [%d, %d]\n", contig_j->property.contigID, read, contigJCutPos1, contigJCutPos2,
					contig_k->property.contigID, read, contigKCutPos1, contigKCutPos2);
			//int type = 0;

			//int junctionCov = checkJunctionCoverage(actx, read, numContigs); // returns true if everything is ok, false otherwise

			// both contigs are haploid
			//if (contig_j->fClass == HAPLOID && contig_k->fClass == HAPLOID)
			if (1)
			{
				// check if there is a valuable dead end, otherwise its a false join

				// todo remove this after testing
				//	assert(contig_j->property.len > 200000 || contig_k->property.contigID > 200000);

				if (contigJCutPos1 == 0)
				{
					if (contigKCutPos1 > actx->TIP_LEN && contig_k->property.len - contigKCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("1 DEAD END ----> check number of read that follow both paths\n");
					}
				}
				else if (contigJCutPos2 == contig_j->property.len)
				{
					if (contigKCutPos1 > actx->TIP_LEN && contig_k->property.len - contigKCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("2 DEAD END ----> check number of read that follow both paths\n");
					}
				}

				else if (contigKCutPos1 == 0)
				{
					if (contigJCutPos1 > actx->TIP_LEN && contig_j->property.len - contigJCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("3 DEAD END ----> check number of read that follow both paths\n");
					}
				}
				else if (contigKCutPos2 == contig_k->property.len)
				{
					if (contigJCutPos1 > actx->TIP_LEN && contig_j->property.len - contigJCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("4 DEAD END ----> check number of read that follow both paths\n");

					}
				}

			}
			// only one contig is haploid
			else if ((contig_j->property.cflag & CLASS_CONTIG_PRIMARY) || (contig_k->property.cflag & CLASS_CONTIG_PRIMARY))
			{
				// TODO
			}
			// both contigs are NOT haploid
			else
			{
				// TODO
			}

			/*

			 /////
			 //// check the "normal" problems. Distinguish between
			 ///  1) normal heterozygous bubbles/spurs (i.e. clear ALT contigs)
			 ///  2) putative false joins
			 ///  3) as well as overlapping contigs in heterozygous regions // i.e. dead end
			 ///




			 int stype = 0;

			 if (contigJCutPos1 == 0)
			 {
			 if (contigKCutPos1 > actx->TIP_LEN && contig_k->property.len - contigKCutPos2 > actx->TIP_LEN) // well inside contigK
			 {
			 if (contig_j->property.len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGRIGHT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINRIGHT;
			 }
			 }
			 else // at tips of contigK --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDRIGHT;
			 }
			 }
			 else if (contigKCutPos1 == 0)
			 {
			 if (contigJCutPos1 > actx->TIP_LEN && contig_j->property.len - contigJCutPos2 > actx->TIP_LEN) // well inside contigJ
			 {
			 if (contig_k->property.len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGRIGHT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINRIGHT;
			 }
			 }
			 else // at tips of contigJ --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDRIGHT;
			 }
			 }

			 else if (contigJCutPos2 == contig_j->property.len)
			 {
			 if (contigKCutPos1 > actx->TIP_LEN && contig_k->property.len - contigKCutPos2 > actx->TIP_LEN) // well inside contigK
			 {
			 if (contig_j->property.len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGLEFT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINLEFT;
			 }
			 }
			 else // at tips of contigK --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDLEFT;
			 }
			 }
			 else if (contigKCutPos2 == contig_k->property.len)
			 {
			 if (contigJCutPos1 > actx->TIP_LEN && contig_j->property.len - contigJCutPos2 > actx->TIP_LEN) // well inside contigJ
			 {
			 if (contig_k->property.len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGLEFT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINLEFT;
			 }
			 }
			 else // at tips of contigK --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDLEFT;
			 }
			 }
			 else
			 {
			 printf("SHIT no classification\n");
			 exit(1);
			 }

			 printf("############## remaining duplicate Read: %d: contig_j %d len %d pos %d, %d TYPE: %d and contig_k %d len %d pos %d, %d TYPE %d\n",
			 aread, contig_j->property.contigID, contig_j->property.len, contigJCutPos1, contigJCutPos2, stype, contig_k->property.contigID, contig_k->property.len, contigKCutPos1,
			 contigKCutPos2, stype);

			 // add split position to contig_j
			 if (contig_j->maxSplits == contig_j->numSplits)
			 {
			 contig_j->maxSplits = (contig_j->maxSplits * 1.2) + 2;
			 contig_j->split = (SplitEvent*) realloc(contig_j->split, sizeof(SplitEvent) * contig_j->maxSplits);
			 assert(contig_j->split != NULL);
			 }
			 assert(contigJCutPos1 >= 0 && contigJCutPos1 < contigJCutPos2 && contigJCutPos2 <= contig_j->property.len);
			 {
			 // find contig reads and their positions for corresponding cutPositions
			 int leftProperRead, rightProperRead;
			 leftProperRead = rightProperRead = -1;
			 int leftProperReadPos1, leftProperReadPos2, rightProperReadPos1, rightProperReadPos2;
			 leftProperReadPos1 = leftProperReadPos2 = rightProperReadPos1 = rightProperReadPos2 = 0;
			 for (k = 0; k < contig_j->numContigReadIds; k++)
			 {
			 if (leftProperRead < 0 && contig_j->reads[k]->cEnd >= contigJCutPos1 && contigJCutPos1 > 0)
			 {
			 leftProperRead = abs(contig_j->reads[k]->id);
			 leftProperReadPos1 = contig_j->reads[k]->beg;
			 if (contig_j->reads[k]->beg > contig_j->reads[k]->end) // reverse complement orientation
			 {
			 leftProperReadPos2 = contig_j->reads[k]->end + (contig_j->reads[k]->cEnd - contigJCutPos1);
			 printf("RV %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 > leftProperReadPos2);
			 }
			 else // forward orientation
			 {
			 leftProperReadPos2 = contig_j->reads[k]->end - (contig_j->reads[k]->cEnd - contigJCutPos1);
			 printf("FW %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 < leftProperReadPos2);
			 }
			 }
			 if (contig_j->reads[k]->cEnd > contigJCutPos2 && contigJCutPos2 < contig_j->property.len)
			 {
			 rightProperRead = abs(contig_j->reads[k]->id);
			 rightProperReadPos2 = contig_j->reads[k]->end;
			 if (contig_j->reads[k]->beg > contig_j->reads[k]->end) // reverse complement orientation
			 {
			 rightProperReadPos1 = contig_j->reads[k]->end + (contig_j->reads[k]->cEnd - contigJCutPos2);
			 printf("FW %d c(%d, %d) r(%d,%d) right(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 > rightProperReadPos2);
			 }
			 else // forward orientation
			 {
			 rightProperReadPos1 = contig_j->reads[k]->end - (contig_j->reads[k]->cEnd - contigJCutPos2);
			 printf("RV %d c(%d, %d) r(%d,%d) right(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 < rightProperReadPos2);
			 }
			 break;
			 }
			 }

			 SplitEvent *split = contig_j->split + contig_j->numSplits;
			 split->duplicRead = aread; // TODO do we need thise aread info
			 split->contigPos1 = contigJCutPos1;
			 split->contigPos2 = contigJCutPos2;
			 split->type = stype;

			 split->leftProperRead = leftProperRead;
			 split->leftProperReadPos1 = leftProperReadPos1;
			 split->leftProperReadPos2 = leftProperReadPos2;

			 split->rightProperRead = rightProperRead;
			 split->rightProperReadPos1 = rightProperReadPos1;
			 split->rightProperReadPos2 = rightProperReadPos2;

			 printSplitEvent(contig_j, split);
			 contig_j->numSplits++;
			 }

			 // add split position to contig
			 if (contig_k->maxSplits == contig_k->numSplits)
			 {
			 contig_k->maxSplits = (contig_k->maxSplits * 1.2) + 2;
			 contig_k->split = (SplitEvent*) realloc(contig_k->split, sizeof(SplitEvent) * contig_k->maxSplits);
			 assert(contig_k->split != NULL);
			 }
			 assert(contigKCutPos1 >= 0 && contigKCutPos1 < contigKCutPos2 && contigKCutPos2 <= contig_k->property.len);
			 {
			 // find contig reads and their positions for corresponding cutPositions
			 int leftProperRead, rightProperRead;
			 leftProperRead = rightProperRead = -1;
			 int leftProperReadPos1, leftProperReadPos2, rightProperReadPos1, rightProperReadPos2;
			 leftProperReadPos1 = leftProperReadPos2 = rightProperReadPos1 = rightProperReadPos2 = 0;
			 for (k = 0; k < contig_k->numContigReadIds; k++)
			 {
			 if (leftProperRead < 0 && contig_k->reads[k]->cEnd >= contigKCutPos1 && contigKCutPos1 > 0)
			 {
			 leftProperRead = abs(contig_k->reads[k]->id);
			 leftProperReadPos1 = contig_k->reads[k]->beg;
			 if (contig_k->reads[k]->beg > contig_k->reads[k]->end) // reverse complement orientation
			 {
			 leftProperReadPos2 = contig_k->reads[k]->end + (contig_k->reads[k]->cEnd - contigKCutPos1);
			 printf("RV %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 > leftProperReadPos2);
			 }
			 else // forward orientation
			 {
			 leftProperReadPos2 = contig_k->reads[k]->end - (contig_k->reads[k]->cEnd - contigKCutPos1);
			 printf("FW %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 < leftProperReadPos2);
			 }
			 }
			 if (contig_k->reads[k]->cEnd > contigKCutPos2 && contigKCutPos2 < contig_k->property.len)
			 {
			 rightProperRead = abs(contig_k->reads[k]->id);
			 rightProperReadPos2 = contig_k->reads[k]->end;
			 if (contig_k->reads[k]->beg > contig_k->reads[k]->end) // reverse complement orientation
			 {
			 rightProperReadPos1 = contig_k->reads[k]->end + (contig_k->reads[k]->cEnd - contigKCutPos2);
			 printf("RV %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 > rightProperReadPos2);
			 }
			 else // forward orientation
			 {
			 rightProperReadPos1 = contig_k->reads[k]->end - (contig_k->reads[k]->cEnd - contigKCutPos2);
			 printf("FW %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 < rightProperReadPos2);
			 }
			 break;
			 }
			 }

			 SplitEvent *split = contig_k->split + contig_k->numSplits;
			 split->duplicRead = aread; // TODO do we need thise aread info
			 split->contigPos1 = contigKCutPos1;
			 split->contigPos2 = contigKCutPos2;
			 split->type = stype;

			 split->leftProperRead = leftProperRead;
			 split->leftProperReadPos1 = leftProperReadPos1;
			 split->leftProperReadPos2 = leftProperReadPos2;

			 split->rightProperRead = rightProperRead;
			 split->rightProperReadPos1 = rightProperReadPos1;
			 split->rightProperReadPos2 = rightProperReadPos2;

			 printSplitEvent(contig_k, split);
			 contig_k->numSplits++;
			 }*/

		}

	}
}

void createSplitEvent(AnalyzeContext * actx, Contig *contig, int read, int sType)
{
	int contigCutPos1, contigCutPos2;

	contigCutPos1 = contigCutPos2 = -1;

	getFullReadPositionInContig(actx, contig, read, &contigCutPos1, &contigCutPos2);

	assert(contigCutPos1 != -1 && contigCutPos2 != -1);

	// add split position to contig
//	if (contig->maxSplits == contig->numSplits)
//	{
//		contig->maxSplits = (contig->maxSplits * 1.2) + 2;
//		contig->split = (SplitEvent*) realloc(contig->split, sizeof(SplitEvent) * contig->maxSplits);
//		assert(contig->split != NULL);
//	}
	assert(contigCutPos1 >= 0 && contigCutPos1 < contigCutPos2 && contigCutPos2 <= contig->property.len);
	{
		// find contig reads and their positions for corresponding cutPositions
		int leftProperRead, rightProperRead;
		leftProperRead = rightProperRead = -1;
		int leftProperReadPos1, leftProperReadPos2, rightProperReadPos1, rightProperReadPos2;
		leftProperReadPos1 = leftProperReadPos2 = rightProperReadPos1 = rightProperReadPos2 = 0;

		///// TODO fix this
//      int k;
//		for (k = 0; k < contig->numContigReadIds; k++)
//		{
//			if (leftProperRead < 0 && contig->reads[k]->cEnd >= contigCutPos1 && contigCutPos1 > 0)
//			{
//				leftProperRead = abs(contig->reads[k]->id);
//				leftProperReadPos1 = contig->reads[k]->beg;
//				if (contig->reads[k]->beg > contig->reads[k]->end) // reverse complement orientation
//				{
//					leftProperReadPos2 = contig->reads[k]->end + (contig->reads[k]->cEnd - contigCutPos1);
//					assert(leftProperReadPos1 > leftProperReadPos2);
//				}
//				else // forward orientation
//				{
//					leftProperReadPos2 = contig->reads[k]->end - (contig->reads[k]->cEnd - contigCutPos1);
//					assert(leftProperReadPos1 < leftProperReadPos2);
//				}
//			}
//			if (contig->reads[k]->cEnd > contigCutPos2 && contigCutPos2 < contig->property.len)
//			{
//				rightProperRead = abs(contig->reads[k]->id);
//				rightProperReadPos2 = contig->reads[k]->end;
//				if (contig->reads[k]->beg > contig->reads[k]->end) // reverse complement orientation
//				{
//					rightProperReadPos1 = contig->reads[k]->end + (contig->reads[k]->cEnd - contigCutPos2);
//					assert(rightProperReadPos1 > rightProperReadPos2);
//				}
//				else // forward orientation
//				{
//					rightProperReadPos1 = contig->reads[k]->end - (contig->reads[k]->cEnd - contigCutPos2);
//					assert(rightProperReadPos1 < rightProperReadPos2);
//				}
//				break;
//			}
//		}

//		SplitEvent *split = contig->split + contig->numSplits;
//		split->duplicRead = read; // TODO do we need thise aread info
//		split->contigPos1 = contigCutPos1;
//		split->contigPos2 = contigCutPos2;
//		split->type = sType;
//
//		split->leftProperRead = leftProperRead;
//		split->leftProperReadPos1 = leftProperReadPos1;
//		split->leftProperReadPos2 = leftProperReadPos2;
//
//		split->rightProperRead = rightProperRead;
//		split->rightProperReadPos1 = rightProperReadPos1;
//		split->rightProperReadPos2 = rightProperReadPos2;
//
//		contig->numSplits++;
	}
}

void finalContigValidation(AnalyzeContext *actx)
{
	// reset CONTIG_READ_CHECKED flag as it will be reused
	int i, j;
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;
		for (j = 0; j < contig->numcReads; j++)
			contig->cReads[j].type &= ~( CLASS_READ_CHECKED);
	}

	int MIN_OUT_COV = 3;
	int MIN_COM_OUT_COV = 1;

	// first: check coverage drops in all contigs, and mark contig read as putative "BREAK"
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;
		int j;
		for (j = 0; j < contig->numcReads - 1; j++)
		{
			ContigRead *cread = contig->cReads + j;
			if (cread->nOutReads < MIN_OUT_COV || cread->nComReadsWithNextContigRead < MIN_COM_OUT_COV)
			{
				printf("CONTIG %d FOUND LOW  COVERAGE CP: read %d, # in %d, # out %d, #comRead %d, cPos [%d, %d], rPos [%d, %d], avgCov %.2f, covDrop %d\n",
						contig->property.contigID, cread->patchedID, cread->nInReads, cread->nOutReads, cread->nComReadsWithNextContigRead, cread->patchedContigPosBeg,
						cread->patchedContigPosEnd, cread->ovlReads->beg, cread->ovlReads->end, cread->avgCov, cread->lowCov);

				cread->type |= (CLASS_READ_BREAK_UNKNOWN);
			}
			cread->type |= (CLASS_READ_CHECKED);
		}
	}

	// todo really necessary?
	for (i = 0; i < actx->numContigs; i++)
		actx->contigs[i].property.cflag &= ~( CLASS_CONTIG_VISITED);

	int k;
	Contig **secondaryContigs;
	int maxSecondaryContigs = 10;
	int curSecondaryContigs = 0;
	secondaryContigs = malloc(sizeof(Contig*) * maxSecondaryContigs);
	assert(secondaryContigs != NULL);

	// second: full analysis of contigs, in each step consider all contigs that are somehow related to each other
	// 2. loop over length-sorted contigs until all contigs are visited:
	//        a) take next longest contig
	//        b) add all contigs that are related to 2a)
	//        c) analyze all relations of contigs in 2b): heterozygous, repeat induced, false joins, and split contigs if necessary
	for (i = 0; i < actx->numContigs; i++)
	{
		if (actx->contigs_sorted[i]->property.cflag & CLASS_CONTIG_VISITED)
		{
			continue;
		}

		Contig * primaryContig = actx->contigs_sorted[i];
		curSecondaryContigs = 0;

		// check all reads, that open a junction
		for (j = 0; j < primaryContig->numcReads; j++)
		{
			int aread = abs(primaryContig->cReads[j].patchedID);

			int nOcc = 0;
			// check only reads that were used in touring
			while (nOcc < actx->maxVReadMask && actx->vreadMask[nOcc][aread] < 0)
			{
				nOcc++;
			}

			if (nOcc > 1)
			{
				for (k = 0; k < nOcc; k++)
				{
					int contigIdx = abs(actx->vreadMask[k][aread]) - 1;
					if (contigIdx == primaryContig->property.contigID)
						continue;

					int add = 1;
					int l;
					for (l = 0; l < curSecondaryContigs; l++)
					{
						if (secondaryContigs[l]->property.contigID == contigIdx)
						{
							add = 0;
							break;
						}
					}
					if (add)
					{
						if (curSecondaryContigs == maxSecondaryContigs)
						{
							maxSecondaryContigs = (maxSecondaryContigs * 1.2 + 5);
							secondaryContigs = (Contig**) realloc(secondaryContigs, sizeof(Contig*) * maxSecondaryContigs);
							assert(secondaryContigs != NULL);
						}
						secondaryContigs[curSecondaryContigs] = actx->contigs + contigIdx;
						secondaryContigs[curSecondaryContigs]->property.tmp = primaryContig->cReads[j].patchedContigPosBeg;
						curSecondaryContigs++;
					}
				}
			}
		}
		// check for remaining repeat-induced contigs, that don't share a common contig read
		// check contig overlap groups
//		for (j = 0; j < primaryContig->numOvlGrps; j++)
//		{
//			if (primaryContig->ovlGrps[j]->flag & OVLGRP_DISCARD)
//				continue;
//
//			if (primaryContig->ovlGrps[j]->flag & OVLGRP_AREAD_IS_CONTAINED)
//			{
//				int add = 1;
//				for (k = 0; k < curSecondaryContigs; k++)
//				{
//					if (secondaryContigs[k]->property.contigID == primaryContig->ovlGrps[j]->bread)
//					{
//						add = 0;
//						break;
//					}
//				}
//				if (add)
//				{
//					if (curSecondaryContigs == maxSecondaryContigs)
//					{
//						maxSecondaryContigs = (maxSecondaryContigs * 1.2 + 5);
//						secondaryContigs = (Contig**) realloc(secondaryContigs, sizeof(Contig*) * maxSecondaryContigs);
//						assert(secondaryContigs != NULL);
//					}
//					secondaryContigs[curSecondaryContigs] = actx->contigs + primaryContig->ovlGrps[j]->bread;
//					secondaryContigs[curSecondaryContigs]->tmp = primaryContig->ovlGrps[j]->first_abpos;
//					curSecondaryContigs++;
//				}
//			}
//		}

		// check contig graph classification
		for (j = 0; j < primaryContig->numReadRelations; j++)
		{
			if (primaryContig->readRelations[j].flag & REL_CONTIG_HAS_ALT)
			{
				int add = 1;
				for (k = 0; k < curSecondaryContigs; k++)
				{
					if (secondaryContigs[k]->property.contigID == primaryContig->readRelations[j].correspID)
					{
						add = 0;
						break;
					}
				}
				if (add)
				{
					if (curSecondaryContigs == maxSecondaryContigs)
					{
						maxSecondaryContigs = (maxSecondaryContigs * 1.2 + 5);
						secondaryContigs = (Contig**) realloc(secondaryContigs, sizeof(Contig*) * maxSecondaryContigs);
						assert(secondaryContigs != NULL);
					}
					secondaryContigs[curSecondaryContigs] = actx->contigs + primaryContig->readRelations[j].correspID;
					secondaryContigs[curSecondaryContigs]->property.tmp = primaryContig->readRelations[j].coveredIntervals[0];
					curSecondaryContigs++;
				}
			}
		}

		// sort secondary contigs according to mapping position in primary contig
		qsort(secondaryContigs, curSecondaryContigs, sizeof(Contig*), cmp_contig_bytmp);

		printf("Primary Contig: %d LEN %d class %d \n", primaryContig->property.contigID, primaryContig->property.len, primaryContig->property.cflag);

		// analyze contig relations
		for (j = 0; j < curSecondaryContigs; j++)
		{
			Contig *secondaryContig = secondaryContigs[j];
			printf("Secondary Contig: %d LEN %d start at: %d class %d\n", secondaryContig->property.contigID, secondaryContig->property.len,
					secondaryContig->property.tmp, secondaryContig->property.cflag);
		}

		primaryContig->property.cflag |= CLASS_CONTIG_VISITED;
		for (j = 0; j < curSecondaryContigs; j++)
			secondaryContigs[j]->property.cflag |= CLASS_CONTIG_VISITED;
	}

	// clean up temp contig buffer
	free(secondaryContigs);

}

//void analyzeRepeatContent(AnalyzeContext *actx)
//{
//	int i;
//
//	for (i = 0; i < actx->numContigs; i++)
//	{
//		Contig *contig = actx->contigs + i;
//
////		printf("analyzeRepeatContent of contig %d (%d), len %d", i, contig->property.contigID, contig->property.len);
//
//		track_anno* ranno = (track_anno*) (actx->contigRepeats_track->anno);
//		track_data* rdata = (track_data*) (actx->contigRepeats_track->data);
//
//		track_anno ob = ranno[i];
//		track_anno oe = ranno[i + 1];
//
//		// no repeats, nothing to do
//
//		if (ob >= oe)
//		{
////			printf(" repeat 0 0 0 REALLY?\n");
//			continue;
//		}
//
//		int nrepeat = 0;
//		int repMerged = 0;
//		int repBases = 0;
//		int MERGE_BASES = 50;
//		int tipL = 0;
//		int tipR = 0;
//
//		ob /= sizeof(track_data);
//		oe /= sizeof(track_data);
//
//		track_anno numIntervals = oe - ob;
//		if (actx->VERBOSE > 2)
//			printf(" number of repeats: (%llu)\n", numIntervals);
//
//		int prvb = rdata[ob];
//		int prve = rdata[ob + 1];
//		ob += 2;
//
//		while (ob < oe)
//		{
//			int b = rdata[ob];
//			int e = rdata[ob + 1];
//
//			if (actx->VERBOSE > 2)
//				printf(" [%d, %d]\n", b, e);
//
//			if (b - prve < MERGE_BASES)
//			{
//				// merge intervals
//				if (actx->VERBOSE > 2)
//					printf("merge intervals : [%d,%d] [%d,%d]", prvb, prve, b, e);
//				prve = e;
//				repMerged++;
//				if (actx->VERBOSE > 2)
//					printf(" --> [%d,%d] (merged: %d)\n", prvb, prve, repMerged);
//			}
//			else
//			{
//				nrepeat++;
//				repBases += prve - prvb;
//				if (prvb < 25000)
//				{
//					tipL += MIN(25000, prve) - prvb;
//				}
//				if (prve > contig->property.len - 25000)
//				{
//					tipR += prve - MAX(prvb, contig->property.len - 25000);
//				}
//				prvb = b;
//				prve = e;
//
//			}
//			ob += 2;
//		}
//		nrepeat++;
//		repBases += prve - prvb;
//
//		if (prvb < 25000)
//		{
//			tipL += MIN(25000, prve) - prvb;
//		}
//		if (prve > contig->property.len - 25000)
//		{
//			tipR += prve - MAX(prvb, contig->property.len - 25000);
//		}
//
//		contig->repBases = repBases;
//		contig->repBasesTipLeft = tipL;
//		contig->repBasesTipRight = tipR;
//		contig->numRepeats = nrepeat;
////		printf(" repeat %d %d %d (%.2f%%) tipL %d, tipR %d\n", nrepeat, repBases, repMerged, repBases * 100.0 / contig->property.len, tipL, tipR);
//	}
//}

void getCorrespondingPositionInARead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2)
{
	int cBeg = *pos1;
	int cEnd = *pos2;

	int abeg = ovl->path.abpos;
	int aend = ovl->path.aepos;

	int bbeg = ovl->path.bbpos;
	int bend = ovl->path.bepos;

	int twidth = actx->twidth;
#ifdef DEBUG_STEP1A
	printf("getCorrespondingPositionInARead %d x %d, c(%d, %d), a(%d, %d) b(%d, %d)", ovl->aread, ovl->bread, cBeg, cEnd, abeg, aend, bbeg, bend);
	if (ovl->flags & OVL_COMP)
	printf(" COMP");
	printf("\n");
#endif

#ifdef DEBUG_STEP1A
	printf("comp: bbeg: %d, bend: %d\n", bbeg, bend);
#endif
	// find corresponding begin
	if (cBeg <= bbeg)
	{
#ifdef DEBUG_STEP1A
		printf("cBeg <= bbeg, %d <= %d\n", cBeg, bbeg);
#endif
		*pos1 = abeg;
	}
	else
	{
		if (ovl->path.tlen)
		{
#ifdef DEBUG_STEP1A
			printf("ccBeg: %d\n", cBeg);
#endif
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = ovl->path.bbpos;

			int j = 0;
			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\n", apos, bpos);
#endif
				if (bpos + twidth / 2 > cBeg)
				{
					apos += cBeg - bpos;
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\n", apos, bpos);
#endif
					*pos1 = -apos;
					break;
				}
				j += 2;
			}
		}
		else
		{
			if ((bend - cBeg) < (cBeg - bbeg))
			{
				*pos1 = -(abeg + (bend - cBeg));
			}
			else
			{
				*pos1 = -(aend - (cBeg - bbeg));
			}
		}

	}
	// find corresponding end
	if (cEnd >= bend)
	{
		*pos2 = aend;
	}
	else
	{
		if (ovl->path.tlen)
		{
#ifdef DEBUG_STEP1A
			printf("ccBeg: %d\n", cEnd);
#endif
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = ovl->path.bbpos;

			int j = 0;
			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\n", apos, bpos);
#endif
				if (bpos + twidth / 2 > cEnd)
				{
					apos += cEnd - bpos;
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\n", apos, bpos);
#endif

					*pos2 = -apos;
					break;
				}
				j += 2;
			}
		}
		else
		{
			if ((bend - cEnd) < (cEnd - bbeg))
			{
				*pos2 = -(abeg + (bend - cEnd));
			}
			else
			{
				*pos2 = -(aend - (cEnd - bbeg));
			}
		}
	}

#ifdef DEBUG_STEP1A
	printf("final range: %d, %d\n", *pos1, *pos2);
#endif
}

void getCorrespondingPositionInBRead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2)
{
	int cBeg = *pos1;
	int cEnd = *pos2;

	int abeg = ovl->path.abpos;
	int aend = ovl->path.aepos;

	int bbeg = ovl->path.bbpos;
	int bend = ovl->path.bepos;

	int twidth = actx->twidth;
	int transCBeg;
	int transCEnd;

#ifdef DEBUG_STEP1A
	{
		printf("getCorrespondingPositionInBRead %d x %d, c(%d, %d), a(%d, %d) b(%d, %d)", ovl->aread, ovl->bread, cBeg, cEnd, abeg, aend, bbeg, bend);
		if (ovl->flags & OVL_COMP)
		printf(" COMP");
		printf("\n");
	}
#endif
// find begin
	if (cBeg <= abeg)
	{
		transCBeg = bbeg;
#ifdef DEBUG_STEP1A
		printf("begin: %d <= %d, set transCBeg to %d\n", cBeg, abeg, transCBeg);
#endif
	}
	else
	{
		if (ovl->path.tlen)
		{
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = bbeg;

			int bb = 0;

			int j = 0;

			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];

				bb = bpos;
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif
				if (apos + twidth / 2 > cBeg)
				{
					bb += (cBeg - apos);
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif
					transCBeg = -bb;
					break;
				}
				j += 2;
			}
		}
		else
		{
#ifdef DEBUG_STEP1A
			printf("path.tlen is 0. set transCBeg\n");
#endif
			if ((cBeg - abeg) < (aend - cBeg))
			{
				transCBeg = -(bbeg + (cBeg - abeg));
			}
			else
			{
				transCBeg = -(bend - (aend - cBeg));
			}
		}

	}

	// find end
	if (cEnd >= aend)
	{
#ifdef DEBUG_STEP1A
		printf("end: %d >0 %d\n", cEnd, aend);
#endif
		transCEnd = bend;
	}
	else
	{
		if (ovl->path.tlen)
		{
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = bbeg;

			int bb = 0;

			int j = 0;

			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];

				bb = bpos;
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif

				if (apos + twidth / 2 > cEnd)
				{
					bb += (cEnd - apos);
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif
					transCEnd = -bb;
					break;
				}
				j += 2;
			}
		}
		else
		{
#ifdef DEBUG_STEP1A
			printf("transCEnd path.tlen is 0\n");
#endif

			if ((aend - cEnd) < (cEnd - abeg))
			{
				transCEnd = -(bend - (aend - cEnd));
			}
			else
			{
				transCEnd = -(bbeg + (cEnd - abeg));
			}

		}
	}

	if (ovl->flags & OVL_COMP)
	{
#ifdef DEBUG_STEP1A
		printf("switch positions transCBeg: %d, transCEnd: %d\n", transCBeg, transCEnd);
#endif
		*pos1 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - abs(transCBeg);
		*pos2 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - abs(transCEnd);

		if (transCBeg < 0)
			(*pos1) = -(*pos1);

		if (transCEnd < 0)
			(*pos2) = -(*pos2);
	}
	else
	{
		*pos1 = transCBeg;
		*pos2 = transCEnd;
	}
#ifdef DEBUG_STEP1A
	printf("getCorrespondingPositionInBRead: final range: %d, %d\n", *pos1, *pos2);
#endif
}

int processReadOverlapsAndMapThemToContigs(void* _ctx, Overlap* ovls, int novl)
{
	AnalyzeContext* actx = (AnalyzeContext*) _ctx;

	int i;

	int aread = ovls->aread;
	int bread;

	for (i = 0; i < novl; i++)
	{
		Overlap *ovl = ovls + i;
		bread = ovl->bread;

		if (ovl->flags & OVL_DISCARD)
		{
#ifdef DEBUG_STEP1A
			printf("skip overlap %d x %d, DISCARD\n", aread, bread);
#endif
			continue;
		}

		if (aread == bread)
		{
#ifdef DEBUG_STEP1A
			printf("skip identity overlap %d == %d\n", aread, bread);
#endif
			continue;
		}

		/// TODO ugly hack, to ignore reads that overlap multiple times with the current a-read, actually they should be removed in LAfilter, why are they still there????
		{
			int j = i;
			while (j + 1 < novl && (ovl->bread == ovls[j + 1].bread))
			{
				printf("ovl[%d] and ovl[%d] have same bread %d? 1(%d vs %d [%d, %d][%d,%d]) and 2(%d vs %d [%d, %d][%d,%d])\n", i, j + 1, ovl->bread, ovl->aread,
						ovl->bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, ovls[j + 1].aread, ovls[j + 1].bread, ovls[j + 1].path.abpos,
						ovls[j + 1].path.aepos, ovls[j + 1].path.bbpos, ovls[j + 1].path.bepos);
				j++;
			}

			if (j > i)
			{
#ifdef DEBUG_STEP1A
				printf("FUCK skip read %d. maps multiple times (%d) with current read %d\n", aread, j-i, bread);
#endif
				i = j;
				continue;
			}
		}

		int oLen = MAX(ovl->path.aepos - ovl->path.abpos, ovl->path.bepos - ovl->path.bbpos);

		if (oLen < actx->min_olen)
		{
#ifdef DEBUG_STEP1A
			printf("skip overlap %d x %d, len: %d\n", aread, bread, oLen);
#endif
			continue;
		}

		if (DB_READ_LEN(actx->patchedReadDB, ovl->bread) < actx->min_rlen)
		{
#ifdef DEBUG_STEP1A
			printf("skip overlap %d x %d, rlen: %d\n", aread, bread, DB_READ_LEN(actx->patchedReadDB, ovl->bread));
#endif
			continue;
		}

		// check if a-read is present in any contig, if yes than add corresponding breads to the contigs
		int j, k, l;
		for (j = 0; actx->vreadMask[j][aread] < 0; j++)
		{
			int conId = abs(actx->vreadMask[j][aread]);

			Contig * contig = actx->contigs + (conId - 1);

#ifdef DEBUG_STEP1A
			printf("found aread %d in contig %d: (len: %d, rep: %d) bread %d\n", aread, conId - 1, contig->property.len,
					getRepeatBasesFromInterval(actx, contig->property.contigID, 0, contig->property.len), bread);
#endif
			// check if corresponding b-read is already visited
			int found = 0;
			for (k = 0; k < actx->maxVReadMask; k++)
			{
				if (abs(actx->vreadMask[k][bread]) == conId)
					found = 1;
				if (actx->vreadMask[k][bread] == 0)
					break;
			}

			if (k == actx->maxVReadMask)
			{
				//printf("TERROR b-read %d is part of %d contigs!!!!\n", bread, k);
				l = 0;
				while (l < k)
				{
					//	printf("read %d is part of contig %d\n", bread, actx->vreadMask[l][bread]);
					l++;
				}
				//fflush(stdout);
				actx->maxVReadMask++;
				actx->vreadMask = (int**) realloc(actx->vreadMask, sizeof(int*) * actx->maxVReadMask);
				actx->vreadMask[k] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
				bzero(actx->vreadMask[k], sizeof(int) * DB_NREADS(actx->patchedReadDB));
			}

			// which coordinates of the a-read are used within the contig
			ContigRead *cread = getContigRead(contig, aread);

			assert(cread != NULL);

			int areadContigBeg = cread->patchedContigPosBeg;
			int areadContigEnd = cread->patchedContigPosEnd;

			int areadReadBeg = cread->ovlReads->beg;
			int areadReadEnd = cread->ovlReads->end;

#ifdef DEBUG_STEP1A
			printf("OVERLAP: %d vs %d [%8d, %8d] %c [%8d, %8d] l [%8d, %8d]\n",ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos, (ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos,
					DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread));
			printf("Aread %d in contig %d: cPos [%d, %d] rPos: [%d, %d] %c\n", aread, contig->property.contigID, areadContigBeg, areadContigEnd, areadReadBeg, areadReadEnd, (areadReadBeg > areadReadEnd) ? 'C' : 'N');
			printf("convert overlap coordinates from %d %d [%d, %d] %c [%d, %d]\n", ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos,
					(ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos);
#endif
			// complement sequence
			int cReadInComplement = 0;
			if (areadReadBeg > areadReadEnd)
			{
				int tmp = areadReadBeg;
				areadReadBeg = areadReadEnd;
				areadReadEnd = tmp;
				cReadInComplement = 1;
			}

			// do the a-read contig coordinates have an intersection with the b-read?
			int intersection = intersect(areadReadBeg, areadReadEnd, ovl->path.abpos, ovl->path.aepos);
			if (intersection == 0)
			{
#ifdef DEBUG_STEP1A
				printf("skip overlap %d x %d, No intersection\n", aread, bread);
#endif
				continue;
			}

			// mark b-read for current contig as visited
			if (!found)
				actx->vreadMask[k][bread] = conId;

			// check if OvlRead is already included due to transitivity
			int check = 1;
			{
				int u;

				for (u = 0; u < cread->numOvlReads; u++)
				{
					OvlRead *ovlU = cread->ovlReads + u;
					if (ovlU->patchedID == bread)
					{
						check = 0;
						break;
					}
				}
			}

			if (check)
			{
				// add valid bread to set of OvlReads (--> proper read for correction)
				if (cread->numOvlReads == cread->maxOvlReads)
				{
					cread->maxOvlReads = (int) (cread->maxOvlReads * 1.2) + 10;
					cread->ovlReads = (OvlRead*) realloc(cread->ovlReads, cread->maxOvlReads * sizeof(OvlRead));

					if (cread->ovlReads == NULL)
					{
						fprintf(stderr, "[ERROR]: Unable to increase OvlRead buffer of contig %d at read: %d, from %lu to %lu!\n", conId, aread,
								cread->numOvlReads * sizeof(OvlRead), cread->maxOvlReads * sizeof(OvlRead));
						exit(1);
					}
				}

				OvlRead *curOvlRead = cread->ovlReads + cread->numOvlReads;
				curOvlRead->patchedID = bread;

				int pos1, pos2;
				int abeg, aend;

				abeg = ovl->path.abpos;
				aend = ovl->path.aepos;

				pos1 = ovl->path.bbpos;
				pos2 = ovl->path.bepos;

				if (abeg < areadReadBeg)
					abeg = areadReadBeg;

				if (aend > areadReadEnd)
					aend = areadReadEnd;

				// contig read sticks in reverse complements in contig
				if (cReadInComplement)
				{
					if (aend > areadReadEnd)
						curOvlRead->cBeg = areadContigBeg;
					else
						curOvlRead->cBeg = areadContigBeg + (areadReadEnd - aend);

					if (abeg < areadReadBeg)
						curOvlRead->cEnd = areadContigEnd;
					else
						curOvlRead->cEnd = areadContigBeg + (areadReadEnd - abeg);
				}
				else
				{
					if (abeg < areadReadBeg)
						curOvlRead->cBeg = areadContigBeg;
					else
						curOvlRead->cBeg = areadContigBeg + (abeg - areadReadBeg);

					if (aend > areadReadEnd)
						curOvlRead->cEnd = areadContigEnd;
					else
						curOvlRead->cEnd = areadContigEnd - (areadReadEnd - aend);

				}

				if (abeg != ovl->path.abpos || aend != ovl->path.aepos)
				{
					pos1 = abeg;
					pos2 = aend;
#ifdef DEBUG_STEP1A
					printf("getCorrespondingPositionInBRead(actx, ovl, &pos1, &pos2);\n");
					getCorrespondingPositionInBRead(actx, ovl, &pos1, &pos2);
#endif
				}

#ifdef DEBUG_STEP1A
				printf("set pos1 to %d and pos2 to %d\n", pos1, pos2);
#endif

				if (cReadInComplement)
				{
					curOvlRead->beg = pos2;
					curOvlRead->end = pos1;
				}
				else
				{
					curOvlRead->beg = pos1;
					curOvlRead->end = pos2;
				}
				cread->numOvlReads++;

				assert(curOvlRead->cBeg < curOvlRead->cEnd);
				assert(curOvlRead->cBeg >= cread->patchedContigPosBeg);
				assert(curOvlRead->cEnd <= cread->patchedContigPosEnd);

#ifdef DEBUG_STEP1A
				ContigRead *falk = getContigRead(contig, aread);
				assert(falk != NULL);
				printf("add overlap %d x %d, a(%d, %d), b(%d, %d), l(%d, %d), %c at ConPos [%d, %d] in Read [%d, %d]\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos,
						ovl->path.bepos, DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread), (ovl->flags & OVL_COMP) ? 'c' : 'n',
						falk->patchedContigPosBeg, falk->patchedContigPosEnd, falk->ovlReads->beg, falk->ovlReads->end);

				printf("=======<<<< added OVL-B-read: %d %d-%d, in c%d-%d-%d %c %c crLen %d oLen %d ocrLen %d ocLen %d\n", curOvlRead->patchedID, curOvlRead->beg, curOvlRead->end, contig->property.contigID, curOvlRead->cBeg, curOvlRead->cEnd,
						(cReadInComplement) ? 'C': 'N', (ovl->flags & OVL_COMP) ? 'C': 'N', cread->ovlReads->cEnd - cread->ovlReads->cBeg, ovl->path.aepos - ovl->path.abpos, curOvlRead->cEnd - curOvlRead->cBeg, abs(abs(curOvlRead->beg)-abs(curOvlRead->end)));

#endif
			}

		}

		/// As there is no proper transitivity given, we have to check if one of the b-reads from all overlaps of the current
		/// aread is part of a contig as well, and was not already incorporated as valid overlap
		// check if b-read is present in any contig, if yes than add corresponding breads to the contigs
		for (j = 0; actx->vreadMask[j][bread] < 0; j++)
		{

			int conId = abs(actx->vreadMask[j][bread]);

			Contig * contig = actx->contigs + (conId - 1);

#ifdef DEBUG_STEP1A
			printf("found bread %d in contig %d: (len: %d, rep: %d)\n", bread, conId - 1, contig->property.len,
					getRepeatBasesFromInterval(actx, contig->property.contigID, 0, contig->property.len));
#endif
			// check if corresponding a-read is already visited
			int found = 0;
			for (k = 0; k < actx->maxVReadMask; k++)
			{
				if (abs(actx->vreadMask[k][aread]) == conId)
					found = 1;
				if (actx->vreadMask[k][aread] == 0)
					break;
			}

			if (k == actx->maxVReadMask)
			{
				//printf("TERROR a-read %d is part of %d contigs!!!!\n", aread, k);
				l = 0;
				while (l < k)
				{
					//printf("read %d is part of contig %d\n", aread, actx->vreadMask[l][aread]);
					l++;
				}
				//fflush(stdout);
				actx->maxVReadMask++;
				actx->vreadMask = (int**) realloc(actx->vreadMask, sizeof(int*) * actx->maxVReadMask);
				actx->vreadMask[k] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
				bzero(actx->vreadMask[k], sizeof(int) * DB_NREADS(actx->patchedReadDB));
			}

			ContigRead *cread = getContigRead(contig, bread);
			assert(cread != NULL);

			int breadContigBeg = cread->patchedContigPosBeg;
			int breadContigEnd = cread->patchedContigPosEnd;

			int breadReadBeg = cread->ovlReads->beg;
			int breadReadEnd = cread->ovlReads->end;

#ifdef DEBUG_STEP1A
			printf("OVERLAP: %d vs %d [%8d, %8d] %c [%8d, %8d] l [%8d, %8d]\n",ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos, (ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos,
					DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread));
			printf("Bread %d in contig %d: cPos [%d, %d] rPos: [%d, %d]\n", bread, contig->property.contigID, breadContigBeg, breadContigEnd, breadReadBeg, breadReadEnd);
			printf("convert overlap coordinates from %d %d [%d, %d] %c [%d, %d]\n", ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos,
					(ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos);
#endif
			int bBeg, bEnd;

			int cReadInComplement = 0;

			// complement sequence
			if (breadReadBeg > breadReadEnd)
			{
				int tmp = breadReadBeg;
				breadReadBeg = breadReadEnd;
				breadReadEnd = tmp;
				cReadInComplement = 1;
			}
#ifdef DEBUG_STEP1A
			printf("cBeg: %d, cEnd: %d\n", breadContigBeg, breadContigEnd);
#endif
			if (ovl->flags & OVL_COMP)
			{
				bBeg = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - ovl->path.bepos;
				bEnd = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - ovl->path.bbpos;
			}
			else
			{
				bBeg = ovl->path.bbpos;
				bEnd = ovl->path.bepos;
			}
#ifdef DEBUG_STEP1A
			printf("get intersection from %d in contig (%d, %d) and b(%d, %d)\n", bread, breadReadBeg, breadReadEnd, bBeg, bEnd);
#endif
			int intersection = intersect(breadReadBeg, breadReadEnd, bBeg, bEnd);
			if (intersection == 0)
			{
				continue;
			}

			// mark a-read for current contig as visited
			if (!found)
				actx->vreadMask[k][aread] = conId;

			// check if we added the same OvlRead twice due to symmetry
			int check = 1;
			{
				int u;

				for (u = 0; u < cread->numOvlReads; u++)
				{
					OvlRead *ovlU = cread->ovlReads + u;
					if (ovlU->patchedID == aread)
					{
						check = 0;
						break;
					}

				}
			}

			if (check)
			{
				if (cread->numOvlReads == cread->maxOvlReads)
				{
					cread->maxOvlReads = (int) (cread->maxOvlReads * 1.2) + 10;
					cread->ovlReads = (OvlRead*) realloc(cread->ovlReads, sizeof(OvlRead) * cread->maxOvlReads);

					if (cread->ovlReads == NULL)
					{
						fprintf(stderr, "[ERROR]: Unable to increase OvlRead buffer of contig %d at read: %d, from %d to %d!\n", contig->property.contigID, bread,
								cread->numOvlReads, cread->maxOvlReads);
						exit(1);
					}
				}

				// add valid a-read to set of OvlReads (--> proper read for correction)
				OvlRead *curOvlRead = cread->ovlReads + cread->numOvlReads;
				curOvlRead->patchedID = aread;

				int pos1, pos2;

				if (ovl->flags & OVL_COMP)
				{
					if (bEnd > breadReadEnd)
						pos1 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bEnd + (bEnd - breadReadEnd);
					else
						pos1 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bEnd;

					if (bBeg < breadReadBeg)
						pos2 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bBeg - (breadReadBeg - bBeg);
					else
						pos2 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bBeg;

				}
				else
				{
					if (bBeg < breadReadBeg)
						pos1 = breadReadBeg;
					else
						pos1 = bBeg;

					if (bEnd > breadReadEnd)
						pos2 = breadReadEnd;
					else
						pos2 = bEnd;
				}
#ifdef DEBUG_STEP1A
				printf("set pos1 to %d and pos2 to %d\n", pos1, pos2);
#endif
				if (cReadInComplement)
				{
					if (bEnd > breadReadEnd)
						curOvlRead->cBeg = breadContigBeg;
					else
						curOvlRead->cBeg = breadContigBeg + (breadReadEnd - bEnd);

					if (bBeg < breadReadBeg)
						curOvlRead->cEnd = breadContigEnd;
					else
						curOvlRead->cEnd = breadContigBeg + (breadReadEnd - bBeg);
				}
				else
				{
					if (bBeg < breadReadBeg)
						curOvlRead->cBeg = breadContigBeg;
					else
						curOvlRead->cBeg = breadContigBeg + (bBeg - breadReadBeg);

					if (bEnd > breadReadEnd)
						curOvlRead->cEnd = breadContigEnd;
					else
						curOvlRead->cEnd = breadContigEnd - (breadReadEnd - bEnd);
				}

				getCorrespondingPositionInARead(actx, ovl, &pos1, &pos2);

				if (cReadInComplement != (ovl->flags & OVL_COMP))
				{
					curOvlRead->beg = pos2;
					curOvlRead->end = pos1;
				}
				else
				{
					curOvlRead->beg = pos1;
					curOvlRead->end = pos2;
				}

				cread->numOvlReads++;

#ifdef DEBUG_STEP1A
				ContigRead *falk = getContigRead(contig, bread);
				assert(falk != NULL);

				printf("add overlap %d x %d, a(%d, %d), b(%d, %d), l(%d, %d) %c at ConPos [%d, %d] in Read [%d, %d]\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos,
						ovl->path.bepos, DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread), (ovl->flags & OVL_COMP) ? 'c' : 'n',
						falk->patchedContigPosBeg, falk->patchedContigPosEnd, falk->ovlReads->beg, falk->ovlReads->end);

				printf("=======>>>> added OVL-A-read: %d %d-%d, in c%d-%d-%d %c %c crLen %d oLen %d ocrLen %d ocLen %d\n", curOvlRead->patchedID, curOvlRead->beg, curOvlRead->end, contig->property.contigID, curOvlRead->cBeg, curOvlRead->cEnd,
						(cReadInComplement) ? 'C': 'N', (ovl->flags & OVL_COMP) ? 'C': 'N', cread->ovlReads->cEnd - cread->ovlReads->cBeg, ovl->path.aepos - ovl->path.abpos, curOvlRead->cEnd - curOvlRead->cBeg, abs(abs(curOvlRead->beg)-abs(curOvlRead->end)));
#endif
			}
		}
	}
	return 1;
}

int cmpOVLreadById(const void *a, const void *b)
{
	return ((const OvlRead*) a)->patchedID - ((const OvlRead*) b)->patchedID;
}

int getRepeatBasesFromInterval(AnalyzeContext* actx, int contigDB, int readID, int beg, int end)
{
	track_anno* rep_anno;
	track_data* rep_data;

	if (contigDB)
	{
		rep_anno = actx->corContigRepeats_track->anno;
		rep_data = actx->corContigRepeats_track->data;

		if (readID < 0 || readID >= DB_NREADS(actx->corContigDB))
		{
			fprintf(stderr, "[ERROR] - getRepeatBasesFromInterval readID: %d out of bounds [0, %d]\n", readID, DB_NREADS(actx->corContigDB) - 1);
			fflush(stderr);
			exit(1);
		}
	}
	else // patched reads db
	{
		rep_anno = actx->patchedReadRepeat_track->anno;
		rep_data = actx->patchedReadRepeat_track->data;

		if (readID < 0 || readID >= DB_NREADS(actx->patchedReadDB))
		{
			fprintf(stderr, "[ERROR] - getRepeatBasesFromInterval readID: %d out of bounds [0, %d]\n", readID, DB_NREADS(actx->patchedReadDB) - 1);
			fflush(stderr);
			exit(1);
		}
	}

	track_anno rb, re;

	int repBases = 0;
	int rBeg, rEnd;

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

void analyzeContigCoverageOfMappedReads(AnalyzeContext *actx)
{
	int i, j, k;

	// sort b-reads according read id and remove duplicates
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;

		if (contig->property.cflag & CLASS_CONTIG_DISCARD)
			continue;

		for (j = 0; j < contig->numcReads; j++)
		{
			qsort(contig->cReads[j].ovlReads, contig->cReads[j].numOvlReads, sizeof(OvlRead), cmpOVLreadById);
		}
	}

	// TODO pass variable, set by user ?
	int COV_BIN_SIZE = 500;
	int MIN_COV = 3;

	bit *inReads = ba_new(DB_NREADS(actx->patchedReadDB));
	bit *outReads = ba_new(DB_NREADS(actx->patchedReadDB));
	int *readCovHist = (int *) malloc(sizeof(int) * (DB_READ_MAXLEN(actx->patchedReadDB) / COV_BIN_SIZE + 1));

	// update cBeg, End intervals for each B-read according global contig positions
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;

		if (contig->property.cflag & CLASS_CONTIG_DISCARD)
			continue;

#ifdef DEBUG_STEP1B

		printf("Contig %d, reads: %d\n", contig->property.contigID, contig->numcReads);
#endif
		for (j = 0; j < contig->numcReads; j++)
		{
			ContigRead *cread = contig->cReads + j;

#ifdef DEBUG_STEP1B

			printf("prev: %d cbe %d %d, be %d, %d", cread->patchedID, cread->patchedContigPosBeg, cread->patchedContigPosEnd, cread->ovlReads->beg,cread->ovlReads->end);
#endif

#ifdef DEBUG_STEP1B
			printf(" comp? %d, crlen: %d, nOvlReads %d\n", crComp, DB_READ_LEN(actx->patchedReadDB, abs(cread->patchedID)), cread->numOvlReads);
#endif
			for (k = 1; k < cread->numOvlReads; k++)
			{
#ifdef DEBUG_STEP1B
				printf("\t %d beg,end [%d, %d], cBeg,cEnd [%d, %d]\n", cread->ovlReads[k].patchedID, cread->ovlReads[k].beg, cread->ovlReads[k].end,
						cread->ovlReads[k].cBeg, cread->ovlReads[k].cEnd);
#endif
			}

			{
				int l, m;
				// 1. analyze in-coverage and mark them as visited
				ba_all_assign(inReads, DB_NREADS(actx->patchedReadDB), FALSE);

				for (l = 1; l < cread->numOvlReads; l++)
				{
					if (cread->patchedContigPosBeg == cread->ovlReads[l].cBeg)
					{
						ba_assign(inReads, cread->ovlReads[l].patchedID, TRUE);
						cread->nInReads++;
					}
				}

				// 2. analyze intersection of out-reads from prev reads and in-reads of current read
				if (j) // skip first contig read
				{
					/////
					/////                   r_i    r_j     r_k     r_l     r_m
					/////    contigA:    -------|-------|-------|-------|-------
					/////                ---    |       |       |       |                fixed reads that overlap with each r_i, ... , r_m
					/////                 ---   |       |       |       |
					///                     --- |       |       |       |
					/////                      1 2     1 2     1 2     1 2               1 ... OUT reads (arrival rate) and 2 ... IN reads (start rate)

					int intersect = 0;
					/// loop over IN-reads of current contig read
					for (l = 1; l < cread->numOvlReads; l++)
					{
						if (ba_value(outReads, cread->ovlReads[l].patchedID))
							intersect++;
					}
					contig->cReads[j - 1].nComReadsWithNextContigRead = intersect;
				}

				// 3. analyze
				//  a) out-coverage and mark them as visited
				//  b) average coverage
				//  c) low coverage
				ba_all_assign(outReads, DB_NREADS(actx->patchedReadDB), FALSE);
				bzero(readCovHist, sizeof(int) * (DB_READ_MAXLEN(actx->patchedReadDB) / COV_BIN_SIZE + 1));
				// todo add contig read as well --> i.e. in most cases the number of common exiting/entering reads is not zero
//				ba_assign(outReads, cread->id, TRUE);
				//3a) + 3b)
				for (l = 1; l < cread->numOvlReads; l++)
				{
					if (cread->patchedContigPosEnd == cread->ovlReads[l].cEnd)
					{
						ba_assign(outReads, cread->ovlReads[l].patchedID, TRUE);
						cread->nOutReads++;
					}
					cread->avgCov += cread->ovlReads[l].cEnd - cread->ovlReads[l].cBeg;
				}
				contig->avgCov += cread->avgCov;
				cread->avgCov /= (cread->patchedContigPosEnd - cread->patchedContigPosBeg);

				// 3c) i.e. create coverage histogram  ...
				//printf("Coverage Histogram for read: %d PreOut-CurIn-Cov %d,%d\n", contig->reads[k][0].id, curCheckPoint[curNum-1].nComReadsWithNextRead, curCheckPoint[curNum].nComReadsWithprevRead);
				m = cread->patchedContigPosBeg;
				while (m < cread->patchedContigPosEnd)
				{
					int cFrom, cTo;
					cFrom = m;
					if (m + COV_BIN_SIZE > cread->patchedContigPosEnd)
						cTo = cread->patchedContigPosEnd;
					else
						cTo = cFrom + COV_BIN_SIZE;

					int bin = (cFrom - cread->patchedContigPosBeg) / COV_BIN_SIZE;

					// increase coverage count only if read fully spans histogram bin
					for (l = 1; l < cread->numOvlReads; l++)
					{
						if (cread->ovlReads[l].cBeg <= cFrom && cread->ovlReads[l].cEnd >= cTo)
							readCovHist[bin]++;
					}

					//				printf("<%9d, %9d> c %3d\n", cFrom, cTo, readCovHist[bin]);
					m = cTo;
				}

				// 3c) ... and check coverage drops
				for (l = 0, m = cread->patchedContigPosBeg; m < cread->patchedContigPosEnd; l++, m += COV_BIN_SIZE)
				{
					if (readCovHist[l] < MIN_COV)
					{
						cread->lowCov = 1;
						break;
					}
				}
			}
		}
		// calculate average coverage of contig
		contig->avgCov /= contig->property.len;

//		if(contig->property.contigID == 0)
		{
			printf("STATS of contig %d\n", contig->property.contigID);
			printf("avgCov: %.2f\n", contig->avgCov);
			for (j = 0; j < contig->numcReads; j++)
			{
				ContigRead *cread = contig->cReads + j;
				printf("%9d - %9d, id %10d, c[%9d, %9d] r[%6d,%6d] inR %10d, outR %d, allR %d, avgCov %.2f, comWNR %d\n", cread->patchedContigPosBeg,
						cread->patchedContigPosEnd, cread->patchedID, cread->ovlReads->cBeg, cread->ovlReads->cEnd, cread->ovlReads->beg, cread->ovlReads->end,
						cread->nInReads, cread->nOutReads, cread->numOvlReads, cread->avgCov, cread->nComReadsWithNextContigRead);
			}
		}
	}

	// cleanup read coverage buffers
	free(readCovHist);
	free(inReads);
	free(outReads);
}

int cmpContigLength(const void *a, const void *b)
{
	return ((const ContigIDAndLen*) b)->len - ((const ContigIDAndLen *) a)->len;
}

int getPathID(AnalyzeContext *actx, int contigID)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getPathID: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	if (Check_Track(actx->corContigDB, "path") != 0)
		return -1;

	HITS_TRACK * pathTrack = track_load(actx->corContigDB, "path");

	assert(pathTrack != NULL);

	track_anno* path_anno = pathTrack->anno;
	track_data* path_data = pathTrack->data;

	track_anno b, e;

	b = path_anno[contigID] / sizeof(track_data);
	e = path_anno[contigID + 1] / sizeof(track_data);

	assert(b < e);

	return (int) path_data[b];
}

int getFileID(FileNamesAndOffsets *fileAndNameOffset, int contigID)
{
	assert(fileAndNameOffset != NULL);
	assert(contigID >= 0 && contigID <= fileAndNameOffset->toDbIdx[fileAndNameOffset->numFileNames - 1]);

	int i;
	for (i = 0; i < fileAndNameOffset->numFileNames; i++)
	{
		if (contigID >= fileAndNameOffset->fromDbIdx[i] && contigID <= fileAndNameOffset->toDbIdx[i])
		{
			break;
		}
	}
	return i;
}

int getNumberOfSubgraphContigs(FileNamesAndOffsets *fileAndNameOffset, int contigID)
{
	assert(fileAndNameOffset != NULL);
	assert(contigID >= 0 && contigID <= fileAndNameOffset->toDbIdx[fileAndNameOffset->numFileNames - 1]);

	int i;
	for (i = 0; i < fileAndNameOffset->numFileNames; i++)
	{
		if (contigID >= fileAndNameOffset->fromDbIdx[i] && contigID <= fileAndNameOffset->toDbIdx[i])
		{
			break;
		}
	}

	return (fileAndNameOffset->toDbIdx[i] - fileAndNameOffset->fromDbIdx[i]) + 1;
}

void getContigsEndIDs(AnalyzeContext *actx, int contigID, int* beg, int *end)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getContigsEndIDs: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	if (Check_Track(actx->corContigDB, "ends") != 0)
		return;

	HITS_TRACK * pathTrack = track_load(actx->corContigDB, "ends");

	assert(pathTrack != NULL);

	track_anno* path_anno = pathTrack->anno;
	track_data* path_data = pathTrack->data;

	track_anno b, e;

	b = path_anno[contigID] / sizeof(track_data);
	e = path_anno[contigID + 1] / sizeof(track_data);

	assert(b < e);

	(*beg) = path_data[b];
	(*end) = path_data[b + 1];
}

void preclassifyContigsByBReadsAndPath(AnalyzeContext *actx)
{
	int j, k, l, m;

	ContigIDAndLen *contigIDAndLength = (ContigIDAndLen*) malloc(sizeof(ContigIDAndLen) * actx->numContigs);
	int numCL = 0;

	// todo define a proper struct?
	bit * contigJReads, *contigKReads, *contigIntersectionJKReads;
	int * contigJBegRange, *contigJEndRange;
	int * contigKBegRange, *contigKEndRange;

	int * contigJCovHist, *contigKCovHist;
	int binSize = 100;
	int nBins = (DB_READ_MAXLEN(actx->corContigDB) / binSize) + 1;
	contigJCovHist = (int*) malloc(sizeof(int) * nBins);
	contigKCovHist = (int*) malloc(sizeof(int) * nBins);

	int numReads = DB_NREADS(actx->patchedReadDB);
	contigJReads = ba_new(numReads); // all raw reads of contig J
	contigKReads = ba_new(numReads); // all raw reads of contig K
	contigIntersectionJKReads = ba_new(numReads); // all reads that are common in both contigs K and J

	contigJBegRange = (int*) malloc(sizeof(int) * numReads);
	contigJEndRange = (int*) malloc(sizeof(int) * numReads);
	contigKBegRange = (int*) malloc(sizeof(int) * numReads);
	contigKEndRange = (int*) malloc(sizeof(int) * numReads);

	int faId;

	//for (faId = 0; faId < actx->contigFileNamesAndOffsets->numFileNames; faId++)
	{
//#ifdef DEBUG_STEP1C
//		printf("i %d, %d %d\n", faId, actx->contigFileNamesAndOffsets->fromDbIdx[faId], actx->contigFileNamesAndOffsets->toDbIdx[faId]);
//#endif

//		if (actx->contigFileNamesAndOffsets->fromDbIdx[faId] == actx->contigFileNamesAndOffsets->toDbIdx[faId])
//		{
//#ifdef DEBUG_STEP1C
//			printf("fasta file %d: %s --> has only one contig\n", faId, actx->contigFileNamesAndOffsets->fileNames[faId]);
//#endif
//			actx->contigs[actx->contigFileNamesAndOffsets->fromDbIdx[faId]].property.readRelationFlags |= READ_IS_UNIQUE;
//		}
//		else // at this point, current compound contains at least 2 contigs
		{
			numCL = 0;

			// add all contigs from current connected compound (same file == subgraph)
//			for (j = actx->contigFileNamesAndOffsets->fromDbIdx[faId]; j <= actx->contigFileNamesAndOffsets->toDbIdx[faId]; j++)
			for (j = 0; j < actx->numContigs; j++)
			{
				Contig *contig = actx->contigs + j;

				if (contig->property.cflag & CLASS_CONTIG_DISCARD)
					continue;
				contigIDAndLength[numCL].idx = contig->property.contigID;
				contigIDAndLength[numCL].len = contig->property.len;
#ifdef DEBUG_STEP1C
				printf("%d len: %d\n", contig->property.contigID, contig->property.len);
#endif
				numCL++;
			}

			// sort contigs according length
			qsort(contigIDAndLength, numCL, sizeof(contigIDAndLength), cmpContigLength);

			// check
			for (j = 0; j < numCL; j++)
			{
				Contig * contig_j = actx->contigs + contigIDAndLength[j].idx;
				int jBins = (contig_j->property.len / binSize) + 1;

				if (contig_j->property.cflag & CLASS_CONTIG_DISCARD)
					continue;

#ifdef DEBUG_STEP1C
				printf("Check against contig_j: %d, len %d\n", contig_j->property.contigID, contig_j->property.len);
#endif

//				if (contig_j->gClassificFlag & CONTIG_IS_CONTAINED)
//				{
//#ifdef DEBUG_STEP1C
//					printf("Skip contig %d, already contained!\n", contig_j->property.contigID);
//#endif
//					continue;
//				}

				ba_all_assign(contigJReads, numReads, FALSE);
				memset(contigJBegRange, -1, sizeof(int) * numReads);
				memset(contigJEndRange, -1, sizeof(int) * numReads);

				// fill breads and positions for contig_j
				for (k = 0; k < contig_j->numcReads; k++)
				{
					ContigRead *cread = contig_j->cReads + k;

#ifdef DEBUG_STEP1C
					printf("contig_j %d cread %d [nBreads %d]\n", contig_j->property.contigID, cread->patchedID, cread->numOvlReads);
#endif
					for (l = 0; l < cread->numOvlReads; l++)
					{
						int readID = abs(cread->ovlReads[l].patchedID);

						ba_assign(contigJReads, readID, TRUE);
						if (contigJBegRange[readID] < 0)
						{
							contigJBegRange[readID] = cread->ovlReads[l].cBeg;
							contigJEndRange[readID] = cread->ovlReads[l].cEnd;
						}
						else
						{
							if (contigJBegRange[readID] > cread->ovlReads[l].cBeg)
								contigJBegRange[readID] = cread->ovlReads[l].cBeg;
							if (contigJEndRange[readID] < cread->ovlReads[l].cEnd)
								contigJEndRange[readID] = cread->ovlReads[l].cEnd;
						}
					}
				}

				// compare all other contigs from current subgraph against contig_j
				for (k = j + 1; k < numCL; k++)
				{
					Contig * contig_k = actx->contigs + contigIDAndLength[k].idx;
					int kBins = (contig_k->property.len / binSize) + 1;

					if (contig_k->property.cflag & CLASS_CONTIG_DISCARD)
						continue;

					if ((contig_k->property.rflag & REL_CONTIG_IS_ALT))
						continue;

					if (contig_k->property.len > 300000)
						continue;

//					if (contig_k->gClassificFlag & CONTIG_IS_CONTAINED)
//					{
//#ifdef DEBUG_STEP1C
//					printf("Skip contig %d, already contained!\n", contig_k->property.contigID);
//#endif
//						continue;
//					}

					ba_all_assign(contigKReads, numReads, FALSE);
					memset(contigKBegRange, -1, sizeof(int) * numReads);
					memset(contigKEndRange, -1, sizeof(int) * numReads);

					for (l = 0; l < contig_k->numcReads; l++)
					{
						ContigRead *cread = contig_k->cReads + l;

#ifdef DEBUG_STEP1C
						printf("contig_k %d add breads for cread %d\n", contig_k->property.contigID, cread->patchedID);
#endif
						for (m = 0; m < cread->numOvlReads; m++)
						{
							int readID = abs(cread->ovlReads[m].patchedID);

							ba_assign(contigKReads, readID, TRUE);
							if (contigKBegRange[readID] < 0)
							{
								contigKBegRange[readID] = cread->ovlReads[m].cBeg;
								contigKEndRange[readID] = cread->ovlReads[m].cEnd;
							}
							else
							{
								if (contigKBegRange[readID] > cread->ovlReads[m].cBeg)
									contigKBegRange[readID] = cread->ovlReads[m].cBeg;
								if (contigKEndRange[readID] < cread->ovlReads[m].cEnd)
									contigKEndRange[readID] = cread->ovlReads[m].cEnd;
							}
						}
					}

					// get read intersection between contig_k and contig_j
					ba_intersection(contigKReads, contigJReads, &contigIntersectionJKReads, numReads, numReads);

					int numComReads = ba_count(contigIntersectionJKReads, numReads);
					int numkReads = ba_count(contigKReads, numReads);

					if (numComReads * 100.0 / numkReads > 1.0)
					{
						memset(contigJCovHist, 0, sizeof(int) * (nBins));
						memset(contigKCovHist, 0, sizeof(int) * (nBins));

						int numCovIvlJ = 0;
						int *covIvlJ = (int*) malloc(sizeof(int) * ((numCovIvlJ + 1) * 3));
						int numCovIvlK = 0;
						int *covIvlK = (int*) malloc(sizeof(int) * ((numCovIvlK + 1) * 3));

						for (l = 0; l < numReads; l++)
						{
							if (ba_value(contigIntersectionJKReads, l))
							{
								int from = contigJBegRange[l] / binSize;
								int to = contigJEndRange[l] / binSize;
//								if (from < 0)
//								{
//									printf("contigIntersectionJKReads Jrange failed for read %d: from %d to %d\n", l, from, to);
//								}

//								if (to >= nBins)
//								{
//									printf("contigIntersectionJKReads Jrange failed for read %d: from %d to %d\n", l, from, to);
//								}

								while (from <= to)
								{
									contigJCovHist[from]++;
									from++;
								}
								from = contigKBegRange[l] / binSize;
								to = contigKEndRange[l] / binSize;

//								if (from < 0)
//								{
//									printf("contigIntersectionJKReads Krange failed for read %d: from %d to %d\n", l, from, to);
//								}
//
//								if (to >= nBins)
//								{
//									printf("contigIntersectionJKReads Krange failed for read %d: from %d to %d\n", l, from, to);
//								}
								while (from <= to)
								{
									contigKCovHist[from]++;
									from++;
								}
							}
						}
#ifdef DEBUG_STEP1C
						printf("analyze contigs: %d p%d (cov %.2f) %d p%d (cov %.2f) len %d %d --> numOfcommonReads: %d of %d: %.3f\n", contig_j->property.contigID,
								contig_j->property.pathID, contig_j->avgCov, contig_k->property.contigID, contig_k->property.pathID, contig_k->avgCov, contig_j->property.len,
								contig_k->property.len, numComReads, numkReads, numComReads * 100.0 / numkReads);
						printf("contigKCovHist:\n");
#endif
						int firstKbin, lastKbin, firstJbin, lastJbin;
						int cumHistKreads, cumHistJreads;
						cumHistKreads = cumHistJreads = 0;
						firstKbin = lastKbin = firstJbin = lastJbin = -1;

						for (l = 0; l < kBins; l++)
						{
							if (contigKCovHist[l] > 0)
							{
								if (firstKbin == -1)
									firstKbin = l;

								lastKbin = l;
								cumHistKreads += contigKCovHist[l];
//#ifdef DEBUG_STEP1C
//								printf("%10d-%10d: %4d\n", l * binSize, (l * binSize) + binSize - 1, contigKCovHist[l]);
//#endif
							}
							else if (firstKbin != -1)  // add coverage block
							{
								covIvlK[numCovIvlK * 3] = firstKbin * binSize;
								covIvlK[numCovIvlK * 3 + 1] = lastKbin * binSize + (binSize - 1);
								covIvlK[numCovIvlK * 3 + 2] = cumHistKreads / (lastKbin - firstKbin + 1);
								numCovIvlK++;
								covIvlK = (int*) realloc(covIvlK, sizeof(int) * ((numCovIvlK + 1) * 3));
								firstKbin = -1;
								cumHistKreads = 0;
							}
						}

						if (firstKbin != -1)
						{
							covIvlK[numCovIvlK * 3] = firstKbin * binSize;
							covIvlK[numCovIvlK * 3 + 1] = lastKbin * binSize + (binSize - 1);
							covIvlK[numCovIvlK * 3 + 2] = cumHistKreads / (lastKbin - firstKbin + 1);
							numCovIvlK++;
						}

//#ifdef DEBUG_STEP1C
//						printf("contigJCovHist:\n");
//#endif
						for (l = 0; l < jBins; l++)
						{
							if (contigJCovHist[l] > 0)
							{
								if (firstJbin == -1)
									firstJbin = l;
								cumHistJreads += contigJCovHist[l];
								lastJbin = l;
//#ifdef DEBUG_STEP1C
//								printf("%10d-%10d: %4d\n", l * binSize, (l * binSize) + binSize - 1, contigJCovHist[l]);
//#endif
							}
							else if (firstJbin != -1)
							{
								covIvlJ[numCovIvlJ * 3] = firstJbin * binSize;
								covIvlJ[numCovIvlJ * 3 + 1] = lastJbin * binSize + (binSize - 1);
								covIvlJ[numCovIvlJ * 3 + 2] = cumHistJreads / (lastJbin - firstJbin + 1);
								numCovIvlJ++;
								covIvlJ = (int*) realloc(covIvlJ, sizeof(int) * ((numCovIvlJ + 1) * 3));
								firstJbin = -1;
								cumHistJreads = 0;
							}
						}
						if (firstJbin != -1)
						{
							covIvlJ[numCovIvlJ * 3] = firstJbin * binSize;
							covIvlJ[numCovIvlJ * 3 + 1] = lastJbin * binSize + (binSize - 1);
							covIvlJ[numCovIvlJ * 3 + 2] = cumHistJreads / (lastJbin - firstJbin + 1);
							numCovIvlJ++;
						}

#ifdef DEBUG_STEP1C
						printf("common reads for both contigs\n");

						m = 0;
						for (l = 0; l < numReads; l++)
						{
							if (ba_value(contigIntersectionJKReads, l))
							{
								printf("%4d: c %7d vs c %7d read: %7d [%7d, %7d] --> [%7d, %7d]\n", m++, contig_j->property.contigID, contig_k->property.contigID, l,
										contigJBegRange[l], contigJEndRange[l], contigKBegRange[l], contigKEndRange[l]);
							}
						}
#endif

						int cumKGapBases = covIvlK[0] + contig_k->property.len - covIvlK[3 * numCovIvlK - 2];
						int cumKCoveredBases = 0;
						for (l = 0; l < numCovIvlK; l++)
						{
#ifdef DEBUG_STEP1C
							printf("CoveredBasesInK [%8d, %8d, cov %3d]\n", covIvlK[l * 3], covIvlK[l * 3 + 1], covIvlK[l * 3 + 2]);
#endif
							cumKCoveredBases += (covIvlK[l * 3 + 1] - covIvlK[l * 3]);

							if (l + 1 < numCovIvlK)
								cumKGapBases += (covIvlK[(l + 1) * 3] - covIvlK[l * 3 + 1]);
						}
#ifdef DEBUG_STEP1C
						printf("cumKGapBases %d (%.2f%%) cumKCoveredBases %d (%.2f%%)\n", cumKGapBases, cumKGapBases * 100.0 / contig_k->property.len, cumKCoveredBases,
								cumKCoveredBases * 100.0 / contig_k->property.len);
#endif
						int cumJGapBases = 0;
						int cumJCoveredBases = 0;
						for (l = 0; l < numCovIvlJ; l++)
						{
#ifdef DEBUG_STEP1C
							printf("CoveredBasesInJ [%8d, %8d, cov %3d]\n", covIvlJ[l * 3], covIvlJ[l * 3 + 1], covIvlJ[l * 3 + 2]);
#endif
							cumJCoveredBases += (covIvlJ[l * 3 + 1] - covIvlJ[l * 3]);

							if (l + 1 < numCovIvlJ)
								cumJGapBases += (covIvlJ[(l + 1) * 3] - covIvlJ[l * 3 + 1]);
						}
#ifdef DEBUG_STEP1C
						printf("cumJGapBases %d (%.2f%%) cumJCoveredBases %d (%.2f%%)\n", cumJGapBases, cumJGapBases * 100.0 / (covIvlJ[3 * numCovIvlJ - 2] - covIvlJ[0]),
								cumJCoveredBases, cumJCoveredBases * 100.0 / (covIvlJ[3 * numCovIvlJ - 2] - covIvlJ[0]));
#endif
						// preliminary classification
						// contig k is contained in J
						if (numComReads * 100.0 / numkReads > actx->contByReads_CommonReadFraction || (cumKCoveredBases * 100.0 / contig_k->property.len >= actx->contByReads_CoveredLenPerc))
						{
							int preClassFlagK = REL_READ_IS_ALT;
							int preClassFlagJ = REL_READ_HAS_ALT;

#ifdef DEBUG_STEP1C
							printContigClassification(stdout, preClassFlagK, ' ');
							printf("contig %d vs contig %d [%d, %d] [%d,%d]\n", contig_k->property.contigID, contig_j->property.contigID, covIvlK[0],
									covIvlK[numCovIvlK * 3 - 2], covIvlJ[0], covIvlJ[numCovIvlJ * 3 - 2]);

							printContigClassification(stdout, preClassFlagJ, ' ');
							printf("contig %d vs contig %d [%d, %d] [%d,%d]\n", contig_j->property.contigID, contig_k->property.contigID, covIvlJ[0],
									covIvlJ[numCovIvlJ * 3 - 2], covIvlK[0], covIvlK[numCovIvlK * 3 - 2]);
#endif

							if (contig_k->numReadRelations == contig_k->maxReadRelations)
							{
								contig_k->maxReadRelations = (int) (contig_k->maxReadRelations * 1.2) + 5;
								contig_k->readRelations = (ReadRelation*) realloc(contig_k->readRelations, sizeof(ReadRelation) * contig_k->maxReadRelations);
							}

							int num = contig_k->numReadRelations;

							contig_k->property.rflag |= preClassFlagK;
							contig_k->readRelations[num].flag = preClassFlagK;
							contig_k->readRelations[num].numCoveredIntervals = numCovIvlK;
							contig_k->readRelations[num].coveredIntervals = covIvlK;
							contig_k->readRelations[num].correspID = contig_j->property.contigID;
							contig_k->readRelations[num].nJointReads = numComReads;
							contig_k->numReadRelations++;

							// store corresponding info for contig J

							if (contig_j->numReadRelations == contig_j->maxReadRelations)
							{
								contig_j->maxReadRelations = (int) (contig_j->maxReadRelations * 1.2) + 5;
								contig_j->readRelations = (ReadRelation*) realloc(contig_j->readRelations, sizeof(ReadRelation) * contig_j->maxReadRelations);
							}

							num = contig_j->numReadRelations;

							contig_j->property.rflag |= preClassFlagJ;
							contig_j->readRelations[num].flag = preClassFlagJ;
							contig_j->readRelations[num].numCoveredIntervals = numCovIvlJ;
							contig_j->readRelations[num].coveredIntervals = covIvlJ;
							contig_j->readRelations[num].correspID = contig_k->property.contigID;
							contig_j->readRelations[num].nJointReads = numComReads;

							contig_j->numReadRelations++;
						}
					}
				}
			}
		}
	}

// cleanup
	free(contigIDAndLength);

	free(contigJReads);
	free(contigKReads);
	free(contigIntersectionJKReads);

	free(contigJBegRange);
	free(contigJEndRange);
	free(contigKEndRange);
	free(contigKBegRange);
}

int compareInt(const void * a, const void * b)
{
	return (*(int*) a - *(int*) b);
}

int contained(int ab, int ae, int bb, int be)
{
	if (ab >= bb && ae <= be)
	{
		return 1;
	}

	return 0;
}

static int getRepeatBasesOfContigRange(AnalyzeContext *ctx, Overlap *ovl, int read)
{
	if (ctx->corContigRepeats_track == NULL)
	{
		return 0;
	}

	assert(ovl->aread == read || ovl->bread == read);

	int bLen = ovl->path.bepos - ovl->path.bbpos;

	// get repeats track
	track_anno* rep_anno = ctx->corContigRepeats_track->anno;
	track_data* rep_data = ctx->corContigRepeats_track->data;

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
			if (rEnd > ovl->path.aepos)
				break;
		}
		else
		{
			if (ovl->flags & OVL_COMP)
			{
				nrep += intersect(bLen - ovl->path.bepos, bLen - ovl->path.bbpos, rBeg, rEnd);
				if (rEnd > bLen - ovl->path.bbpos)
					break;
			}
			else
			{
				nrep += intersect(ovl->path.bbpos, ovl->path.bepos, rBeg, rEnd);
				if (rEnd > ovl->path.bepos)
					break;
			}
		}
		rb += 2;
	}
	return nrep;
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

int analyzeContigVsContigOverlaps(void* _ctx, Overlap* ovls, int novl)
{
	AnalyzeContext* actx = (AnalyzeContext*) _ctx;

	if (actx->contigs[ovls->aread].property.cflag & CLASS_CONTIG_DISCARD)
		return 1;

	if (actx->VERBOSE)
		printf("BEGIN -- Analyze overlaps for contig: %d numOvls: %d\n", ovls->aread, novl);

	Contig *conA = actx->contigs + ovls->aread;

	// speed up things for long contigs

	int i, j, k;
	j = k = 0;
	while (j < novl)
	{
		while (k < novl - 1 && ovls[j].bread == ovls[k + 1].bread)
			k++;

		Overlap *o1, *o2;
		o1 = ovls + j;
		int fail = 0;
		int properBegA = 1;
		int properEndA = 1;
		int properBegB = 1;
		int properEndB = 1;

		int overlapBasesA = 0;
		int overlapBasesB = 0;
		int properGapLen = 1;

		if (o1->path.abpos > actx->nFuzzBases)
		{
			properBegA = 0;
		}
		if (o1->path.bbpos > actx->nFuzzBases)
		{
			properBegB = 0;
		}

		if (ovls[j].path.aepos + actx->nFuzzBases < DB_READ_LEN(actx->corContigDB, ovls[j].aread))
		{
			properEndA = 0;
		}

		if (ovls[j].path.bepos + actx->nFuzzBases < DB_READ_LEN(actx->corContigDB, ovls[j].bread))
		{
			properEndB = 0;
		}

		overlapBasesA = o1->path.aepos - o1->path.abpos;
		overlapBasesB = o1->path.bepos - o1->path.bbpos;

		for (i = j + 1; i <= k; i++)
		{
			o2 = ovls + i;

			overlapBasesA += o2->path.aepos - o2->path.abpos;
			overlapBasesB += o2->path.bepos - o2->path.bbpos;

			// validate oriantation
			if ((o1->flags & OVL_COMP) != (o2->flags & OVL_COMP))
			{
				fail = 1;
				break;
			}

			// validate direction
			if (o1->path.abpos > o2->path.abpos || o1->path.bbpos > o2->path.bbpos || o1->path.aepos > o2->path.aepos || o1->path.bepos > o2->path.bepos)
			{
				fail = 2;
				break;
			}

			// validate overlap between overlaps - A-coordinates
			int len = o1->path.aepos - o2->path.abpos;
			if (len > 0) // there is an overlap of neighboring alignments
			{
				if (len > MIN(o1->path.aepos - o1->path.abpos, o2->path.aepos - o2->path.abpos))
				{
					fail = 3;
					break;
				}
				overlapBasesA -= len;
			}
			else // there is a gap between neighboring alignments
			{
				if (abs(len) > actx->nFuzzBases)
				{
					fail = 4;
					break;
				}
			}
			// validate overlap between overlaps - A-coordinates
			len = o1->path.bepos - o2->path.bbpos;
			if (len > 0) // there is an overlap of neighboring alignments
			{
				if (len > MIN(o1->path.bepos - o1->path.bbpos, o2->path.bepos - o2->path.bbpos))
				{
					fail = 5;
					break;
				}
				overlapBasesB -= len;
			}
			else // there is a gap between neighboring alignments
			{
				if (abs(len) > actx->nFuzzBases)
				{
					fail = 6;
					break;
				}
			}

			o1 = o2;
		}

		int validContainment = 0;
		int validBridge = 0;

		if(properGapLen)
		{
			if(MAX(overlapBasesA, overlapBasesB) >= (int) (actx->contByContigs_CoveredLenPerc / 100.0 * MIN(DB_READ_LEN(actx->corContigDB, ovls[j].aread), DB_READ_LEN(actx->corContigDB, ovls[j].bread))))
			{
				validContainment = 1;
			}

			//      contigA         ----------------
			//      contigB	---------
			//			min 50Kb overlap, both overhangs min 100Kb
			else if(properBegA && !properEndA && properEndB && !properBegB && MAX(overlapBasesA, overlapBasesB) >= MIN(50000, 3*actx->nFuzzBases) && ovls[k].path.aepos + 100000 < DB_READ_LEN(actx->corContigDB, ovls[j].aread)
					&& ovls[j].path.bbpos - 100000 > 0)
			{
				validBridge = 1;
			}
			//      contigA         ----------------
			//      								contigB			------------
			//			min 50Kb overlap, both overhangs min 100Kb
			else if(!properBegA && properEndA && !properEndB && properBegB && MAX(overlapBasesA, overlapBasesB) >= MIN(50000, 3*actx->nFuzzBases) && ovls[j].path.abpos - 100000 > 0 && ovls[k].path.bepos + 100000 < DB_READ_LEN(actx->corContigDB, ovls[j].bread))
			{
				validBridge = 1;
			}
		}

		if (!fail && (validBridge || validContainment))
		{
			Contig *conB = actx->contigs + ovls[j].bread;

			if (conA->numContigRelations == conA->maxContigRelations)
			{
				conA->maxContigRelations = conA->maxContigRelations * 1.2 + 5;
				conA->contigRelations = (ContigRelation*) realloc(conA->contigRelations, sizeof(ContigRelation) * conA->maxContigRelations);
				bzero(conA->contigRelations + conA->numContigRelations, sizeof(ContigRelation) * (conA->maxContigRelations - conA->numContigRelations));
			}

			ContigRelation *crel = conA->contigRelations + conA->numContigRelations;
			crel->corContigIdx = conB->property.contigID;
			crel->numPos = (k - j + 1);
			crel->abpos = (int*) malloc(sizeof(int) * (crel->numPos) * 4);
			crel->aepos = crel->abpos + crel->numPos;
			crel->bbpos = crel->aepos + crel->numPos;
			crel->bepos = crel->bbpos + crel->numPos;

			int cumAaln = ovls[j].path.aepos - ovls[j].path.abpos;

			// add first coordinates
			crel->abpos[0] = ovls[j].path.abpos;
			crel->aepos[0] = ovls[j].path.aepos;

			// add further coordinates, but subtract ovlapping neighboring overlaps
			for (i = 1; i < crel->numPos; i++)
			{
				crel->abpos[i] = ovls[j+i].path.abpos;
				crel->aepos[i] = ovls[j+i].path.aepos;

				cumAaln += crel->aepos[i] - crel->abpos[i];
				if (crel->abpos[i] < crel->aepos[i-1])
				{
					cumAaln -=  (crel->aepos[i-1] - crel->abpos[i]);
				}
			}

			int cumBaln;
			if (ovls[k].flags & OVL_COMP)
			{
				cumBaln = ovls[j].path.bepos - ovls[j].path.bbpos;

				crel->bbpos[0] = conB->property.len - ovls[j].path.bbpos;
				crel->bepos[0] = conB->property.len - ovls[j].path.bepos;

				for (i = 1; i < crel->numPos; i++)
				{
					crel->bbpos[i] = conB->property.len - ovls[j + i].path.bbpos;
					crel->bepos[i] = conB->property.len - ovls[j + i].path.bepos;

					cumBaln += crel->bbpos[i] - crel->bepos[i];

					if (crel->bepos[i-1] < crel->bbpos[i])
					{
						cumBaln -= (crel->bbpos[i] - crel->bepos[i-1]);
					}
				}
			}
			else
			{
				cumBaln = ovls[j].path.bepos - ovls[j].path.bbpos;

				crel->bbpos[0] = ovls[j].path.bbpos;
				crel->bepos[0] = ovls[j].path.bepos;

				for (i = 1; i < crel->numPos; i++)
				{
					crel->bbpos[i] = ovls[j + i].path.bbpos;
					crel->bepos[i] = ovls[j + i].path.bepos;

					cumBaln += crel->bepos[i] - crel->bbpos[i];
					if (crel->bbpos[i] < crel->bepos[i-1])
					{
						cumBaln -=  (crel->bepos[i-1] - crel->bbpos[i]);
					}
				}
			}

			printf("  ADD ContigRelation %d (l %d) vs %d (l %d) nOvls: %d alignedA %d / %d alignedB %d / %d (validBridge %d, validContainment %d)\n", conA->property.contigID, conA->property.len, conB->property.contigID, conB->property.len, crel->numPos, cumAaln, conA->property.len, cumBaln, conB->property.len, validBridge, validContainment);
		}
		else
		{
			Contig *conB = actx->contigs + ovls[j].bread;
			printf("  FAILED [%d,b%d,c%d,g%d] ContigRelation %d (l %d) vs %d (l %d) nOvls: %d\n", fail, validBridge, validContainment,properGapLen, conA->property.contigID, conA->property.len, conB->property.contigID, conB->property.len, (k - j + 1));
		}
		k++;
		j = k;
	}

	printf("END -- Analyze overlaps for contig: %d numOvls: %d\n", ovls->aread, novl);

	return 1;
}

void preClassifyContigsByContigOverlaps(AnalyzeContext *actx)
{
	int i, j;

	for (i = 0; i < actx->numContigs; i++)
	{
		printf("%5d contig: %d, len: %d\n", i, actx->contigs[i].property.contigID, actx->contigs[i].property.len);

		Contig *contigA = actx->contigs + i;

		if (contigA->property.cflag & CLASS_CONTIG_DISCARD)
			continue;

		// no valid overlap chain found --> its unique
		if (contigA->numContigRelations == 0)
		{
			contigA->property.rflag |= (REL_CONTIG_UNIQUE);
			continue;
		}

		// analyze all overlap groups
		int overAllNotClassified = 1;

		/// TODO: only analyze containment relationships for now
		/// later include overlapping contigs --> trim back !!!

		for(j=0; j< contigA->numContigRelations; j++)
		{
			ContigRelation *cRel = contigA->contigRelations;

			Contig *contigB = actx->contigs + cRel->corContigIdx;

			printf("check Relationship %d vs %d\n", contigA->property.contigID, contigB->property.contigID);

			int covBasesInA, covBasesInB;
			covBasesInA = covBasesInB = 0;
			int k;
			covBasesInA = cRel->aepos[0] - cRel->abpos[0];
			covBasesInB = abs(cRel->bepos[0] - cRel->bbpos[0]);

			int comp = (cRel->bbpos[0] > cRel->bepos[0]) ? 1 : 0;

			for (k=1; k < cRel->numPos; k++)
			{
				covBasesInA += cRel->aepos[i] - cRel->abpos[i];
				if (cRel->abpos[i] < cRel->aepos[i-1])
				{
					covBasesInA -= (cRel->aepos[i-1] - cRel->abpos[i]);
				}

				if (comp)
				{
					covBasesInB += cRel->bbpos[i] - cRel->bepos[i];
					if (cRel->bepos[i-1] < cRel->bbpos[i])
					{
						covBasesInB -= (cRel->bbpos[i] - cRel->bepos[i-1]);
					}
				}
				else
				{
					covBasesInB += cRel->bepos[i] - cRel->bbpos[i];
					if (cRel->bbpos[i] < cRel->bepos[i-1])
					{
						covBasesInB -= (cRel->bepos[i-1] - cRel->bbpos[i]);
					}
				}
			}

			int leftUnalnBasesInA, rightUnalnBasesInA;
			int leftUnalnBasesInB, rightUnalnBasesInB;

			leftUnalnBasesInA = cRel->abpos[0];
			rightUnalnBasesInA = contigA->property.len - cRel->aepos[cRel->numPos - 1];

			if(comp)
			{
				leftUnalnBasesInB = contigB->property.len - cRel->bbpos[0];
				rightUnalnBasesInB = cRel->bepos[cRel->numPos - 1];
			}
			else
			{
				leftUnalnBasesInB = cRel->bbpos[0];
				rightUnalnBasesInB = contigB->property.len - cRel->bepos[cRel->numPos - 1];
			}


			// check if there is a containment relation

			if(covBasesInA * 100.0 / contigA->property.len >= actx->contByContigs_CoveredLenPerc)
			{
				cRel->flag |= REL_CONTIG_IS_ALT;
				contigA->property.rflag |= REL_CONTIG_IS_ALT;
				printf("CONTIGA %d IS ALT [covBases %d, leftUnA %d, right Una %d]", contigA->property.contigID, covBasesInA, leftUnalnBasesInA, rightUnalnBasesInA);
			}

			if(covBasesInB * 100.0 / contigB->property.len >= actx->contByContigs_CoveredLenPerc)
			{
				cRel->flag |= REL_CONTIG_HAS_ALT;
				contigA->property.rflag |= REL_CONTIG_HAS_ALT;
				contigB->property.rflag |= REL_CONTIG_IS_ALT;
				printf("CONTIGB %d IS ALT [covBases %d, leftUnA %d, right Una %d]", contigB->property.contigID, covBasesInB, leftUnalnBasesInB, rightUnalnBasesInB);
			}
			// TODO check for putative joins
			// TODO check for putative misjoins
		}
	}
}

int createOutDir(char *out)
{
	int res = 0;

	if (out)
	{
		struct stat s;

		int err = stat(out, &s);

		if (err == -1)
		{
			if (errno == ENOENT)
			{
				err = mkdir(out, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				if (err != 0)
				{
					fprintf(stderr, "Cannot create output directory: %s\n", out);
					res = 1;
				}
			}
			else
			{
				fprintf(stderr, "Cannot create output directory: %s\n", out);
				res = 1;
			}
		}
		else
		{
			if (!S_ISDIR(s.st_mode))
			{
				fprintf(stderr, "Output directory name: \"%s\" exist - but its not a directory\n", out);
				res = 1;
			}
		}
	}
	return res;
}

static void writeContigInfoFile(AnalyzeContext *actx, Contig* contig, FILE* contigFile, int cBegPos, int cEndPos, bit *allReads, int *contigCov)
{
	ba_all_assign(allReads, DB_NREADS(actx->patchedReadDB), FALSE);
	bzero(contigCov, sizeof(int) * DB_READ_MAXLEN(actx->corContigDB));

// determine all overlapping reads and coverage of contig
	int i, j, k;
	for (i = 0; i < contig->numcReads; i++)
	{
		ContigRead *cread = contig->cReads + i;

		if (intersect(cread->ovlReads->cBeg, cread->ovlReads->cEnd, cBegPos, cEndPos))
		{
			for (j = 0; j < cread->numOvlReads; j++)
			{
				if (intersect(cread->ovlReads[j].cBeg, cread->ovlReads[j].cEnd, cBegPos, cEndPos))
				{
					ba_assign(allReads, abs(cread->ovlReads[j].patchedID), TRUE);

					int covBeg = MAX(cread->ovlReads[j].cBeg, cBegPos);
					int covEnd = MIN(cread->ovlReads[j].cEnd, cEndPos);
					for (k = covBeg; k < covEnd; k++)
						contigCov[k]++;
				}
			}
		}
	}

// print all patched reads that overlap with current contig
	fprintf(contigFile, "PATCH_READS");
	for (i = 0; i < DB_NREADS(actx->patchedReadDB); i++)
		if (ba_value(allReads, i))
			fprintf(contigFile, " %d", i);
	fprintf(contigFile, "\n");

// print all source reads that should overlap with current contig
	if (actx->patchedReadSource_track)
	{
		// fixed source read ids from contig
		track_anno* sreadanno = actx->patchedReadSource_track->anno;
		track_data* sreaddata = actx->patchedReadSource_track->data;

		fprintf(contigFile, "SOURCE_READS");
		for (i = 0; i < DB_NREADS(actx->patchedReadDB); i++)
			if (ba_value(allReads, i))
			{
				track_anno rb = sreadanno[i] / sizeof(track_data);
				track_anno re = sreadanno[i + 1] / sizeof(track_data);

				assert(rb < re);

				fprintf(contigFile, " %d", sreaddata[rb]);
			}
		fprintf(contigFile, "\n");
	}

// print coverage histogram
	int binSize = 1000;
	fprintf(contigFile, "COV_HIST_B%d", binSize);

	int overallCov = 0;
	for (i = cBegPos; i + binSize < cEndPos; i++)
		overallCov += contigCov[i];

	fprintf(contigFile, " %.1f", overallCov * 1.0 / (cEndPos - cBegPos));

	int cumCov;
	for (i = cBegPos; i + binSize < cEndPos; i += binSize)
	{
		cumCov = 0;
		for (k = i; k < i + binSize; k++)
			cumCov += contigCov[k];
		fprintf(contigFile, " %.1f", cumCov * 1.0 / binSize);
	}
	if (i < cEndPos)
	{
		cumCov = 0;
		binSize = cEndPos - i;
		for (k = i; k < cEndPos; k++)
			cumCov += contigCov[k];
		fprintf(contigFile, " %.1f", cumCov * 1.0 / binSize);
	}
	fprintf(contigFile, "\n");

// repeat annotation:
	fprintf(contigFile, "REPEAT_%s", actx->corContigRepeats_track->name);
	fprintf(contigFile, " %d", getRepeatBasesFromInterval(actx, 1, contig->property.contigID, cBegPos, cEndPos));
	{
		track_anno* rep_anno = actx->corContigRepeats_track->anno;
		track_data* rep_data = actx->corContigRepeats_track->data;

		track_anno rb, re;
		// repeat bases in a-read
		rb = rep_anno[contig->property.contigID] / sizeof(track_data);
		re = rep_anno[contig->property.contigID + 1] / sizeof(track_data);

		while (rb < re)
		{
			int rBeg = rep_data[rb];
			int rEnd = rep_data[rb + 1];

			if (rEnd < cBegPos)
			{
				continue;
				rb += 2;
			}
			if (rBeg > cEndPos)
				break;

			if (rBeg < cBegPos)
				rBeg = cBegPos;

			if (rEnd > cEndPos)
				rEnd = cEndPos;

			fprintf(contigFile, " %d,%d", rBeg, rEnd);

			rb += 2;
		}
	}
	fprintf(contigFile, "\n");

// report splits for raw contigs
//	if (contig->property.len == cEndPos - cBegPos && contig->numSplits && !(contig->split->type & SPLIT_IGNORE))
//	{
//		fprintf(contigFile, "SPLITS");
//
//		for (i = 0; i < contig->numSplits; i++)
//		{
//			if (contig->split[i].type & SPLIT_IGNORE)
//				break;
//			fprintf(contigFile, " %d,%d,%d", contig->split[i].type, contig->split[i].contigPos1, contig->split[i].contigPos2);
//		}
//		fprintf(contigFile, "\n");
//	}
}

// todo check if cEndPos is really exclusive
void writeFasta(AnalyzeContext *actx, Contig* contig, FILE* contigFile, int splitIdx, int cBegPos, int cEndPos)
{

	int i;
	int WIDTH = 100;

	int path = getPathID(actx, contig->property.contigID);
	int end1, end2;
	getContigsEndIDs(actx, contig->property.contigID, &end1, &end2);
	char *InFastaName = getContigFastaFileName(actx, contig);

///// create header
// new name = fasta file name + Contig idx
	fprintf(contigFile, ">%s_%d", InFastaName, contig->property.contigID);

///// add some track info
// path track
	fprintf(contigFile, " path=%d", path);
// length track
	fprintf(contigFile, " length=%d", cEndPos - cBegPos);

	if (splitIdx == 0) // otherwise the ends can be invalid
	{
		// ends track
		fprintf(contigFile, " ends=%d,%d", end1, end2);
		// repeat track
		//fprintf(contigFile, " repeat=%d,%d,%d", contig->repBases, contig->repBasesTipLeft, contig->repBasesTipRight);
	}
	else
	{
		fprintf(contigFile, " splitIdx=%d", splitIdx);
	}
// reads track
	if (contig->numcReads)
	{
		fprintf(contigFile, " reads=");
		int curPos = 0;
		for (i = 0; i < contig->numcReads; i++)
		{
			if (curPos < cBegPos)
			{
				curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);
				continue;
			}
			if (curPos > cEndPos)
				break;

			// first read of split contig
			if (cBegPos > curPos && cBegPos < curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg))
			{
				if (contig->cReads[i].ovlReads->beg < contig->cReads[i].ovlReads->end) // forward orientation
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg + (cBegPos - curPos),
							contig->cReads[i].ovlReads->end);
				}
				else
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg - (cBegPos - curPos),
							contig->cReads[i].ovlReads->end);
				}
			}
			// last read of split contig
			else if (curPos < cEndPos && cEndPos < curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg))
			{
				if (contig->cReads[i].ovlReads->beg < contig->cReads[i].ovlReads->end) // forward orientation
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg,
							contig->cReads[i].ovlReads->end - (curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg) - cEndPos));
				}
				else
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg,
							contig->cReads[i].ovlReads->end + (curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg) - cEndPos));
				}
				break;
			}
			else
			{
				fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg, contig->cReads[i].ovlReads->end);
			}

			curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);
			if (curPos < cEndPos)
				fprintf(contigFile, ",");

		}
	}
// sreads track
	if (actx->corContigRawReads_track)
	{
		// fixed source read ids from contig
		track_anno* sreadanno = actx->corContigRawReads_track->anno;
		track_data* sreaddata = actx->corContigRawReads_track->data;

		track_anno rb = sreadanno[contig->property.contigID] / sizeof(track_data);
		track_anno re = sreadanno[contig->property.contigID + 1] / sizeof(track_data);

		assert(rb < re && (int) (re - rb) == contig->numcReads);

		fprintf(contigFile, " sreads=");
		int curPos = 0;
		for (i = 0; i < contig->numcReads; i++, rb++)
		{
			if (curPos < cBegPos)
			{
				curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);
				continue;
			}
			if (curPos > cEndPos)
				break;

			fprintf(contigFile, "%d", sreaddata[rb]);
			curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);

			if (curPos < cEndPos)
				fprintf(contigFile, ",");
		}

	}
// classification tracks
// 1)  based on fixed read overlaps
//	if (contig->gClassificFlag & CONTIG_HAS_CONTAINED)
//	{
//		int multi = 0;
//		ContigReadClassification *gclass;
//		for (i = 0; i < contig->numGClassific; i++)
//		{
//			gclass = contig->gClassific + i;
//
//			if (gclass->flag & CONTIG_HAS_CONTAINED)
//			{
//				if (multi)
//				{
//					fprintf(contigFile, ",%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals * 3 - 2]);
//				}
//				else
//				{
//					fprintf(contigFile, " gContains=%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0],
//							gclass->coveredIntervals[gclass->numCoveredIntervals * 3 - 2]);
//					multi = 1;
//				}
//			}
//		}
//	}
//	if (contig->gClassificFlag & CONTIG_IS_CONTAINED)
//	{
//		int multi = 0;
//		ContigReadClassification *gclass;
//		for (i = 0; i < contig->numGClassific; i++)
//		{
//			gclass = contig->gClassific + i;
//
//			if (gclass->flag & CONTIG_IS_CONTAINED)
//			{
//				if (multi)
//				{
//					fprintf(contigFile, ",%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals * 3 - 2]);
//				}
//				else
//				{
//					fprintf(contigFile, " gContainedIn=%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0],
//							gclass->coveredIntervals[gclass->numCoveredIntervals * 3 - 2]);
//					multi = 1;
//				}
//			}
//		}
//	}
//// 2)  based on contig overlap groups
//	if (contig->flag & CONTIG_HAS_CONTAINED)
//	{
//		int multi = 0;
//		OverlapGroup *ogr;
////		for (i = 0; i < contig->numOvlGrps; i++)
////		{
////			ogr = contig->ovlGrps[i];
////
////			if ((ogr->flag & CONTIG_HAS_CONTAINED) && contig->property.len > DB_READ_LEN(actx->corContigDB, ogr->bread))
////			{
////				if (multi)
////				{
////					fprintf(contigFile, ",%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
////				}
////				else
////				{
////					fprintf(contigFile, " cContains=%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
////					multi = 1;
////				}
////			}
////		}
//	}
//	if (contig->flag & CONTIG_IS_CONTAINED)
//	{
//		int multi = 0;
//		OverlapGroup *ogr;
////		for (i = 0; i < contig->numOvlGrps; i++)
////		{
////			ogr = contig->ovlGrps[i];
////
////			if ((ogr->flag & CONTIG_IS_CONTAINED) && contig->property.len < DB_READ_LEN(actx->corContigDB, ogr->bread))
////			{
////				if (multi)
////				{
////					fprintf(contigFile, ",%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
////				}
////				else
////				{
////					fprintf(contigFile, " cContainedIn=%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
////					multi = 1;
////				}
////			}
////		}
//	}
//
//	if (contig->flag & CONTIG_AL50PCTCOVERED)
//	{
//		int multi = 0;
//
//		for (i = 0; i < actx->numContigs; i++)
//		{
////			if ((int) contig->pctCoveredBasesInOtherContigs[i] >= 50)
////			{
////				if (multi)
////				{
////					fprintf(contigFile, ",%d,%d", i, (int) contig->pctCoveredBasesInOtherContigs[i]);
////				}
////				else
////				{
////					fprintf(contigFile, " cCoveredIn=%d,%d", i, (int) contig->pctCoveredBasesInOtherContigs[i]);
////					multi = 1;
////				}
////			}
//		}
//	}

// todo sort splits according positions ?
// todo write out only real splits, i.e. skip heterozygoues ones, and print out only one valid split position either left or, or right or sometimes both
//	if (contig->numSplits && !(contig->split->type & SPLIT_IGNORE))
//	{
//		int multi = 0;
//		// print absolute contig Positions
//		for (i = 0; i < contig->numSplits; i++)
//		{
//			if (multi)
//			{
//				fprintf(contigFile, ",%d,%d,%d", contig->split[i].type, contig->split[i].contigPos1, contig->split[i].contigPos2);
//			}
//			else
//			{
//				fprintf(contigFile, " splitAbsPos=%d,%d,%d", contig->split[i].type, contig->split[i].contigPos1, contig->split[i].contigPos2);
//				multi = 1;
//			}
//		}
//		multi = 0;
//		// print split reads and their positions
//		for (i = 0; i < contig->numSplits; i++)
//		{
//			if (multi)
//			{
//				fprintf(contigFile, ",%d,%d,%d,%d,%d,%d", contig->split[i].leftProperRead, contig->split[i].leftProperReadPos1,
//						contig->split[i].leftProperReadPos2, contig->split[i].rightProperRead, contig->split[i].rightProperReadPos1,
//						contig->split[i].rightProperReadPos2);
//			}
//			else
//			{
//				fprintf(contigFile, " split=%d,%d,%d,%d,%d,%d", contig->split[i].leftProperRead, contig->split[i].leftProperReadPos1,
//						contig->split[i].leftProperReadPos2, contig->split[i].rightProperRead, contig->split[i].rightProperReadPos1,
//						contig->split[i].rightProperReadPos2);
//				multi = 1;
//			}
//		}
//	}

	fprintf(contigFile, "\n");

///// add sequence
	Load_Read(actx->corContigDB, contig->property.contigID, actx->readSeq, 1);

	for (i = cBegPos; i + WIDTH < cEndPos; i += WIDTH)
		fprintf(contigFile, "%.*s\n", WIDTH, actx->readSeq + i);
	if (i < cEndPos)
		fprintf(contigFile, "%.*s\n", cEndPos - i, actx->readSeq + i);

	fflush(contigFile);
}

static int cmpInts(const void * a, const void * b)
{
	return (*(int*) b - *(int*) a);
}

//static int cmpSplitEvents(const void * a, const void * b)
//{
//	SplitEvent *sa = (SplitEvent*) a;
//	SplitEvent *sb = (SplitEvent*) b;
//
//	if ((sa->type & SPLIT_IGNORE) && (sb->type & SPLIT_IGNORE))
//		return (sa->contigPos1 - sb->contigPos1);
//	if (sa->type & SPLIT_IGNORE)
//		return 1;
//	if (sb->type & SPLIT_IGNORE)
//		return -1;
//	return (sa->contigPos1 - sb->contigPos1);
//}

/*
 void printFinalAltContigClassification(FILE *out, FinalContigAltClass facc)
 {
 switch (facc)
 {
 case BUBBLE_ALT:
 fprintf(out, "BUBBLE_ALT");
 break;
 case BUBBLE_HUGEDIF:
 fprintf(out, "BUBBLE_HUGEDIF");
 break;
 case BUBBLE_LOWCOVEREDREPEAT:
 fprintf(out, "BUBBLE_LOWCOVEREDREPEAT");
 break;
 case BUBBLE_REPEAT:
 fprintf(out, "BUBBLE_REPEAT");
 break;
 case BUBBLE_WEIRD:
 fprintf(out, "BUBBLE_WEIRD");
 break;
 case SPUR_ALT:
 fprintf(out, "SPUR_ALT");
 break;
 case SPUR_ATCONTIGEND:
 fprintf(out, "SPUR_ATCONTIGEND");
 break;
 default:
 fprintf(out, "NONE");
 break;
 }
 }

 void printFinalContigClassification(FILE *out, FinalContigClass fcc)
 {
 switch (fcc)
 {
 case HAPLOID:
 fprintf(out, "HAPLOID");
 break;
 case CONTAINED_BY_BOTH:
 fprintf(out, "CONTAINED_BY_BOTH");
 break;
 case CONTAINED_BY_BOTHSTRANGE:
 fprintf(out, "CONTAINED_BY_BOTHSTRANGE");
 break;
 case CONTAINED_BY_CONTIG:
 fprintf(out, "CONTAINED_BY_CONTIG");
 break;
 case CONTAINED_BY_FIXREAD:
 fprintf(out, "CONTAINED_BY_FIXREAD");
 break;
 case CONTAINED_BY_COVEREDBASES:
 fprintf(out, "CONTAINED_BY_COVEREDBASES");
 break;
 default:
 fprintf(out, "UNKNOWN");
 break;
 }
 }

 char getClassification(FinalContigClass fcc)
 {
 switch (fcc)
 {
 case HAPLOID:
 return 'H';
 case CONTAINED_BY_BOTH:
 return 'A';
 case CONTAINED_BY_BOTHSTRANGE:
 return 'W';
 case CONTAINED_BY_CONTIG:
 return 'C';
 case CONTAINED_BY_FIXREAD:
 return 'R';
 case CONTAINED_BY_COVEREDBASES:
 return 'B';
 default:
 return 'U';
 }
 }*/

//char* getSplitType(SplitEvent* se)
//{
//	static char* result = NULL;
//
//	if (result == NULL)
//		result = malloc(20);
//
//	result[0] = '\0';
//
//	if (se->type & SPLIT_IGNORE)
//	{
//		strcat(result, "I");
//	}
//	if (se->type & SPLIT_UNKNOWN)
//	{
//		strcat(result, "u");
//	}
//	if (se->type & SPLIT_MULTIREAD)
//	{
//		strcat(result, "M");
//	}
//	if (se->type & SPLIT_STARTEND)
//	{
//		strcat(result, "S");
//	}
//	if (se->type & SPLIT_FALSEJOINLEFT)
//	{
//		strcat(result, "f");
//	}
//	if (se->type & SPLIT_FALSEJOINRIGHT)
//	{
//		strcat(result, "F");
//	}
//	if (se->type & SPLIT_HETEROZYGLEFT)
//	{
//		strcat(result, "h");
//	}
//	if (se->type & SPLIT_HETEROZYGRIGHT)
//	{
//		strcat(result, "H");
//	}
//	if (se->type & SPLIT_DEADENDLEFT)
//	{
//		strcat(result, "d");
//	}
//	if (se->type & SPLIT_DEADENDRIGHT)
//	{
//		strcat(result, "D");
//	}
//
//	return result;
//}

//static int cmpGraphClassByPos(const void* a, const void* b)
//{
//	ContigReadClassification* c1 = (ContigReadClassification*) a;
//	ContigReadClassification* c2 = (ContigReadClassification*) b;
//
//	if (c1->flag & (CONTIG_UNCLASSIFIED | CONTIG_DISCARD | CONTIG_NORELATION | CONTIG_UNIQUE))
//		return 1;
//
//	if (c2->flag & (CONTIG_UNCLASSIFIED | CONTIG_DISCARD | CONTIG_NORELATION | CONTIG_UNIQUE))
//		return -1;
//
//	return (c1->coveredIntervals - c2->coveredIntervals);
//}

void rawClassification(AnalyzeContext *actx)
{
	/*
	 int i, j, k;
	 int end1_hap, end2_hap; // ends from touring
	 int path_hap;

	 Contig *contig_hap;

	 bit *contigVisited = ba_new(actx->numContigs);

	 int idx = 0;
	 int curStart = 0;
	 char *InFastaName;
	 int path_other, end1_other, end2_other;

	 while (curStart < actx->numContigs)
	 {
	 int process = 0;
	 int curNum = getNumberOfSequencesFromFastaFileName(actx, curStart);

	 int compoundNumAllBases = 0;
	 int compoundNumHapContigs = 0;
	 int compoundNumHapBases = 0;
	 int compoundNumALTcontigs = 0;
	 int compoundNumALTbases = 0;
	 int compoundNumREPcontigs = 0;
	 int compoundNumREPbases = 0;
	 int compoundNumWEIRDcontigs = 0;
	 int compoundNumWEIRDbases = 0;

	 printf("#--------- Result of Connected Compound %5d ---------\n", idx);
	 while (1)
	 {
	 contig_hap = NULL;
	 for (i = curStart; i < curStart + curNum; i++)
	 {
	 if (actx->contigs[i].fClass == HAPLOID && ba_value(contigVisited, i) == FALSE)
	 {
	 contig_hap = actx->contigs + i;
	 path_hap = getPathID(actx, contig_hap->property.contigID);
	 getContigsEndIDs(actx, contig_hap->property.contigID, &end1_hap, &end2_hap);
	 break;
	 }
	 }

	 if (contig_hap != NULL)
	 {
	 ba_assign(contigVisited, i, TRUE);
	 process++;

	 int numALTcontigs = 0;
	 int numALTbases = 0;
	 int numREPcontigs = 0;
	 int numREPbases = 0;
	 int numWEIRDcontigs = 0;
	 int numWEIRDbases = 0;

	 // sort overlap groups and contig classification according to begin position in current haploid contig

	 qsort(contig_hap->ovlGrps, contig_hap->numOvlGrps, sizeof(OverlapGroup*), cmpOVLGroupByPos);
	 qsort(contig_hap->gClassific, contig_hap->numGClassific, sizeof(ContigGraphClassification), cmpGraphClassByPos);

	 for (j = 0, k = 0; j < contig_hap->numOvlGrps && k < contig_hap->numGClassific;)
	 {
	 OverlapGroup *og = contig_hap->ovlGrps[j];
	 ContigGraphClassification *cgc = contig_hap->gClassific + k;

	 assert(og->bread >= curStart && og->bread <= curStart + curNum);
	 assert(cgc->correspID >= curStart && cgc->correspID <= curStart + curNum);

	 if ((og->flag & OVLGRP_DISCARD) || (cgc->flag & (CONTIG_NORELATION | CONTIG_DISCARD | CONTIG_UNCLASSIFIED | CONTIG_UNIQUE)))
	 break;

	 if (og->bread == cgc->correspID)
	 {
	 assert(intersect(og->first_abpos, og->last_aepos, cgc->bpos, cgc->epos) > 0);

	 Contig * contig_other = actx->contigs + og->bread;

	 if (ba_value(contigVisited, contig_other->property.contigID) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->property.contigID);
	 path_other = getPathID(actx, contig_other->property.contigID);
	 getContigsEndIDs(actx, contig_other->property.contigID, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_other, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->property.len);
	 // position in contig --> only meaningful in contained contigs
	 if (og->flag & OVLGRP_COMP)
	 fprintf(stdout, "\t%8d\t%8d", og->first_bepos, og->last_bbpos);
	 else
	 fprintf(stdout, "\t%8d\t%8d", og->first_bbpos, og->last_bepos);
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", og->first_abpos, og->last_aepos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);

	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->property.contigID, TRUE);
	 process++;

	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->property.len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->property.len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->property.len;
	 break;
	 }
	 }
	 else
	 {
	 printf("found ovl group and gclass between contigs %d and %d more then once!!\n", contig_hap->property.contigID, contig_other->property.contigID);
	 }
	 j++;
	 k++;
	 continue;
	 }
	 if (og->first_abpos < cgc->bpos)
	 {
	 Contig * contig_other = actx->contigs + og->bread;

	 if (ba_value(contigVisited, contig_other->property.contigID) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->property.contigID);
	 path_other = getPathID(actx, contig_other->property.contigID);
	 getContigsEndIDs(actx, contig_other->property.contigID, &end1_other, &end2_other);

	 //
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->property.len);
	 // position in contig --> only meaningful in contained contigs
	 if (og->flag & OVLGRP_COMP)
	 fprintf(stdout, "\t%8d\t%8d", og->first_bepos, og->last_bbpos);
	 else
	 fprintf(stdout, "\t%8d\t%8d", og->first_bbpos, og->last_bepos);
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", og->first_abpos, og->last_aepos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);

	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->property.contigID, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->property.len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->property.len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->property.len;
	 break;
	 }
	 }
	 else
	 {
	 printf("found ovl group and gclass between contigs %d and %d more then once!!\n", contig_hap->property.contigID, contig_other->property.contigID);
	 }
	 j++;

	 continue;
	 }
	 else
	 {
	 Contig * contig_other = actx->contigs + cgc->correspID;

	 if (ba_value(contigVisited, contig_other->property.contigID) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->property.contigID);
	 path_other = getPathID(actx, contig_other->property.contigID);
	 getContigsEndIDs(actx, contig_other->property.contigID, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->property.len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", cgc->bpos, cgc->epos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);
	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->property.contigID, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->property.len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->property.len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->property.len;
	 break;
	 }

	 }
	 else
	 {
	 printf("found ovl group and gclass between contigs %d and %d more then once!!\n", contig_hap->property.contigID, contig_other->property.contigID);
	 }
	 k++;
	 continue;
	 }
	 }
	 // check if there are remaining bubbles
	 for (; j < contig_hap->numOvlGrps; j++)
	 {
	 OverlapGroup *og = contig_hap->ovlGrps[j];

	 assert(og->bread >= curStart && og->bread <= curStart + curNum);
	 if (og->flag & OVLGRP_DISCARD)
	 break;
	 Contig * contig_other = actx->contigs + og->bread;

	 if (ba_value(contigVisited, contig_other->property.contigID) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->property.contigID);
	 path_other = getPathID(actx, contig_other->property.contigID);
	 getContigsEndIDs(actx, contig_other->property.contigID, &end1_other, &end2_other);

	 //
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->property.len);
	 // position in contig --> only meaningful in contained contigs
	 if (og->flag & OVLGRP_COMP)
	 fprintf(stdout, "\t%8d\t%8d", og->first_bepos, og->last_bbpos);
	 else
	 fprintf(stdout, "\t%8d\t%8d", og->first_bbpos, og->last_bepos);
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", og->first_abpos, og->last_aepos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);
	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->property.contigID, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->property.len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->property.len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->property.len;
	 break;
	 }
	 }
	 }
	 for (; k < contig_hap->numGClassific; k++)
	 {
	 ContigGraphClassification *cgc = contig_hap->gClassific + k;
	 assert(cgc->correspID >= curStart && cgc->correspID <= curStart + curNum);

	 if ((cgc->flag & (CONTIG_NORELATION | CONTIG_DISCARD | CONTIG_UNCLASSIFIED | CONTIG_UNIQUE)))
	 break;

	 Contig * contig_other = actx->contigs + cgc->correspID;

	 if (ba_value(contigVisited, contig_other->property.contigID) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->property.contigID);
	 path_other = getPathID(actx, contig_other->property.contigID);
	 getContigsEndIDs(actx, contig_other->property.contigID, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "   %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->property.len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", cgc->bpos, cgc->epos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);
	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->property.contigID, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->property.len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->property.len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->property.len;
	 break;
	 }
	 }
	 }

	 InFastaName = getFastaFileNameFromDB(actx, contig_hap->property.contigID);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "%c", getClassification(contig_hap->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_hap);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_hap->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_hap->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_hap->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_hap->repBases * 100.0 / contig_hap->property.len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_hap, end2_hap);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_hap->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_hap->fAltClass);
	 fprintf(stdout, "\n");
	 // haploid contig summary: print number of bubbles/spurs and percentage of bases according to contig length
	 fprintf(stdout, "# HAPLOID-Summary\n");
	 fprintf(stdout, "# HAPcontig  \tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", 1, contig_hap->property.len, 100.00);
	 fprintf(stdout, "# ALTcontig  \tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", numALTcontigs, numALTbases, numALTbases * 100.0 / contig_hap->property.len);
	 fprintf(stdout, "# REPcontig  \tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", numREPcontigs, numREPbases, numREPbases * 100.0 / contig_hap->property.len);
	 fprintf(stdout, "# WEIRDcontig\tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", numWEIRDcontigs, numWEIRDbases, numWEIRDbases * 100.0 / contig_hap->property.len);

	 // log splits
	 printf("splits: \n");
	 for (k = 0; k < contig_hap->numSplits; k++)
	 {
	 printSplitEvent(contig_hap, contig_hap->split + k);
	 }

	 printf("#-----------------------------------------------------\n");

	 // update compound numbers
	 compoundNumAllBases += (numALTbases + numREPbases + numWEIRDbases);
	 compoundNumALTcontigs += numALTcontigs;
	 compoundNumALTbases += numALTbases;
	 compoundNumREPcontigs += numREPcontigs;
	 compoundNumREPbases += numREPbases;
	 compoundNumWEIRDcontigs += numWEIRDcontigs;
	 compoundNumWEIRDbases += numWEIRDbases;
	 compoundNumHapContigs++;
	 compoundNumHapBases += contig_hap->property.len;
	 continue;
	 }
	 else // print out remaining contigs
	 {
	 if (process < curNum)
	 {

	 for (i = curStart; i < curStart + curNum; i++)
	 {
	 Contig * contig_other = actx->contigs + i;

	 if (ba_value(contigVisited, contig_other->property.contigID) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->property.contigID);
	 path_other = getPathID(actx, contig_other->property.contigID);
	 getContigsEndIDs(actx, contig_other->property.contigID, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "   %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->property.contigID);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->property.len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->property.len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);

	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->property.contigID, TRUE);
	 process++;

	 // update compound numbers
	 compoundNumAllBases += contig_other->property.len;

	 }
	 }
	 //printf("#-----------------------------------------------------\n");x
	 }
	 break;
	 }
	 }
	 printf("#--------- Summary of Connected Compound %5d ---------\n", idx);
	 fprintf(stdout, "# HAPcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumHapContigs, compoundNumHapBases, compoundNumHapBases * 100.0 / compoundNumAllBases);
	 fprintf(stdout, "# ALTcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumALTcontigs, compoundNumALTbases, compoundNumALTbases * 100.0 / compoundNumAllBases);
	 fprintf(stdout, "# REPcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumREPcontigs, compoundNumREPbases, compoundNumREPbases * 100.0 / compoundNumAllBases);
	 fprintf(stdout, "# WEIRDcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumWEIRDcontigs, compoundNumWEIRDbases,
	 compoundNumWEIRDbases * 100.0 / compoundNumAllBases);

	 curStart = curStart + curNum + 1;
	 idx++;
	 break;
	 }*/
}

/// evaluate all three relationship information: TourRelation, ContigAlnRelation, ReadIntersectionRelation
void classify(AnalyzeContext *actx)
{
	int i, j, k;
	int end1, end2; // ends from touring
	int path;

	int classified = 0;
	for (i = 0;  i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs_sorted[i];

		printf("Classify contig %d l %d nRel (T: %d C: %d R: %d)", contig->property.contigID, contig->property.len, contig->numTourRelations, contig->numContigRelations, contig->numReadRelations);

		// check for ALT contigs first
		if (contig->property.rflag & (REL_TOUR_IS_ALT))
		{
			int tourRelIdx = -1;
			int contRelIdx = -1;
			int readRelIdx = -1;

			for (j=0; j<contig->numTourRelations; j++)
			{
				TourRelation *tRel = contig->tourRelations + j;
				if (tRel->flag & REL_TOUR_IS_ALT)
				{
					assert(tourRelIdx == -1);
					tourRelIdx = j;
				}
			}

			assert(tourRelIdx != -1);

			for (j=0; j<contig->numContigRelations; j++)
			{
				ContigRelation *cRel = contig->contigRelations + j;
				if ((cRel->flag & REL_CONTIG_IS_ALT) && cRel->corContigIdx == contig->tourRelations[tourRelIdx].contigID1)
				{
					contRelIdx = j;
					break;
				}
			}

			for (j=0; j<contig->numReadRelations; j++)
			{
				ReadRelation *rRel = contig->readRelations + j;
				if ((rRel->flag & REL_READ_IS_ALT) && rRel->correspID == contig->tourRelations[tourRelIdx].contigID1)
				{
					readRelIdx = j;
					break;
				}
			}

			if (contig->tourRelations[tourRelIdx].contigID1 == contig->contigRelations[contRelIdx].corContigIdx)
			{
				contig->property.cflag |= (CLASS_CONTIG_ALT | CLASS_CONTIG_CLASSIFIED);
				contig->classif->correspTourRelation   = contig->tourRelations + tourRelIdx;
				contig->classif->correspContigRelation = contig->contigRelations + contRelIdx;
			}

			if (contig->tourRelations[tourRelIdx].contigID1 == contig->readRelations[readRelIdx].correspID)
			{
				contig->property.cflag |= CLASS_CONTIG_ALT | CLASS_CONTIG_CLASSIFIED;
				contig->classif->correspTourRelation = contig->tourRelations + tourRelIdx;
				contig->classif->correspReadRelation = contig->readRelations + readRelIdx;
			}
		}
		else
		{
			int primaryLenValid = contig->property.len > actx->minPrimContigLen;
			int primaryRepeatsValid = MAX(contig->property.repBasesFromContigLAS, contig->property.repBasesFromReadLAS) * 100.0/ contig->property.len <= actx->maxPrimContigRepeatPerc;
			int primareNumCreadsValid = contig->numcReads > actx->minPrimContigReads;



		}








		// check for containment
		if ((contig->property.rflag & (REL_TOUR_HAS_ALT | REL_TOUR_UNIQUE)) && (contig->property.rflag & (REL_CONTIG_HAS_ALT | REL_CONTIG_UNIQUE)) && (contig->property.rflag & (REL_READ_HAS_ALT | REL_READ_UNIQE)))
		{
			printf(" 3_PRIM");
			classified++;
		}
		else if ((contig->property.rflag & (REL_TOUR_IS_ALT)) && (contig->property.rflag & (REL_CONTIG_IS_ALT)) && (contig->property.rflag & (REL_READ_IS_ALT)))
		{
			printf(" 3_ALT");
			classified++;
		}

		printf("  -- classified (%d / %d)\n", classified, actx->numContigs);
	}
}

int checkJunctionCoverage(AnalyzeContext *actx, Contig *contig, int read) // returns true if everything is ok, false otherwise
{
	int i, j;

	int result = 1;

// todo pass variable, set by user ?
	int NUM_NODES_TOCHECK = 3; // number or reads to check before and after duplicated reads

	printf("checkJunctionCoverage for read %d and contig %d\n", read, contig->property.contigID);

	j = 0;
	while (j < contig->numcReads)
	{
		if (abs(contig->cReads[j].patchedID) == read)
			break;

		j++;
	}

	assert(j < contig->numcReads);

// try to start NUM_NODES_TOCHECK before read
	int from, to;
	from = to = j;
	int k = NUM_NODES_TOCHECK;
	int k_len_before = 0;

	while (k > 0)
	{
		if (from - 1 > 0)
			from--;
		if (to + 1 < contig->numcReads)
			to++;

		k--;
	}

	for (k = from; k <= to; k++)
	{
		ContigRead *cread = contig->cReads + k;

		if (cread->nComReadsWithNextContigRead == 0 || cread->lowCov)
			result = 0;

		printf("CP: read %d, # in %d, # out %d, #comRead %d, cPos [%d, %d], rPos [%d, %d], avgCov %.2f, covDrop %d\n", cread->patchedID, cread->nInReads,
				cread->nOutReads, cread->nComReadsWithNextContigRead, cread->patchedContigPosBeg, cread->patchedContigPosEnd, cread->ovlReads->beg,
				cread->ovlReads->end, cread->avgCov, cread->lowCov);

	}
	return result;
}

static void usage()
{
	fprintf(stderr,
			" [-v] [-clxsfeLNP <int>] [-d <dir>] [-rt <Track>] [-o <file>] -C <contigDB> <ContigOverlaps> -F <fixedReadDB> <fixedReadOverlaps> -D <correctedReadDB> \n");
	fprintf(stderr, "options: -v         ... verbose\n");
	fprintf(stderr, "         -x         ... min read length (default: 0)\n");
	fprintf(stderr, "         -l         ... min alignment length (default: 1000)\n");
	fprintf(stderr, "         -c         ... expected coverage\n");
	fprintf(stderr, "         -s         ... maximum spur/bubble length (default: %d)\n", DEF_SPUR_LEN);
	fprintf(stderr, "         -e         ... maximum contig end length (tips) to consider false joins (default: %d)\n", DEF_TIP_LEN);
	fprintf(stderr, "         -r         ... repeats track for contigs (default: repeats)\n");
	fprintf(stderr, "         -R         ... repeats track for patached reads (default: repeats)\n");
	fprintf(stderr, "         -t         ... trim track for fixed reads (default: none)\n");
	fprintf(stderr, "         -C DB OVL  ... contig database and contig overlaps, required to classify contigs\n");
	fprintf(stderr, "         -F DB OVL  ... fixed read database and fixed read overlaps, required to extract all b-reads for blasr mapping\n");
	fprintf(stderr, "         -D DB 		 ... corrected read database\n");
	fprintf(stderr, "         -d         ... write classified contigs ands stats file into -d directoryName (default: cwd)\n");
	fprintf(stderr, "         -o         ... write out filtered chained overlaps (default: cwd)\n");
	fprintf(stderr, "         -f         ... allow maximum of -f bases of structural variations between two contig overlaps of a chain, (default %d)\n", DEF_ARG_F);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "EPERIMENTAL - contig filter + classification options:\n");
	fprintf(stderr, "         -L         ... minimum primary contig length, (default %d)\n", DEF_ARG_L);
	fprintf(stderr, "         -N         ... minimum number of contig reads for a primary contig (default %d)\n", DEF_ARG_N);
	fprintf(stderr, "         -P         ... maximum repeat percentage of a primary contig, (default %d)\n", DEF_ARG_P);
}

int main(int argc, char* argv[])
{
	HITS_DB correctedContigDB;
	HITS_DB patchedReadDB;
	HITS_DB correctedReadDB;

	FILE* patchedReadLAS;
	FILE* correctedContigLAS;

	PassContext* contig_pctx;
	PassContext* patched_pctx;
	AnalyzeContext actx;

	char *patchedReadTrimTrack = NULL;
	char *patchedReadRepeatTrack = NULL;
	char *corContigRepeatsTrack = DEF_REPEATS_TRACK;

	bzero(&actx, sizeof(AnalyzeContext));
	opterr = 0;

	char *filteredContigOvlsName = NULL;
	FILE *filteredContigOvlsFile = NULL;

	actx.VERBOSE = 0;
	actx.min_rlen = 0;
	actx.min_olen = 1000;
	actx.SPUR_LEN = DEF_SPUR_LEN;
	actx.TIP_LEN = DEF_TIP_LEN;
	actx.nFuzzBases = DEF_ARG_F;
	actx.contByReads_CommonReadFraction = DEF_ARG_CRF;
	actx.contByReads_CoveredLenPerc = DEF_ARG_CL;
	actx.contByContigs_CoveredLenPerc = DEF_ARG_CL;

	actx.minPrimContigLen = DEF_ARG_L;
	actx.minPrimContigReads = DEF_ARG_N;
	actx.maxPrimContigRepeatPerc = DEF_ARG_P;

	int c;
	while ((c = getopt(argc, argv, "vx:l:c:d:t:r:C:F:s:e:D:o:R:f:L:N:P:")) != -1)
	{
		switch (c)
		{
		case 'o':
			filteredContigOvlsName = optarg;
			break;
		case 'v':
			actx.VERBOSE++;
			break;
		case 'L':
			actx.minPrimContigLen = atoi(optarg);
			break;
		case 'N':
			actx.minPrimContigReads = atoi(optarg);
			break;
		case 'P':
			actx.maxPrimContigRepeatPerc = atoi(optarg);
			break;
		case 'f':
			actx.nFuzzBases = atoi(optarg);
			break;
		case 'x':
			actx.min_rlen = atoi(optarg);
			break;
		case 'l':
			actx.min_olen = atoi(optarg);
			break;
		case 'e':
			actx.TIP_LEN = atoi(optarg);
			break;
		case 's':
			actx.SPUR_LEN = atoi(optarg);
			break;
		case 'c':
		{
			actx.exp_cov = atoi(optarg);
			if (actx.exp_cov <= 0)
			{
				fprintf(stderr, "[ERROR] - analyzeContigs: expected coverage must be positive (>= 1)\n");
				usage();
				exit(1);
			}
		}
			break;
		case 'd':
			actx.outDir = optarg;
			break;
		case 't':
			patchedReadTrimTrack = optarg;
			break;
		case 'r':
			corContigRepeatsTrack = optarg;
			break;
		case 'R':
			patchedReadRepeatTrack = optarg;
			break;
		case 'C':
		{
			actx.corContigDBName = optarg;
			if (optind < argc && *argv[optind] != '-')
			{
				actx.corContigLASName = argv[optind];
				optind++;
			}
			else
			{
				fprintf(stderr, "\n-C option require TWO arguments <contigDB> <ContigOverlaps>\n\n");
				usage();
				exit(1);
			}
			break;
		}
		case 'F':
		{
			actx.patchedReadDBName = optarg;
			if (optind < argc && *argv[optind] != '-')
			{
				actx.patchedReadLASName = argv[optind];
				optind++;
			}
			else
			{
				fprintf(stderr, "\n-F option require TWO arguments <fixedReadDB> <fixedReadOverlaps>\n\n");
				usage();
				exit(1);
			}
			break;
		}
		case 'D':
			actx.corReadDBName = optarg;
			break;
		default:
			fprintf(stderr, "Unknown option: %s\n", argv[optind - 1]);
			usage();
			exit(1);
		}
	}

	if (actx.corReadDBName == NULL)
	{
		fprintf(stderr, "[ERROR] - corrected read database is required!\n");
		exit(1);
	}

	if (actx.corContigDBName == NULL && actx.corContigLASName == NULL)
	{
		fprintf(stderr, "[ERROR] - At least a contig dabase and the corresponding overlaps are required to analyze contigs!\n");
		exit(1);
	}

	if ((correctedContigLAS = fopen(actx.corContigLASName, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.corContigLASName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Open database: %s\n", actx.corReadDBName);

	if (Open_DB(actx.corReadDBName, &correctedReadDB))
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.corReadDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Open database: %s\n", actx.corContigDBName);

	if (Open_DB(actx.corContigDBName, &correctedContigDB))
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.corContigDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Load track %s for db %s\n", corContigRepeatsTrack, actx.corContigDBName);
	if ((actx.corContigRepeats_track = track_load(&correctedContigDB, corContigRepeatsTrack)) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'. Create a repeat track with LArepeat first!\n", corContigRepeatsTrack,
				actx.corContigDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Load track %s for db %s\n", "creads", actx.corReadDBName);
	if ((actx.corContigCorrectReadsAndPos_track = track_load(&correctedContigDB, "creads")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "creads", actx.corContigDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Load track %s for db %s\n", "rreads", actx.corContigDBName);
	if ((actx.corContigRawReads_track = track_load(&correctedContigDB, "rreads")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "rreads", actx.corContigDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Load track %s for db %s\n", "preads", actx.corContigDBName);
	if ((actx.corContigPatchReadsAndPos_track = track_load(&correctedContigDB, "preads")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "preads", actx.corContigDBName);
		exit(1);
	}

	if ((patchedReadLAS = fopen(actx.patchedReadLASName, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.patchedReadLASName);
		exit(1);
	}

///todo make this optional
	if (actx.VERBOSE)
		printf("Open database: %s\n", actx.patchedReadDBName);

	if (Open_DB(actx.patchedReadDBName, &patchedReadDB))
	{
		fprintf(stderr, "could not open '%s'\n", actx.patchedReadDBName);
		exit(1);
	}
	if (patchedReadTrimTrack != NULL)
	{
		if (actx.VERBOSE)
			printf("Load track %s for db %s\n", patchedReadTrimTrack, actx.patchedReadDBName);
		if ((actx.patchedReadTrim_track = track_load(&patchedReadDB, patchedReadTrimTrack)) == NULL)
		{
			fprintf(stderr, "could not open track '%s' for database '%s'\n", patchedReadTrimTrack, actx.patchedReadDBName);
			exit(1);
		}
	}
	else if ((actx.patchedReadTrim_track = track_load(&patchedReadDB,
	TRACK_TRIM)) == NULL)
	{
		fprintf(stderr, "could not open track '%s' for database '%s'\n",
		TRACK_TRIM, actx.patchedReadDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Load track %s for db %s\n", patchedReadRepeatTrack, actx.patchedReadDBName);
	if ((actx.patchedReadRepeat_track = track_load(&patchedReadDB, patchedReadRepeatTrack)) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'!\n", patchedReadRepeatTrack, actx.patchedReadDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Load track %s for db %s\n", "source", actx.patchedReadDBName);
	if ((actx.patchedReadSource_track = track_load(&patchedReadDB, "source")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "source", actx.patchedReadDBName);
		exit(1);
	}

// init

	actx.corContigDB = &correctedContigDB;
	actx.patchedReadDB = &patchedReadDB;
	actx.corReadDB = &correctedReadDB;

	char *cwd = malloc(10);
	if (actx.outDir == NULL)
	{
		sprintf(cwd, ".");
		actx.outDir = cwd;
	}

	{
		if (createOutDir(actx.outDir))
			exit(1);

		char *out = malloc(strlen(actx.outDir) + 30);
		sprintf(out, "%s/raw", actx.outDir);
		if (createOutDir(out))
			exit(1);

		sprintf(out, "%s/split", actx.outDir);
		if (createOutDir(out))
			exit(1);

		free(out);
	}

	if (filteredContigOvlsName)
	{
		char *out = malloc(strlen(actx.outDir) + strlen(filteredContigOvlsName) + 10);
		sprintf(out, "%s/%s", actx.outDir, filteredContigOvlsName);

		if ((filteredContigOvlsFile = fopen(out, "w")) == NULL)
		{
			fprintf(stderr, "could not open %s\n", out);
			exit(1);
		}
		free(out);
	}

// contig pass context
	{
		contig_pctx = pass_init(correctedContigLAS, filteredContigOvlsFile);

		contig_pctx->split_b = 0;
		contig_pctx->load_trace = 0;
		contig_pctx->unpack_trace = 0;
		contig_pctx->data = &actx;

		pre(contig_pctx, &actx);
	}
// fixed read pass context
	{
		patched_pctx = pass_init(patchedReadLAS, NULL);
		patched_pctx->split_b = 0;
		patched_pctx->load_trace = 1;
		patched_pctx->unpack_trace = 1;
		patched_pctx->data = &actx;

		pre(patched_pctx, &actx);
	}
// initialize contigs
	printf("START    ---   STEP0a: initialize AnalyzeContext - START\n");
	initAnalyzeContext(&actx);
	printf("START    ---   STEP0a: initialize AnalyzeContext - DONE\n");
// todo for now do all steps: input format should be a MARVEL assembly run
// analyze overlaps of fixed reads for each contig ( coverage, containments, duplicate reads used in different contigs )
// STEP1
	if (0)
	{
		printf("START    ---   STEP1A: -- use patched read overlaps (from touring) and assign those to corresponding contigs\n");
		pass(patched_pctx, processReadOverlapsAndMapThemToContigs);
		printf("DONE     ---   STEP1A: -- use patched read overlaps (from touring) and assign those to corresponding contigs\n");

		printf("START    ---   STEP1B: -- do coverage analysis for contig-read-overlaps, set absolute contig positions\n");
		analyzeContigCoverageOfMappedReads(&actx);
		printf("DONE     ---   STEP1B: -- do coverage analysis for contig-read-overlaps, set absolute contig positions\n");

		printf("START    ---   STEP1C: -- preclassify Contig-Contig relationships by doing a mapped-read intersection strategy\n");
		preclassifyContigsByBReadsAndPath(&actx);
		printf("DONE     ---   STEP1C: -- preclassify Contig-Contig relationships by doing a mapped-read intersection strategy\n");
	}

// analyze overlaps between contigs
// STEP2
	if (1)
	{
		printf("START    ---   STEP2a: analyze contig vs contig overlaps\n");
		pass(contig_pctx, analyzeContigVsContigOverlaps);
		printf("DONE     ---   STEP2a: analyze contig vs contig overlaps\n");

		printf("START    ---   STEP2b: preclassify contig by contig overlaps\n");
		preClassifyContigsByContigOverlaps(&actx);
		printf("DONE     ---   STEP2b: preclassify contig by contig overlaps\n");
	}

// final classification, based on STEP1 and STEP2
// STEP3
	if (1)
	{
		printf("START    ---   STEP3: refine contig classification (based on STEP1 and STEP2)");
		classify(&actx);
		printf("DONE     ---   STEP3: refine contig classification (based on STEP1 and STEP2)");
	}
// check reads that occur multiple times in contigs, and set proper split positions if required
// STEP4
	{
		printf("STEP4: analyze reads that occur multiple times in contigs, and set proper split positions if required\n");
		// analyze number of shared reads between contigs, trim them back if possible /// split contigs ??
		finalContigValidation(&actx);
	}
// create output: contigs, splitted contigs, stat and bed files with i.e. coverage histogram, repeat tracks, split coordinates, ...
// STEP5
	{
		printf("STEP5: create output\n");
		rawClassification(&actx);
	}

	// clean up - put everything in a method
	pass_free(contig_pctx);
	pass_free(patched_pctx);

	fflush(stdout);
	fflush(stderr);

	Close_DB(&patchedReadDB);
	Close_DB(&correctedReadDB);
	Close_DB(&correctedContigDB);
	free(actx.readSeq - 1);

	int i;
	if (actx.contigFileNamesAndOffsets)
	{
		for (i = 0; i < actx.contigFileNamesAndOffsets->numFileNames; i++)
		{
			free(actx.contigFileNamesAndOffsets->fastaNames[i]);
			free(actx.contigFileNamesAndOffsets->fileNames[i]);
		}
		free(actx.contigFileNamesAndOffsets);
	}

	if (actx.rawReadFileNamesAndOffsets)
	{
		for (i = 0; i < actx.rawReadFileNamesAndOffsets->numFileNames; i++)
		{
			free(actx.rawReadFileNamesAndOffsets->fastaNames[i]);
			free(actx.rawReadFileNamesAndOffsets->fileNames[i]);
		}
		free(actx.rawReadFileNamesAndOffsets);
	}

	for (i = 0; i < actx.maxVReadMask; i++)
	{
		free(actx.vreadMask[i]);
	}
	free(actx.vreadMask);

	free(cwd);
	return 0;
}
