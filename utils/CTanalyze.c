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

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define DEF_REPEATS_TRACK "repeats"

#define	DEF_SPUR_LEN 50000
#define DEF_TIP_LEN 100000
#define DEF_ARG_F 10000
#define DEF_ARG_CRF 50
#define DEF_ARG_CL 50

#define DEF_ARG_L	100000 		// 	... minimum contig length, to be considered as a primary contig
#define DEF_ARG_N	5					//	... minimum number of contig reads,  to be considered as a primary contig
#define DEF_ARG_P	100				// 	... maximum repeat percentage, to be considered as a primary contig, not considered

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

void analyzeContigCreadIntersection(AnalyzeContext *actx)
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

	// int faId;

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

				if (contig_j->property.len < 100000)
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

					if ((contig_k->property.rflag & (REL_READ_IS_ALT)))
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

		if (ovls[k].path.aepos + actx->nFuzzBases < DB_READ_LEN(actx->corContigDB, ovls[j].aread))
		{
			properEndA = 0;
		}

		if (ovls[k].path.bepos + actx->nFuzzBases < DB_READ_LEN(actx->corContigDB, ovls[j].bread))
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

		Contig *conB = actx->contigs + ovls[j].bread;

		int validContainment = 0;
		int validBridge = 0;
		int validRepeatContainment = 0;

		if(properGapLen)
		{
			if(MAX(overlapBasesA, overlapBasesB) >= (int) (actx->contByContigs_CoveredLenPerc / 100.0 * MIN(conA->property.len, conB->property.len)))
			{
				validContainment = 1;
			}

			//      contigA         ----------------
			//      contigB	---------
			//			min 50Kb overlap, both overhangs min 100Kb
			else if(properBegA && !properEndA && properEndB && !properBegB && MAX(overlapBasesA, overlapBasesB) >= MIN(50000, 3*actx->nFuzzBases) && ovls[k].path.aepos + 100000 < conA->property.len && ovls[j].path.bbpos - 100000 > 0)
			{
				validBridge = 1;
			}
			//      contigA         ----------------
			//      								contigB			------------
			//			min 50Kb overlap, both overhangs min 100Kb
			else if(!properBegA && properEndA && !properEndB && properBegB && MAX(overlapBasesA, overlapBasesB) >= MIN(50000, 3*actx->nFuzzBases) && ovls[j].path.abpos - 100000 > 0 && ovls[k].path.bepos + 100000 < conB->property.len)
			{
				validBridge = 1;
			}
		}
		if (!validBridge && !validContainment)
		{
			if ((fail == 0 || fail == 4 || fail == 6) && MAX(overlapBasesA, overlapBasesB) >= (int) (0.4 * MIN(conA->property.len, conB->property.len)) &&
					(MAX(conA->property.repBasesFromReadLAS,conA->property.repBasesFromContigLAS)/100.0 * conA->property.len > actx->maxPrimContigRepeatPerc ||
					 MAX(conB->property.repBasesFromReadLAS,conB->property.repBasesFromContigLAS)/100.0 * conB->property.len > actx->maxPrimContigRepeatPerc))
			{
				validRepeatContainment = 1;
			}
		}

		if ((!fail && (validBridge || validContainment)) || validRepeatContainment)
		{

			if (conA->numContigRelations == conA->maxContigRelations)
			{
				conA->maxContigRelations = conA->maxContigRelations * 1.2 + 5;
				conA->contigRelations = (ContigRelation*) realloc(conA->contigRelations, sizeof(ContigRelation) * conA->maxContigRelations);
				bzero(conA->contigRelations + conA->numContigRelations, sizeof(ContigRelation) * (conA->maxContigRelations - conA->numContigRelations));
			}

			ContigRelation *crel = conA->contigRelations + conA->numContigRelations;

			if (validContainment)
			{
				crel->flag |= REL_CONTIG_IS_ALT;
				conA->property.rflag |= REL_CONTIG_IS_ALT;
				conB->property.rflag |= REL_CONTIG_HAS_ALT;
			}
			else if (validBridge)
			{
				crel->flag |= REL_CONTIG_BRIDGE;
				conA->property.rflag |= REL_CONTIG_BRIDGE;
				conB->property.rflag |= REL_CONTIG_BRIDGE;
			}
			else
			{
				crel->flag |= REL_CONTIG_IS_REPEAT_ALT;
				conA->property.rflag |= REL_CONTIG_IS_REPEAT_ALT;
				conB->property.rflag |= REL_CONTIG_HAS_ALT;
			}

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
			conA->numContigRelations++;
			printf("  ADD ContigRelation %d (l %d) vs %d (l %d) nOvls: %d alignedA %d / %d alignedB %d / %d (validBridge %d, validContainment %d)\n", conA->property.contigID, conA->property.len, conB->property.contigID, conB->property.len, crel->numPos, cumAaln, conA->property.len, cumBaln, conB->property.len, validBridge, validContainment);
		}
		else
		{
			Contig *conB = actx->contigs + ovls[j].bread;
			printf("  FAILED [%d,b%d,c%d,g%d] ContigRelation %d (l %d) vs %d (l %d) nOvls: %d pBegA %d , propEndA %d, propBegB %d, propEndB %d, ovlBasA %d, ovlBasB %d\n", fail, validBridge, validContainment,properGapLen, conA->property.contigID, conA->property.len,
					conB->property.contigID, conB->property.len, (k - j + 1), properBegA, properEndA, properBegB, properEndB, overlapBasesA, overlapBasesB);
		}
		k++;
		j = k;
	}

	printf("END -- Analyze overlaps for contig: %d numOvls: %d\n", ovls->aread, novl);

	return 1;
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

/// evaluate all three relationship information: TourRelation, ContigAlnRelation, ReadIntersectionRelation
void classify(AnalyzeContext *actx)
{
	int i, j;

	// i ... [0, ..., numContigs-1
	// j ... [0, ..., numContigs+1]
	// j=0. position: keep number of distinct ALT relationships
	// j=1. position: keep number of distinct CONTIG_BRIDGE relationships
	// j=2. - (n+1).position: contig relationships

	int ** relTable = (int **)malloc (sizeof(int*)*(actx->numContigs));
	for(i=0; i<actx->numContigs; i++)
	{
		relTable[i] = (int*) malloc (sizeof(int)*(actx->numContigs+2));
		bzero(relTable[i], sizeof(int)*(actx->numContigs+2));
	}

	for (i = actx->numContigs - 1;  i >= 0 ; i--)
	{
		Contig *conA = actx->contigs_sorted[i];

		// add containment touring relations to relation table
		for (j=0; j<conA->numTourRelations; j++)
		{
			TourRelation *tRel = conA->tourRelations + j;
			if (tRel->flag & REL_TOUR_IS_ALT)
			{
				if (relTable[conA->property.contigID][tRel->contigID1+2] == 0)
				{
					relTable[conA->property.contigID][0]++;
				}
				relTable[conA->property.contigID][tRel->contigID1+2] |= REL_TOUR_IS_ALT;
			}
		}

		// add containment contig relations to relation table
		for (j=0; j<conA->numContigRelations; j++)
		{
			ContigRelation *cRel = conA->contigRelations + j;
			if (cRel->flag & (REL_CONTIG_IS_ALT | REL_CONTIG_IS_REPEAT_ALT))
			{
				if (relTable[conA->property.contigID][cRel->corContigIdx+2] == 0)
				{
					relTable[conA->property.contigID][0]++;
				}
				if (cRel->flag & REL_CONTIG_IS_ALT)
					relTable[conA->property.contigID][cRel->corContigIdx+2] |= REL_CONTIG_IS_ALT;
				else
					relTable[conA->property.contigID][cRel->corContigIdx+2] |= REL_CONTIG_IS_REPEAT_ALT;
			}

			if (cRel->flag & (REL_CONTIG_BRIDGE))
			{
				if (relTable[conA->property.contigID][cRel->corContigIdx+2] == 0)
				{
					relTable[conA->property.contigID][1]++;
				}
				relTable[conA->property.contigID][cRel->corContigIdx+2] |= REL_CONTIG_BRIDGE;

				if (relTable[cRel->corContigIdx][conA->property.contigID+2] == 0)
				{
					relTable[cRel->corContigIdx][1]++;
				}
				relTable[cRel->corContigIdx][conA->property.contigID+2] |= REL_CONTIG_BRIDGE;
			}
		}

		// add containment read relations to relation table
		for (j=0; j<conA->numReadRelations; j++)
		{
			ReadRelation *rRel = conA->readRelations + j;
			if ((rRel->flag & REL_READ_IS_ALT))
			{
				if (relTable[conA->property.contigID][rRel->correspID+2] == 0)
				{
					relTable[conA->property.contigID][0]++;
				}
				relTable[conA->property.contigID][rRel->correspID+2] |= REL_READ_IS_ALT;
			}
			else if (rRel->flag & REL_READ_HAS_ALT)
			{
				if (relTable[rRel->correspID][conA->property.contigID+2] == 0)
				{
					relTable[rRel->correspID][0]++;
				}
				relTable[rRel->correspID][conA->property.contigID+2] |= REL_READ_IS_ALT;
			}

		}
	}

	// classify ALT, PRIM, and CRAP relationships
	for(i=0; i<actx->numContigs; i++)
	{
		Contig *conA = actx->contigs_sorted[i];

		int validRepeatPerc = (MAX(conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS)*100.0/conA->property.len < actx->maxPrimContigRepeatPerc) ? 1 : 0;
		int validMinCreads = (conA->numcReads > actx->minPrimContigReads) ? 1 : 0;
		int validMinLen = (conA->property.len > actx->minPrimContigLen) ? 1 : 0;

		if(relTable[conA->property.contigID][0] == 0)
		{
				if(validMinCreads && validRepeatPerc && validMinLen)
				{
					printf("CLASSIFY: PRIM %d l%d r(%d %d) v(%d, %d ,%d)\n",conA->property.contigID, conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);
					conA->property.cflag |= (CLASS_CONTIG_CLASSIFIED | CLASS_CONTIG_PRIMARY);
				}
				else
				{
					printf("CLASSIFY: CRAP %d l%d r(%d %d) v(%d, %d ,%d)\n",conA->property.contigID, conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);
					conA->property.cflag |= (CLASS_CONTIG_CLASSIFIED | CLASS_CONTIG_DISCARD);
					if(!validMinCreads)
						conA->property.cflag |= CLASS_CONTIG_DISCARD_CREADS;
					if(!validMinLen)
						conA->property.cflag |= CLASS_CONTIG_DISCARD_LEN;
					if(!validRepeatPerc)
						conA->property.cflag |= CLASS_CONTIG_DISCARD_REPEAT;
				}
		}
		else if(relTable[conA->property.contigID][0] == 1) // exactly one to one relationship
		{
				for(j=2; j<=actx->numContigs+1; j++)
				{
					if(relTable[conA->property.contigID][j]!=0)
						break;
				}
				assert(j<=actx->numContigs);
				printf("CLASSIFY: %d ALT{ of %d} l%d r(%d %d) v(%d, %d ,%d): ",conA->property.contigID, j-2, conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);
				if(relTable[conA->property.contigID][j] & REL_CONTIG_IS_ALT)
					printf(" REL_CONTIG_IS_ALT");
				if(relTable[conA->property.contigID][j] & REL_CONTIG_IS_REPEAT_ALT)
					printf(" REL_CONTIG_IS_REPEAT_ALT");
				if(relTable[conA->property.contigID][j] & REL_READ_IS_ALT)
					printf(" REL_READ_IS_ALT");
				if(relTable[conA->property.contigID][j] & REL_TOUR_IS_ALT)
					printf(" REL_TOUR_IS_ALT");
				printf("\n");

				conA->property.cflag |= (CLASS_CONTIG_ALT);
		}
	}

	// resolve MULTI CONTAINED relationships
	for(i=0; i<actx->numContigs; i++)
	{
		Contig *conA = actx->contigs_sorted[i];

		int validRepeatPerc = (MAX(conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS)*100.0/conA->property.len < actx->maxPrimContigRepeatPerc) ? 1 : 0;
		int validMinCreads = (conA->numcReads > actx->minPrimContigReads) ? 1 : 0;
		int validMinLen = (conA->property.len > actx->minPrimContigLen) ? 1 : 0;

		if(relTable[conA->property.contigID][0] > 1)
		{
			printf("CLASSIFY: %d MULTI l%d r(%d %d) v(%d, %d ,%d):\n",conA->property.contigID, conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);

			if (relTable[conA->property.contigID][0] > 2)
			{
				conA->property.cflag |= (CLASS_CONTIG_CLASSIFIED | CLASS_CONTIG_DISCARD);
				printf("CLASSIFY: %d MULTIRSOLVED CRAP nrRel%d l%d r(%d %d) v(%d, %d ,%d):\n",conA->property.contigID, relTable[conA->property.contigID][0], conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);
			}
			else
			{
				int c=1;
				int numPrim = 0;
				int numAlt  = 0;
				for(j=2; j<=actx->numContigs+1; j++)
				{
					if(relTable[conA->property.contigID][j]!=0)
					{
						printf("  %d. REL with c%d: ", c++, j-2);
						if(relTable[conA->property.contigID][j] & REL_CONTIG_IS_ALT)
							printf(" REL_CONTIG_IS_ALT");
						if(relTable[conA->property.contigID][j] & REL_CONTIG_IS_REPEAT_ALT)
							printf(" REL_CONTIG_IS_REPEAT_ALT");
						if(relTable[conA->property.contigID][j] & REL_READ_IS_ALT)
							printf(" REL_READ_IS_ALT");
						if(relTable[conA->property.contigID][j] & REL_TOUR_IS_ALT)
							printf(" REL_TOUR_IS_ALT");
						printf("\n");

						if(actx->contigs[j-2].property.cflag & CLASS_CONTIG_ALT)
							numAlt++;
						if(actx->contigs[j-2].property.cflag & CLASS_CONTIG_PRIMARY)
							numAlt++;
					}

				}
				assert(numAlt + numPrim == 2);

				if(numAlt)
				{
					conA->property.cflag |= (CLASS_CONTIG_CLASSIFIED | CLASS_CONTIG_DISCARD);
					printf("CLASSIFY: %d MULTIRSOLVED CRAP ContInCont l%d r(%d %d) v(%d, %d ,%d):\n",conA->property.contigID, conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);
				}
				else
				{
					conA->property.cflag |= (CLASS_CONTIG_CLASSIFIED | CLASS_CONTIG_ALT);
					printf("CLASSIFY: %d MULTIRSOLVED ALT l%d r(%d %d) v(%d, %d ,%d):\n",conA->property.contigID, conA->property.len, conA->property.repBasesFromContigLAS, conA->property.repBasesFromReadLAS, validRepeatPerc, validMinCreads, validMinLen);
				}
			}
		}
	}

	// check bridging contigs and trim back overhang
	for(i=0; i<actx->numContigs; i++)
	{
		Contig *conA = actx->contigs_sorted[i];

		if(relTable[conA->property.contigID][1] > 1)
		{
			for(j=2; j<=actx->numContigs+1; j++)
			{
				if(relTable[conA->property.contigID][j]!=0)
				{

					if(conA->property.contigID < j-2)
					{
						printf("found putative bridge: %d vs %d\n", conA->property.contigID, j-2);
					}
				}
			}
		}
	}
}

static void usage()
{
	fprintf(stderr,
			" [-v] [-csfeLNP <int>] [-d <dir>] [-rt <Track>] -C <contigDB> <ContigOverlaps> -F <fixedReadDB> <fixedReadOverlaps> -D <correctedReadDB> \n");
	fprintf(stderr, "options: -v         ... verbose\n");
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

	actx.VERBOSE = 0;
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
	while ((c = getopt(argc, argv, "vc:d:t:r:C:F:s:e:D:R:f:L:N:P:")) != -1)
	{
		switch (c)
		{
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

// contig pass context
	{
		contig_pctx = pass_init(correctedContigLAS, NULL);

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
	if (1)
	{
		printf("START    ---   STEP1A: -- use patched read overlaps (from touring) and assign those to corresponding contigs\n");
		pass(patched_pctx, processReadOverlapsAndMapThemToContigs);
		printf("DONE     ---   STEP1A: -- use patched read overlaps (from touring) and assign those to corresponding contigs\n");

		printf("START    ---   STEP1B: -- do coverage analysis for contig-read-overlaps, set absolute contig positions\n");
		analyzeContigCoverageOfMappedReads(&actx);
		printf("DONE     ---   STEP1B: -- do coverage analysis for contig-read-overlaps, set absolute contig positions\n");

		printf("START    ---   STEP1C: -- analyzeContigCreadIntersection by applying a read intersection strategy\n");
		analyzeContigCreadIntersection(&actx);
		printf("DONE     ---   STEP1C: -- analyzeContigCreadIntersection by applying a read intersection strategy\n");
	}

// analyze overlaps between contigs
// STEP2
	if (1)
	{
		printf("START    ---   STEP2a: analyze contig vs contig overlaps\n");
		pass(contig_pctx, analyzeContigVsContigOverlaps);
		printf("DONE     ---   STEP2a: analyze contig vs contig overlaps\n");
	}

// final classification, based on STEP1 and STEP2
// STEP3
	if (1)
	{
		printf("START    ---   STEP3: refine contig classification (based on STEP1 and STEP2)");
		classify(&actx);
		printf("DONE     ---   STEP3: refine contig classification (based on STEP1 and STEP2)");
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
