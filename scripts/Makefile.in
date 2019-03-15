
include ../Makefile.settings

ALL =
SCRIPTS = rid2bid.py DBprepare.py paths2rids.py marvelStats.sh n50.py splitDiploidAssembly.py \
		  marvel_slurm.sh createAndSubmitMarvelSlurmJobs.sh createCommandPlan.sh createSubmitScripts.sh \
		  createRepmaskPlans.sh createReadPatchingPlans.sh \
          createRepmaskPlans2.sh createScrubbingPlans.sh createFilteringPlans.sh createTouringPlans.sh \
          createCorrectionPlans.sh createContigAnalyzePlans.sh \
          createPacBioArrowPlans.sh \
          slurmStats.sh assemblyStats.sh exampleConfig.sh \
          createDAScoverPlans.sh createPurgeHaplotigPlans.sh createFreeBayesPolishingPlans.sh \
          createHiCPlans.sh get_stats.pl filter_five_end.pl two_read_bam_combiner.pl \
          createMitoAssemblyPlans.sh createPhasePlans.sh trimLowerCaseTips.py \
          create10XPlans.sh splitReads.py createBionanoPlans.sh \
	  	  generate_site_positions.py \
          createQCandStatsPlans.sh mashPlot.R

all: $(ALL)

install: all
	$(INSTALL_PROGRAM) -m 0755 $(SCRIPTS) $(install_scripts)

clean:

