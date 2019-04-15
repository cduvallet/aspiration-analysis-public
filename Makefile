## Makefile to reproduce data cleaning, analyses, and figures
## for aspiration paper
# Author: Claire Duvallet

# Define some directory shortcuts
RESULTSDIR = data/analysis
SRCDIR = src/analysis


all: DATA ANALYSIS FIGURES TABLES SUPP_FILES #TREE

# These are the analyses that take a long time to do
long_stuff: CLASSIFIERS EXCHANGE_LONG

#############################
#       DATA                #
#############################

# Inputs
meta_dir=data/raw/metadata
raw_otu=data/raw/rosen_mincount10_maxee2_trim200_results/RDP/rosen_mincount10_maxee2_trim200.otu_table.99.denovo.rdp_assigned
batch2014=data/raw/rosen_mincount10_maxee2_trim200_results/2014_batch.txt
batch2016=data/raw/rosen_mincount10_maxee2_trim200_results/2016_batch.txt

# Outputs
wrangled_metadata=data/clean/all_metadata.txt
clean_otu_counts=data/clean/rosen.otu_table.counts.clean
clean_otu=data/clean/rosen.otu_table.rel_abun.clean
clean_metadata=data/clean/rosen.metadata.clean

DATA: $(wrangled_metadata) $(clean_otu) $(clean_metadata)

# Wrangle the metadata in various files into a consolidated, sample-wise
# metadata table.
# **Note** that this rule also depends on the raw files in data/raw/metadata/
$(wrangled_metadata): src/data/metadata_wrangling.py $(raw_otu) $(batch2014) $(batch2016)
	python src/data/metadata_wrangling.py $(meta_dir) $(raw_otu) $(batch2014) $(batch2016) $(wrangled_metadata)

$(clean_otu): src/data/clean_otu_table.py $(raw_otu) $(wrangled_metadata)
	python src/data/clean_otu_table.py $(raw_otu) $(wrangled_metadata) $(clean_otu_counts) $(clean_otu) $(clean_metadata) --samplereads 2000

# The script that makes the clean OTU table also makes the clean metadata
# this rule just checks for the case in which I deleted clean_metadata
# without modifying the clean OTU table or any of its prerequisites
$(clean_metadata): $(clean_otu)
	  @if test -f $@; then :; else \
	    rm -f $(clean_otu); \
	    $(MAKE) $(AM_MAKEFLAGS) $(clean_otu); \
	  fi

#############################
#       TREE                #
#############################

## Make the phylogenetic tree
# Inputs
raw_seqs = data/raw/rosen_mincount10_maxee2_trim200_results/rosen_mincount10_maxee2_trim200.otu_seqs.99.fasta

# Outputs
aligned = data/tree/rosen.otu_seqs.aligned
align_log = data/tree/rosen.otu_seqs.aligned.log
tree = data/tree/rosen.otu_seqs.newick
tree_log = data/tree/rosen.otu_seqs.newick.log
TREE: $(tree)

## Alignment
# Note: PyNAST aligns to a template (e.g. aligned GG), MUSCLE doesn't need
# template.
$(aligned): $(raw_seqs)
	muscle -in $< -out $@ -log $(align_log)

$(tree): $(aligned)
	FastTree -log $(tree_log) -nt $< > $@

#############################
#       ANALYSIS            #
#############################

#univariate=data/analysis/univariate_qvalues.txt

## Beta diversity
# JSD (all files are made in jsd.py):
#    jsd.txt has metadata in individual columns.
#    jsd.wide.txt is in symmetric distance matrix format
#    jsd.tidy_reflux_bile.txt is jsd.txt with the reflux and bile
#    columns tidyfied (var='bile_type' or 'reflux_type', val='bile_value' or
#    'reflux_value')
jsd = $(RESULTSDIR)/jsd.txt
#widejsd = $(RESULTSDIR)/jsd.wide.txt
#refluxjsd = $(RESULTSDIR)/jsd.tidy_reflux_bile.txt

braycurtis = $(RESULTSDIR)/braycurtis.txt

# Unifrac
#weighted_unifrac = $(RESULTSDIR)/weighted_unifrac.txt
#unweighted_unifrac = $(RESULTSDIR)/unweighted_unifrac.txt
beta_div: $(jsd) #$(weighted_unifrac) $(unweighted_unifrac)

# Within vs. between patient calculations
within_vs_between_jsd = $(RESULTSDIR)/within_vs_between_beta.txt
within_vs_between_bc = $(RESULTSDIR)/within_vs_between_beta.braycurtis.txt
within_vs_between_analysis: $(within_vs_between_jsd)

# Alpha diversity
alpha := $(RESULTSDIR)/alpha_diversity.txt

# Differential abundance
qvalues := $(RESULTSDIR)/qvalues.kruskal_wallis.txt
#qvalues_ppi := $(RESULTSDIR)/qvalues.kruskal_wallis.ppi.txt
#qvalues_steroid := $(RESULTSDIR)/qvalues.kruskal_wallis.inhaled_steroids.txt
#qvalues_pneum := $(RESULTSDIR)/qvalues.kruskal_wallis.pneumonia.txt
#qvalues_abx := $(RESULTSDIR)/qvalues.kruskal_wallis.abx.txt

reviewer_analysis: $(alpha) $(qvalues) $(braycurtis)

## Exchanged OTUs
exchange_w_partial = $(RESULTSDIR)/exchange.with_partial_corrs.txt
null_exchange_w_partial = $(RESULTSDIR)/exchange.with_partial_corrs.null_20reps.txt
exchanged_otus = $(RESULTSDIR)/exchange.labeled_otus.txt
exchange_prevalence = $(RESULTSDIR)/prevalence.partial_corrs.nthresh10-qthresh01-rthresh0.txt
exchange: $(exchange_w_partial) $(exchanged_otus) $(exchange_prevalence)

EXCHANGE_LONG: $(null_exchange_w_partial)

## Classifiers
# Based on entire communities
rf_summaries = $(RESULTSDIR)/rf_results.summary_stats.txt
#rf_rocs = $(RESULTSDIR)/rf_results.rocs.txt
rf_predictions = $(RESULTSDIR)/rf_results.predictions.txt
# Based on exchanged OTUs only
rf_summaries_exchange = $(RESULTSDIR)/rf_results.exchanged.summaries.txt
#rf_rocs_exchange = $(RESULTSDIR)/rf_results.exchanged.rocs.txt
rf_predictions_exchange = $(RESULTSDIR)/rf_results.exchanged.predictions.txt

CLASSIFIERS: $(rf_summaries) $(rf_summaries_exchange)

ANALYSIS: beta_div exchange within_vs_between_analysis reviewer_analysis

#######################
#   ANALYSIS RECIPES  #
#######################

$(jsd): $(SRCDIR)/jsd.py $(clean_otu) $(clean_metadata)
	python src/analysis/jsd.py $(clean_otu) $(clean_metadata) $(jsd)

#$(widejsd): $(jsd)
#	@if test -f $@; then :; else \
#	  rm -f $<; \
#	  make $<; \
#	fi

$(braycurtis): $(SRCDIR)/bray_curtis.py $(clean_otu_counts) $(clean_metadata)
	python $< $(clean_otu_counts) $(clean_metadata) $@

#$(weighted_unifrac): $(SRCDIR)/unifrac.py $(clean_otu_counts) $(clean_metadata) $(tree)
#	python $< $(clean_otu_counts) $(tree) $(clean_metadata) $@ --weighted

#$(unweighted_unifrac): $(SRCDIR)/unifrac.py $(clean_otu_counts) $(clean_metadata) $(tree)
#	python $< $(clean_otu_counts) $(tree) $(clean_metadata) $@

# Calculate correlations between sites (AKA exchanged-ness)
# Also run the script that tracks patients used in this analysis
$(exchange_w_partial): $(SRCDIR)/exchange_tidy.py $(clean_otu) $(clean_metadata)
	python $< $(clean_otu) $(clean_metadata) $@
	python $(SRCDIR)/exchange_track_samples.py $(clean_otu) $(clean_metadata) data/patients/exchange_calculation --npatients 10

$(null_exchange_w_partial): $(SRCDIR)/exchange_null.py $(clean_otu) $(clean_metadata)
	python $< $(clean_otu) $(clean_metadata) $@ --nthresh 10 --nshuffle 20

# Label the exchanged OTUs in a separate file
$(exchanged_otus): $(SRCDIR)/label_exchanged_otus.py $(exchange_w_partial)
	python $< $(exchange_w_partial) $@ --nthresh 10 --qthresh 0.1

# Calculate the prevalence of the exchanged OTUs. Note that thresholds should
# match the ones used in the recipe to label them
# Also track the samples used in the prevalence calculations
$(exchange_prevalence): $(SRCDIR)/prevalence_exchange.py $(exchange_w_partial) $(clean_otu) $(clean_metadata)
	python $< $(exchange_w_partial) $(clean_otu) $(clean_metadata) $@ --nthresh 10 --qthresh 0.1
	python $(SRCDIR)/prevalence_exchange_track_patients.py $(exchange_w_partial) $(clean_otu) $(clean_metadata) data/patients/prevalence_calculation --nthresh 10 --qthresh 0.1

# This also makes all of the other RF files
$(rf_summaries): $(SRCDIR)/make_loo_classifiers.py $(clean_otu) $(clean_metadata)
	python $< $(clean_otu) $(clean_metadata) $(rf_summaries) $(rf_predictions)

# RF exchanged
# Note: this python script is *very* similar to the one above, I could
# almost certainly combine them more intelligently. But for now, let's
# just stick with the redundancy, oh well.
$(rf_summaries_exchange): $(SRCDIR)/make_loo_classifiers_exchanged_OTUs.py $(clean_otu) $(clean_metadata) $(exchange_prevalence)
	python $< $(clean_otu) $(clean_metadata) $(exchange_prevalence) $(rf_summaries_exchange) $(rf_predictions_exchange)

# Alpha diversity
$(alpha): $(SRCDIR)/alpha_diversity.py $(clean_otu_counts) $(clean_metadata)
	python $< $(clean_otu_counts) $(clean_metadata) $@

# Within vs. between patient JSD
$(within_vs_between_jsd): $(SRCDIR)/within_beta_vs_between_beta.py
	python $< $(jsd) $(clean_metadata) $@

# Within vs. between patient bray curtis
$(within_vs_between_bc): $(SRCDIR)/within_beta_vs_between_beta.py
	python $< $(braycurtis) $(clean_metadata) $@

# Differential abundance
$(qvalues): $(SRCDIR)/get_qvalues.py $(clean_otu) $(clean_meatadata)
	python $< $(clean_otu) $(clean_metadata) $@ --method kruskal-wallis

# $(qvalues_ppi): $(SRCDIR)/get_qvalues.py $(clean_otu) data/clean/rosen.clinical_metadata.clean
# 	python $< $(clean_otu) data/clean/rosen.clinical_metadata.clean $@ --method kruskal-wallis --metacol ppi_all --caselabel 1 --ctrllabel 0
#
# $(qvalues_steroid): $(SRCDIR)/get_qvalues.py $(clean_otu) data/clean/rosen.clinical_metadata.clean
# 	python $< $(clean_otu) data/clean/rosen.clinical_metadata.clean $@ --method kruskal-wallis --metacol inhaled_steroids_all --caselabel 1 --ctrllabel 0
#
# $(qvalues_pneum): $(SRCDIR)/get_qvalues.py $(clean_otu) data/clean/rosen.clinical_metadata.clean
# 	python $< $(clean_otu) data/clean/rosen.clinical_metadata.clean $@ --method kruskal-wallis --metacol pneum_all --caselabel 1 --ctrllabel 0
#
# $(qvalues_abx): $(SRCDIR)/get_qvalues.py $(clean_otu) data/clean/rosen.clinical_metadata.clean
# 	python $< $(clean_otu) data/clean/rosen.clinical_metadata.clean $@ --method kruskal-wallis --metacol abx_all --caselabel 1 --ctrllabel 0


#qvalues_steroid := $(RESULTSDIR)/qvalues.kruskal_wallis.inhaled_steroid.txt

#############################
#       FIGURES             #
#############################

# Figures are going to be in the notebooks, which I will commit to the repo directly.
# Might want to add something about making the text files with the patients I use in each figure?

# Original directory
EXPDIR = src/exploration
# Final directory, where only final notebooks will live
FIGSRC = src/figures
# Where to write figures
FIGDIR = final/figures
# Where to write text files with patients
PTDIR = final/patients

#############################
####### MAIN FIGURES ########
#############################
# Note: all figure names and patient list file names are hardcoded in their
# respective notebooks

# Figure 1 - community overview
overview_fig = $(FIGDIR)/figure1.overview_barplots.png
overview_pts = $(PTDIR)/figure1.overview_barplots.subjects.txt
# Notebooks
overview_nb = $(EXPDIR)/2018-02-20.community_overview_ordered_by_aspiration_final.ipynb
overview_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(overview_nb))
OVERVIEW_FIG: $(overview_fig)

# Figure 2 - across-patient beta diversity and PCoA plot
pca_fig = $(FIGDIR)/figure2.pcoa.png
acrosspt_fig = $(FIGDIR)/figure2.across_patient_beta_div_same_site.png
across_pts = $(PTDIR)/figure2.across_patient_beta_div.subjects.txt
# Notebooks
across_nb = $(EXPDIR)/2018-01-02.pca_and_between_patient_div_final.ipynb
across_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(across_nb))
ACROSS_FIG: $(pca_fig) $(acrosspt_fig)

# Figure 3 - within-patient beta diversity
within_fig = $(FIGDIR)/figure3.within_patient_beta_div.png
sitecls_fig = $(FIGDIR)/figure3.site_classifiers.png
#within_pts = $(PTDIR)/figure3.within_comparisons.subjects.txt
# Notebooks
within_nb = $(EXPDIR)/2017-12-26.within_patient_beta_div_figure_final.ipynb
within_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(within_nb))
sitecls_nb = $(EXPDIR)/2018-01-26.classify_sites_final.ipynb
sitecls_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(sitecls_nb))
WITHIN_FIG: $(within_fig) $(sitecls_fig)

# Figure 4 - within-patient vs between-patient JSD
within_vs_between_fig = $(FIGDIR)/figure.within_pt_vs_between_pt.png
# Notebook
within_vs_between_nb = $(EXPDIR)/2019-01-16.within_vs_between_figure_final.ipynb
within_vs_between_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(within_vs_between_nb))
WITHIN_VS_BETWEEN_FIG: $(within_vs_between_fig)

# Figure 4 - beta-diversity and exchange in aspirators vs. non-aspirators
aspbeta_fig = $(FIGDIR)/figure4.asp_vs_nonasp_beta.png
aspex_fig = $(FIGDIR)/figure4.asp_vs_nonasp_exchange.png
#asp_pts = $(PTDIR)/figure4.asp_vs_nonasp.subjects.txt
# Notebooks
aspex_nb = $(EXPDIR)/2018-01-05.asp_vs_nonasp_exchange_fig_final.ipynb
aspex_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(aspex_nb))
aspbeta_nb = $(EXPDIR)/2018-01-03.asp_vs_nonasp_beta_fig_final.ipynb
aspbeta_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(aspbeta_nb))
ASP_FIG: $(aspex_fig) $(aspbeta_fig)

# Figure 5 - correlation between reflux and bal-gastric JSD
reflux_fig = $(FIGDIR)/figure5.reflux_correlations.png
reflux_pts = $(PTDIR)/figure5.reflux_correlation.subjects.txt
# Notebooks
reflux_nb = $(EXPDIR)/2018-01-23.reflux_correlations_final.ipynb
reflux_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(reflux_nb))
REFLUX_FIG: $(reflux_fig)

#############################
####### SUPP FIGURES ########
#############################

## No longer making these figures, since we switch to leave-one-out
## classification intead of repeating 100 times
# # Figure 7 and 8 - AUCs and fisher p values for exchanged OTUs classifiers
# exch_aucs_fig = $(FIGDIR)/figure7.auc_asp_classifiers.exchangedOTUs.png
# exch_p_fig = $(FIGDIR)/figure8.fisherp_asp_classifiers.exchangedOTUs.png
# #exch_cls_pts = $(PTDIR)/figure7.exchangedOTUs_classifiers.subjects.txt
# # Notebooks
# cls_nb = $(EXPDIR)/2018-01-19.aspiration_classifiers_results_final.ipynb
# cls_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(cls_nb))
#
# EXCH_CLS_FIG: $(exch_aucs_fig) $(exch_p_fig)
#
# # Figure 9 - AUCs and p values for full-abundance classifiers
# asp_cls_fig = $(FIGDIR)/figure9.auc_fisherp_asp_classifiers.png
# # Notebook: this one is also made by cls_nb
# ASP_CLS_FIG: $(asp_cls_fig)

# Alpha diversity
alpha_fig = $(FIGDIR)/suppfig.alpha_diversity.png
alpha_fig_nb = $(EXPDIR)/2019-01-24.alpha_diversity_reviewer_response_final.ipynb
alpha_fig_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(alpha_fig_nb))
ALPHA_FIG: $(alpha_fig)

## Bray-curtis figures

# Across-patient beta diversity and PCoA plot - bray-curtis
pca_fig_bc = $(FIGDIR)/suppfig.pcoa.braycurtis.png
acrosspt_fig_bc = $(FIGDIR)/suppfig.across_patient_beta_div_same_site.braycurtis.png
# Notebooks
across_nb_bc = $(EXPDIR)/2019-01-27.pca_and_between_patient_div_braycurtis.ipynb
across_final_nb_bc = $(subst $(EXPDIR),$(FIGSRC),$(across_nb_bc))
ACROSS_FIG_BRAYCURTIS: $(pca_fig_bc) $(acrosspt_fig_bc)

# Within-patient beta diversity - bray curtis
within_fig_bc = $(FIGDIR)/suppfig.within_patient_beta_div.braycurtis.png
# Notebooks
within_nb_bc = $(EXPDIR)/2019-01-27.within_patient_beta_div_figure_braycurtis.ipynb
within_final_nb_bc = $(subst $(EXPDIR),$(FIGSRC),$(within_nb_bc))
WITHIN_FIG_BRAYCURTIS: $(within_fig_bc)

# Within-patient vs between-patient Bray-Curtis
within_vs_between_fig_bc = $(FIGDIR)/suppfig.within_pt_vs_between_pt.braycurtis.png
# Notebook
within_vs_between_nb_bc = $(EXPDIR)/2019-01-30.within_vs_between_figure_final_braycurtis.ipynb
within_vs_between_final_nb_bc = $(subst $(EXPDIR),$(FIGSRC),$(within_vs_between_nb_bc))
WITHIN_VS_BETWEEN_FIG_BRAYCURTIS: $(within_vs_between_fig_bc)

# Bray curtis beta-diversity and exchange in aspirators vs. non-aspirators
aspbeta_fig_bc = $(FIGDIR)/suppfig.asp_vs_nonasp_beta.braycurtis.png
# Notebooks
aspbeta_nb_bc = $(EXPDIR)/2019-01-27.asp_vs_nonasp_beta_fig_braycurtis.ipynb
aspbeta_final_nb_bc = $(subst $(EXPDIR),$(FIGSRC),$(aspbeta_nb_bc))
ASP_FIG_BRAYCURTIS: $(aspbeta_fig_bc)

## PPI vs. lung-gastric
ppi_lunggastric = $(FIGDIR)/reviewer.ppi_vs_lung_gastric.png
ppi_lunggastric_nb = $(EXPDIR)/2019-02-03.gastric_lung_ppi_use_reviewer.ipynb
ppi_lunggastric_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(ppi_lunggastric_nb))
PPI_LUNGGASTRIC: $(ppi_lunggastric)

## Reflux vs lung-gastric JSD, colored by PPI
ppi_reflux = $(FIGDIR)/reviewer.reflux_correlations_ppi_status.png
# This figure is made by the same notebook as above, so don't include
# if in make

total_reads = $(FIGDIR)/suppfig.total_reads.png
total_reads_nb = $(EXPDIR)/2019-02-02.total_reads_reviewer_response.ipynb
total_reads_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(total_reads_nb))
TOTAL_READS: $(total_reads)

BRAYCURTIS_FIGS: ACROSS_FIG_BRAYCURTIS WITHIN_FIG_BRAYCURTIS WITHIN_VS_BETWEEN_FIG_BRAYCURTIS ASP_FIG_BRAYCURTIS

FIGURES: OVERVIEW_FIG ACROSS_FIG WITHIN_FIG ASP_FIG REFLUX_FIG WITHIN_VS_BETWEEN_FIG ALPHA_FIG BRAYCURTIS_FIGS PPI_LUNGGASTRIC TOTAL_READS

################################
#######  FIG RECIPES    ########
################################

# Figure 1 - community overview
$(overview_fig): $(overview_final_nb) $(clean_otu) $(clean_metadata)
	jupyter nbconvert --execute $< --ExecutePreprocessor.timeout=-1

$(overview_final_nb): $(overview_nb)
	cp $< $@

# Figure 2 - across-patient beta diversity and PCoA plot
# Note: this notebook actually uses the jsd.wide.txt file, but I only
# include the rule to make jsd.txt (which automatically makes
# jsd.wide.txt), hence why this depends on just $(jsd)
$(pca_fig): $(across_final_nb) $(clean_otu) $(clean_metadata) $(jsd)
	jupyter nbconvert --execute $< --ExecutePreprocessor.timeout=-1

$(acrosspt_fig): $(pca_fig)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

$(across_final_nb): $(across_nb)
	cp $< $@

# Figure 3 - within-patient beta diversity
$(within_fig): $(within_final_nb) $(jsd)
	jupyter nbconvert --execute $<

$(within_final_nb): $(within_nb)
	cp $< $@

$(sitecls_fig): $(sitecls_final_nb) $(clean_otu) $(clean_metadata)
	jupyter nbconvert --execute $< --ExecutePreprocessor.timeout=-1

$(sitecls_final_nb): $(sitecls_nb)
	cp $< $@

# Figure 4 - within-patient vs between-patient beta diversity
$(within_vs_between_fig): $(within_vs_between_final_nb) $(within_vs_between_jsd)
	jupyter nbconvert --execute $<

$(within_vs_between_final_nb): $(within_vs_between_nb)
	cp $< $@

# Figure 4 - beta-diversity and exchange in aspirators vs. non-aspirators
$(aspex_fig): $(aspex_final_nb) $(exchange_prevalence)
	jupyter nbconvert --execute $<

$(aspex_final_nb): $(aspex_nb)
	cp $< $@

$(aspbeta_fig): $(aspbeta_final_nb) $(jsd)
	jupyter nbconvert --execute $<

$(aspbeta_final_nb): $(aspbeta_nb)
	cp $< $@

# Figure 5 - correlation between reflux and bal-gastric JSD
$(reflux_fig): $(reflux_final_nb) $(jsd)
	jupyter nbconvert --execute $<

$(reflux_final_nb): $(reflux_nb)
	cp $< $@

# Figure 6 is the exchange schematic

# # Figure 7 and 8 - AUCs and fisher p values for exchanged OTUs classifiers
# # Also update the notebook with the tables of exchanged OTUs
# $(exch_aucs_fig): $(cls_final_nb) $(rf_summaries_exchange)
# 	jupyter nbconvert --execute $<
#
# $(exch_p_fig): $(exch_aucs_fig)
# 	@if test -f $@; then :; else \
# 	  rm -f $<; \
# 	  make $<; \
# 	fi
#
# $(cls_final_nb): $(cls_nb)
# 	cp $< $@
#
# # Figure 9 - AUCs and p values for full-abundance classifiers
# $(asp_cls_fig): $(cls_final_nb) $(rf_summaries)
# 	jupyter nbconvert --execute $<

# Alpha diversity
$(alpha_fig): $(alpha_fig_final_nb) $(alpha)
	jupyter nbconvert --execute $<

$(alpha_fig_final_nb): $(alpha_fig_nb)
	cp $< $@

# Bray curtis figures
$(pca_fig_bc): $(across_final_nb_bc) $(braycurtis)
	jupyter nbconvert --execute $< --ExecutePreprocessor.timeout=-1

$(acrosspt_fig_bc): $(pca_fig_bc)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

$(across_final_nb_bc): $(across_nb_bc)
	cp $< $@

$(within_fig_bc): $(within_final_nb_bc) $(braycurtis)
	jupyter nbconvert --execute $<

$(within_final_nb_bc): $(within_nb_bc)
	cp $< $@

$(within_vs_between_fig_bc): $(within_vs_between_final_nb_bc) $(within_vs_between_bc)
	jupyter nbconvert --execute $<

$(within_vs_between_final_nb_bc): $(within_vs_between_nb_bc)
	cp $< $@

$(aspbeta_fig_bc): $(aspbeta_final_nb_bc) $(braycurtis)
	jupyter nbconvert --execute $<

$(aspbeta_final_nb_bc): $(aspbeta_nb_bc)
	cp $< $@

# PPI vs. lung-gastric
# this also makes reflux vs. lung gastric, colored by PPI status
$(ppi_lunggastric): $(ppi_lunggastric_final_nb) $(jsd) $(clinical_metadata)
	jupyter nbconvert --execute $<

$(ppi_lunggastric_final_nb): $(ppi_lunggastric_nb)
	cp $< $@

## Total reads

$(total_reads): $(total_reads_final_nb) $(clean_otu) $(clean_meta)
	jupyter nbconvert --execute $<

$(total_reads_final_nb): $(total_reads_nb)
	cp $< $@

###########################
####      TABLES      #####
###########################

# Where to write code that makes tables
TABSRC = src/tables

# Note: allt hese tables get printed as latex string
# in their respective notebooks, so they don't actually
# produce any files. That's why the notebook is the target.

## Table 1 and S1, patient demographics
demographics_nb = $(EXPDIR)/2019-01-26.wrangle_new_metadata_table1_s1.ipynb
demographics_final_nb = $(subst $(EXPDIR),$(TABSRC),$(demographics_nb))

$(demographics_final_nb): $(demographics_nb) $(clean_meta) $(clean_otu) $(clinical_metadata)
	cp $< $@
	jupyter nbconvert --execute $@

# Table 2 and supp table 2: sample sites sequenced
sites_nb = $(EXPDIR)/2019-01-27.sample_info_table_2.ipynb
sites_final_nb = $(subst $(EXPDIR),$(TABSRC),$(sites_nb))

$(sites_final_nb): $(sites_nb) $(patient_list)
	cp $< $@
	jupyter nbconvert --execute $@

### Reviewer response update: we do this for all metadata now
## Notebook that quickly checks for PPI confounding
#ppi_nb = $(EXPDIR)/2018-02-14.ppi_confounding_check_final.ipynb
#ppi_final_nb = $(subst $(EXPDIR),$(TABSRC),$(ppi_nb))
#
#$(ppi_final_nb): $(ppi_nb)
#	cp $< $@
#	jupyter nbconvert --execute $@

# Within vs between comparison (former main table, that's replace
# by a figure now)
# This one is made in the same notebook as the figure
within_btw_table_final_nb = $(subst $(FIGSRC),$(TABSRC),$(within_vs_between_final_nb))

$(within_btw_table_final_nb): $(within_vs_between_final_nb)
	cp $< $@
	# Don't run it, because that re-makes the figure

# Tables of exchanged bugs
exchanged_tables = $(EXPDIR)/2018-01-25.exchanged_OTUs_tables_final.ipynb
exchanged_tables_final = $(subst $(EXPDIR),$(TABSRC),$(exchanged_tables))

$(exchanged_tables_final): $(exchanged_tables)
	cp $< $@
	jupyter nbconvert --execute $@

# Classifier results
classifier_results = $(EXPDIR)/2018-01-19.aspiration_classifiers_results_final.ipynb
cls_results_final = $(subst $(EXPDIR),$(TABSRC),$(classifier_results))

$(cls_results_final): $(classifier_results) $(rf_summaries) $(rf_summaries_exchange)
	cp $< $@
	jupyter nbconvert --execute $@

TABLES : $(sites_final_nb) $(demographics_final_nb) $(exchanged_tables_final) $(cls_results_final)

#TABLES: $(ppi_final_nb) $(site_comp_final_nb) $(exchanged_tables_final) #$(patient_tables_final)

###########################
####    SUPP FILES    #####
###########################

SUPPDIR = final/supp_files

################ Patients ################
# We need to keep track of which patients were used in which analyses
patient_list = $(SUPPDIR)/patients_with_sites_sampled.csv
# This notebook combines all the individual patient files and
# writes the final csv
patient_nb = $(EXPDIR)/2019-01-21.collate_all_patients_final.ipynb
patient_nb_final = $(subst $(EXPDIR),$(TABSRC),$(patient_nb))

# Okay, this recipe is going to be wonky because I don't want to re-run
# all of the analysis scripts, but that is what they depend on. Each
# analysis script writes its own (hard-coded) file name with the list
# of samples used in that analysis. Let's just keep this rule as is and
# hope that make runs things in the right order...

# - make_classifiers.py and make_classifiers_exchanged_OTUs.py run the
#   RF analyses and also write the sample IDs to files
# - exchange_track_samples.py only write the sample IDs from a command-line
#   call with the file path stem specified; it's re-run whenever
#   exchange_tidy.py is changed (or any of its dependencies)
# - prevalence_exchange_track_patients.py also writes the sample IDs only
#  from a command-line call (which specifies the file name stem); it
#  depends on prevalence_exchange.py (and its dependencies)
# - all of the jupyter notebooks depend on many other things, but we'll just
#   put the final versions of these notebooks here and hope that's good...
#   - 2017-12-26.within_patient_beta_div_figure_final.ipynb writes
#     figure3.within_patient_beta_div.samples.txt
#   - 2018-01-26.classify_sites_final.ipynb writes
#     figure3.site_classifiers.samples.txt
#   - 2018-01-03.asp_vs_nonasp_beta_fig_final.ipynb writes
#     figure4.asp_vs_nonasp_beta.samples.txt
#   - 2018-01-02.pca_and_between_patient_div_final.ipynb writes
#     figure2.2016_pcoa.samples.txt, figure2.2014_pcoa.samples.txt,
#     figure2.between_patient_jsd.samples.txt,
#     figure2.between_patient_jsd_2016_permanova.samples.txt,
#     figure2.between_patient_jsd_2014_permanova.samples.txt
#   - 2018-02-20.community_overview_ordered_by_aspiration_final.ipynb
#     writes figure1.overview_plots.samples.txt
#   - 2018-01-23.reflux_correlations_final.ipynb writes
#     figure5.reflux_correlation.samples.txt
patient_writing_files = $(SRCDIR)/make_loo_classifiers.py \
	$(SRCDIR)/make_classifiers_exchanged_OTUs.py \
	$(SRCDIR)/exchange_track_samples.py \
	$(SRCDIR)/prevalence_exchange_track_patients.py \
	$(FIGSRC)/2017-12-26.within_patient_beta_div_figure_final.ipynb \
	$(FIGSRC)/2018-01-26.classify_sites_final.ipynb \
	$(FIGSRC)/2018-01-03.asp_vs_nonasp_beta_fig_final.ipynb \
	$(FIGSRC)/2018-01-02.pca_and_between_patient_div_final.ipynb \
	$(FIGSRC)/2018-02-20.community_overview_ordered_by_aspiration_final.ipynb \
	$(FIGSRC)/2018-01-23.reflux_correlations_final.ipynb \
	$(clean_otu) \
	$(clean_meta)

$(patient_list): $(patient_nb_final) $(patient_writing_files)
	jupyter nbconvert --execute $<

$(patient_nb_final): $(patient_nb)
	cp $< $@

################ Clinical metadata ################
# This has the patient-wise cleaned up clinical metadata and sites that were sampled. This is the metadata that Rachel should give anyone who asks for data.
## the clinical_metadata file is patient-wise, the notebook also makes a
# sample-wise version that's in data/clean/rosen.clinical_metadata.clean
clinical_metadata := $(SUPPDIR)/patient_clinical_metadata.csv
clinical_meta_nb := $(EXPDIR)/2019-02-03.write_new_metadata_file.ipynb
clinical_meta_nb_final := $(subst $(EXPDIR),$(TABSRC),$(clinical_meta_nb))

$(clinical_metadata): $(clinical_meta_nb_final) $(clean_meta) $(clean_otu)
	jupyter nbconvert --execute $<

$(clinical_meta_nb_final): $(clinical_meta_nb)
	cp $< $@

################ qvalues ################
qvalues_supp_file = $(subst $(RESULTSDIR),$(SUPPDIR),$(qvalues))

$(qvalues_supp_file): $(qvalues)
	cp $< $@

################ SRA files ################
sra_nb = $(EXPDIR)/2018-03-12.upload_data.ipynb
sra_nb_final = $(subst $(EXPDIR),$(TABSRC),$(sra_nb))
sra_file = $(SUPPDIR)/sra_metadata.SUB3758953.txt

# The notebook makes two SRA metadata files
# It technically depends on the intermediate sample files in
# data/patients/*.samples.txt, but we'll put in the final patient_list
# (in supp_files/) as the dependency instead.
$(sra_file): $(sra_nb_final) $(patient_list)
	jupyter nbconvert --execute $<

$(sra_nb_final): $(sra_nb)
	cp $< $@

SUPP_FILES: $(patient_list) $(qvalues_supp_file) $(sra_file)
