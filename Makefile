## Makefile to reproduce data cleaning, analyses, and figures
## for aspiration paper
# Author: Claire Duvallet

# Define some directory shortcuts
RESULTSDIR = data/analysis
SRCDIR = src/analysis


all: DATA ANALYSIS FIGURES TABLES #TREE

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
widejsd = $(RESULTSDIR)/jsd.wide.txt
refluxjsd = $(RESULTSDIR)/jsd.tidy_reflux_bile.txt

# Unifrac
weighted_unifrac = $(RESULTSDIR)/weighted_unifrac.txt
unweighted_unifrac = $(RESULTSDIR)/unweighted_unifrac.txt
beta_div: $(jsd) #$(weighted_unifrac) $(unweighted_unifrac)

##TODO: check reflux correlation codes
#refluxcorrs = $(RESULTSDIR)/reflux_balgastric_corrs.txt

# Alpha diversity
alpha := $(RESULTSDIR)/alpha_diversity.txt

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
rf_rocs = $(RESULTSDIR)/rf_results.rocs.txt
rf_predictions = $(RESULTSDIR)/rf_results.predictions.txt
# Based on exchanged OTUs only
rf_summaries_exchange = $(RESULTSDIR)/rf_results.exchanged.summaries.txt
rf_rocs_exchange = $(RESULTSDIR)/rf_results.exchanged.rocs.txt
rf_predictions_exchange = $(RESULTSDIR)/rf_results.exchanged.predictions.txt

CLASSIFIERS: $(rf_summaries) $(rf_summaries_exchange)

ANALYSIS: beta_div exchange

## Not used
#$(univariate): src/analysis/get_qvalues.py $(clean_otu) $(clean_metadata)
#	python src/analysis/get_qvalues.py $(clean_otu) $(clean_metadata) $(univariate)

$(jsd): $(SRCDIR)/jsd.py $(clean_otu) $(clean_metadata)
	python src/analysis/jsd.py $(clean_otu) $(clean_metadata) $(jsd)

$(widejsd): $(jsd)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

$(weighted_unifrac): $(SRCDIR)/unifrac.py $(clean_otu_counts) $(clean_metadata) $(tree)
	python $< $(clean_otu_counts) $(tree) $(clean_metadata) $@ --weighted

$(unweighted_unifrac): $(SRCDIR)/unifrac.py $(clean_otu_counts) $(clean_metadata) $(tree)
	python $< $(clean_otu_counts) $(tree) $(clean_metadata) $@

$(exchange_w_partial): $(SRCDIR)/exchange_tidy.py $(clean_otu) $(clean_metadata)
	python $< $(clean_otu) $(clean_metadata) $@

$(null_exchange_w_partial): $(SRCDIR)/exchange_null.py $(clean_otu) $(clean_metadata)
	python $< $(clean_otu) $(clean_metadata) $@ --nthresh 10 --nshuffle 20

# Label the exchanged OTUs in a separate file
$(exchanged_otus): $(SRCDIR)/label_exchanged_otus.py $(exchange_w_partial)
	python $< $(exchange_w_partial) $@ --nthresh 10 --qthresh 0.1

# Calculate the prevalence of the exchanged OTUs. Note that thresholds should
# match the ones used in the recipe to label them
$(exchange_prevalence): $(SRCDIR)/prevalence_exchange.py $(exchange_w_partial) $(clean_otu) $(clean_metadata)
	python $< $(exchange_w_partial) $(clean_otu) $(clean_metadata) $@ --nthresh 10 --qthresh 0.1

# This also makes all of the other RF files
$(rf_summaries): $(SRCDIR)/make_classifiers.py $(clean_otu) $(clean_metadata)
	python $< $(clean_otu) $(clean_metadata) $(rf_summaries) $(rf_rocs) $(rf_predictions) --rfreps 100

# RF exchanged
$(rf_summaries_exchange): $(SRCDIR)/make_classifiers_exchanged_OTUs.py $(clean_otu) $(clean_metadata) $(exchange_prevalence)
	python $< $(clean_otu) $(clean_metadata) $(exchange_prevalence) $(rf_summaries_exchange) $(rf_rocs_exchange) $(rf_predictions_exchange)

# Alpha diversity
$(alpha): $(SRCDIR)/alpha_diversity.py $(clean_otu_counts) $(clean_metadata)
	python $< $(clean_otu_counts) $(clean_metadata) $@

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

# Figure 7 and 8 - AUCs and fisher p values for exchanged OTUs classifiers
exch_aucs_fig = $(FIGDIR)/figure7.auc_asp_classifiers.exchangedOTUs.png
exch_p_fig = $(FIGDIR)/figure8.fisherp_asp_classifiers.exchangedOTUs.png
#exch_cls_pts = $(PTDIR)/figure7.exchangedOTUs_classifiers.subjects.txt
# Notebooks
cls_nb = $(EXPDIR)/2018-01-19.aspiration_classifiers_results_final.ipynb
cls_final_nb = $(subst $(EXPDIR),$(FIGSRC),$(cls_nb))

EXCH_CLS_FIG: $(exch_aucs_fig) $(exch_p_fig)

# Figure 9 - AUCs and p values for full-abundance classifiers
asp_cls_fig = $(FIGDIR)/figure9.auc_fisherp_asp_classifiers.png
# Notebook: this one is also made by cls_nb
ASP_CLS_FIG: $(asp_cls_fig)

FIGURES: OVERVIEW_FIG ACROSS_FIG WITHIN_FIG ASP_FIG REFLUX_FIG EXCH_CLS_FIG ASP_CLS_FIG

################################
#######  FIG RECIPES    ########
################################

# Figure 1 - community overview
$(overview_fig): $(overview_final_nb) $(clean_otu) $(clean_metadata)
	jupyter nbconvert --execute $< --ExecutePreprocessor.timeout=-1

$(overview_final_nb): $(overview_nb)
	cp $< $@

# Figure 2 - across-patient beta diversity and PCoA plot
$(pca_fig): $(across_final_nb) $(clean_otu) $(clean_metadata) $(widejsd)
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

# Figure 7 and 8 - AUCs and fisher p values for exchanged OTUs classifiers
# Also update the notebook with the tables of exchanged OTUs
$(exch_aucs_fig): $(cls_final_nb) $(rf_summaries_exchange)
	jupyter nbconvert --execute $<

$(exch_p_fig): $(exch_aucs_fig)
	@if test -f $@; then :; else \
	  rm -f $<; \
	  make $<; \
	fi

$(cls_final_nb): $(cls_nb)
	cp $< $@

# Figure 9 - AUCs and p values for full-abundance classifiers
$(asp_cls_fig): $(cls_final_nb) $(rf_summaries)
	jupyter nbconvert --execute $<

###########################
####      TABLES      #####
###########################

# Where to write code that makes tables
TABSRC = src/tables

# Notebook that quickly checks for PPI confounding
ppi_nb = $(EXPDIR)/2018-02-14.ppi_confounding_check_final.ipynb
ppi_final_nb = $(subst $(EXPDIR),$(TABSRC),$(ppi_nb))

$(ppi_final_nb): $(ppi_nb)
	cp $< $@
	jupyter nbconvert --execute $@

# Body sites within people more similar than across people
site_comp_nb = $(EXPDIR)/2018-02-08.lung_gastric_vs_lung_lung_final.ipynb
site_comp_final_nb = $(subst $(EXPDIR),$(TABSRC),$(site_comp_nb))

$(site_comp_final_nb): $(site_comp_nb)
	cp $< $@
	jupyter nbconvert --execute $@

# Tables of exchanged bugs
exchanged_tables = $(EXPDIR)/2018-01-25.exchanged_OTUs_tables_final.ipynb
exchanged_tables_final = $(subst $(EXPDIR),$(TABSRC),$(exchanged_tables))

$(exchanged_tables_final): $(exchanged_tables)
	cp $< $@
	jupyter nbconvert --execute $@

# Number of patients
patient_tables = $(EXPDIR)/2018-03-06.number_of_patients_final.ipynb
patient_tables_final = $(subst $(EXPDIR),$(TABSRC),$(patient_tables))

$(patient_tables_final): $(patient_tables)
	cp $< $@
	jupyter nbconvert --execute $@

TABLES: $(ppi_final_nb) $(site_comp_final_nb) $(exchanged_tables_final) $(patient_tables_final)
