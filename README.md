# Aerodigestive sampling reveals altered microbial exchange between lung, oropharyngeal, and gastric microbiomes in children with impaired swallow function

Analyses, files, and figures accompanying Duvallet et al (2019)'s paper "Aerodigestive sampling reveals altered microbial exchange between lung, oropharyngeal, and gastric microbiomes in children with impaired swallow function".

The preprint is available on bioRxiv: https://doi.org/10.1101/476580

## Data

The raw fastq data used in these analyses will be made available from the SRA at [BioProject PRJNA450850](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA450850).

OTU tables and metadata used in these analyses will be made available on Zenodo. Metadata will also be made available in this repo, in `final/supp_files/`.

## Code

For most users, what you want will be in the `src/figures/` and `src/tables/` folders.
These contain the notebooks I used to generate the figures and tables in the paper.
You can view the iPython notebooks directly in github, or download the associated html's to view them that way.
If there are discrepancies between the html and notebooks, the html's have the most up-to-date results.

`src/analysis` may also be of interest: it contains the raw scripts I used for most analyses (e.g. JSD, classifiers, exchanged OTU definition, etc).

The `Makefile` has most of the analyses I did, though I did not go through and make sure that it runs correctly all the way and makes all the files in the final version of the paper (reviews came back after I graduated from my PhD, so I unfortunately didn't have time to vet this as thoroughly as needed).

If you have questions about any of these files or analyses, please email me (cduvallet at gmail.com).

# Directory structure

The structure of this repo loosely follows [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)'s recommended structure.

```
|-- Makefile   
|-- README.md
|
|-- data: data is available upon request
|
|-- final
|    |-- supp_files: supplementary files
|    |
|    |-- tables: empty folder pointing user to src/tables for notebooks
|    |
|    |-- figures: png files used to assemble final figures
|    
|-- src
     |-- analysis: scripts to do all analyses
     |    |-- including JSD, classifiers, exchange, prevalence of exchange
     |    |-- alpha diversity, q-values, and unifrac distances
     |
     |-- data: scripts to wrangle metadata and clean up OTU tables/metadata
     |
     |-- exploration: ipython notebooks with many exploratory analyses. Final
     |   analyses are in notebooks labeled `_final`. These notebooks are also
     |   copied to the `src/figures` folder. Note to myself: final notebooks
     |   need to be manually re-run in this folder so that their content
     |   reflects any updates. (`make` turns these into html's but does not
     |   update their actual content).
     |
     |-- figures: notebooks used to produce all figures
     |
     |-- tables: notebooks used to produce all tables
     |
     |-- util: some useful functions
```
