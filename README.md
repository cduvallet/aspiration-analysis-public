# The oropharyngeal microbiome shapes lung microbial communities in children with oropharyngeal dysphagia and aspiration

Analyses, files, and figures accompanying Duvallet et al (2018)'s paper "The oropharyngeal microbiome shapes lung microbial communities in children with oropharyngeal dysphagia and aspiration".

The OTU tables and metadata used in these analyses are available upon request, and are not included in the public version of this repo.

# Directory structure

The structure of this repo follows loosely follows [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)'s recommended structure.

```
|-- Makefile   
|-- README.md
|
|-- data: data is available upon request
|
|-- final
|    |-- patients: text files with patient IDs used for each figure
|    |
|    |-- tables: empty folder pointing user to src/tables for notebooks
|    |
|    |-- figures: png files used to assemble final figures
|    
|-- src
|    |-- analysis: scripts to do all analyses
|    |    |-- including JSD, classifiers, exchange, prevalence of exchange
|    |    |-- scripts to calculate alpha diversity, q-values, and
|    |        unifrac distances are not used in the paper
|    |
|    |-- data: scripts to wrangle metadata and clean up OTU tables/metadata
|    |
|    |-- exploration: ipython notebooks with many exploratory analyses. Final
|    |   analyses are in notebooks labeled `_final`. These notebooks are also
|    |   copied to the `src/figures` folder. Note to myself: final notebooks
|    |   need to be manually re-run in this folder so that their content
|    |   reflects any updates. (`make` turns these into html's but does not
|    |   update their actual content).
|    |
|    |-- figures: notebooks used to produce all figures
|    |
|    |-- tables: notebooks used to produce all tables
|    |
|    |-- util: some useful functions
```
