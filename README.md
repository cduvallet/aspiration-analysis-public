# The oropharyngeal microbiome shapes lung microbial communities in children with oropharyngeal dysphagia and aspiration

Analyses, files, and figures accompanying Duvallet et al (2018)'s paper "The oropharyngeal microbiome shapes lung microbial communities in children with oropharyngeal dysphagia and aspiration".

The OTU tables and metadata used in these analyses are available upon request, and are not included in the public version of this repo.

For most users, what you want will be in the `src/figures/` and `src/tables/` folders.
These contain the notebooks I used to generate the figures and tables in the paper.
You can view the iPython notebooks directly in github, or download the associated html's to view them that way.
If there are discrepancies between the html and notebooks, the html's have the most up-to-date results.
`src/analysis` may also be of interest: it contains the raw scripts I used for most analyses (e.g. JSD, classifiers, exchanged OTU definition, etc).

# Directory structure

The structure of this repo loosely follows [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)'s recommended structure.

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
     |-- analysis: scripts to do all analyses
     |    |-- including JSD, classifiers, exchange, prevalence of exchange
     |    |-- scripts to calculate alpha diversity, q-values, and
     |        unifrac distances are not used in the paper
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
