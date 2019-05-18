# Aerodigestive sampling reveals altered microbial exchange between lung, oropharyngeal, and gastric microbiomes in children with impaired swallow function

Analyses, files, and figures accompanying Duvallet et al (2019)'s paper "Aerodigestive sampling reveals altered microbial exchange between lung, oropharyngeal, and gastric microbiomes in children with impaired swallow function".

The preprint is available on bioRxiv: https://doi.org/10.1101/476580

## Data

The raw fastq data used in these analyses are available from the SRA at [BioProject PRJNA450850](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA450850).

OTU tables and metadata used in these analyses are available on [Zenodo record 2678107](https://doi.org/10.5281/zenodo.2678107). Metadata and supplementary files are also available in this repo, in `final/supp_files/`.

## Code

For most users, what you want will be in the `src/figures/` and `src/tables/` folders.
These contain the notebooks I used to generate the figures and tables in the paper.
You can view the iPython notebooks directly in github, or download the associated html's to view them that way.
I'm not sure if I went through and re-ran all the notebooks after reviews, so I recommend that you re-run the iPython notebooks yourself if there are any discrepancies between what's in the notebooks and what's in the paper.

`src/analysis` may also be of interest: it contains the raw scripts I used for most analyses (e.g. JSD, classifiers, exchanged OTU definition, etc).

The `Makefile` has most of the analyses I did, though I did not go through and make sure that it runs correctly all the way and makes all the files in the final version of the paper (reviews came back after I graduated from my PhD, so I unfortunately didn't have time to vet this as thoroughly as needed). However, if you need to trace back which script was used to generate which figure(s), looking through the Makefile might come in handy.

If you have questions about any of these files or analyses, please email me (cduvallet at gmail.com).

# Reproducing analyses

There have been some backward-incompatible changes since I started writing this code, so you'll need to run these analyses in a virtual environment with older versions of the packages I used (especially scikit learn).

Do this by creating a virtual environment with the packages specified in `environment.yaml`. If you use conda, you can do this with the following command:

```
conda env create -f environment.yaml
source activate aspiration
```

If you have all the required data, you should be able to run `make` to re-make the paper. (Though, as described above, I am not 100% sure that every analysis gets run by `make`.)

Before doing that, however, you should run `python test_sklearn.py` and double-check that your python is using sklearn version 0.18 or lower.

## Running make from a tmux or screen session

If you want to run `make` from within a [tmux](https://github.com/tmux/tmux) or screen session, you'll need to open the session _first_ and then `source activate aspiration`.
Otherwise (if you activate your environment first and then open your tmux session), your call to `python` will look through the wrong directory for the packages. I have no idea why this happens but found the solution [here](https://stackoverflow.com/questions/50591901/screen-inside-the-conda-environment-doesnt-work).

## Opening notebooks interactively

To be able to run the notebooks interactively (i.e. by opening them after running `jupyter notebook`), you also need to install ipykernel:

```
conda install ipykernel
python -m ipykernel install --user --name aspiration --display-name "Python (aspiration)"
```

Now, when you open jupyter notebook you can change to this kernel from the top menu (`Kernel > Change kernel > Python (aspiration)`).

# Directory structure

The structure of this repo loosely follows [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)'s recommended structure.

```
|-- Makefile   
|-- README.md
|
|-- data: README with link to data files
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
