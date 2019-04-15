This folder contains the supplementary files associated with the paper.

## SRA information

**biosample_attributes.SUB3758953.txt** and **sra_metadata.SUB3758953.txt** contain the metadata submitted to the SRA.

## Patient metadata

**combined_patient_site_and_analyses.csv** describes which patients were used in which analyses and figures, and which patients have which sites sequenced.

**patients_used_in_each_analysis.csv** indicates which patients were used in which analyses and figures (no site indication).

**patients_with_sites_sampled.csv** indicates which patients have which sites sequenced.

**patient_clinical_metadata.csv** is the clinical metadata associated with each patient. It also indicates which sites have sequencing data for each patient. Note that these metadata were consolidated across multiple clinical studies (which is why each column has the `_all` suffix, which I added in the wrangling and consolidating process).

## Analysis results

**qvalues.kruskal_wallis.txt** contains the calculated p- and q-values for the aspiration vs. non-aspiration comparison for each taxa.
