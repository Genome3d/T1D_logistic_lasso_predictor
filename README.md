# T1D_logistic_lasso_regression_model
Type 1 diabetes (T1D) logistic lasso regression predictor model

Python 3.7.3
tsfresh                   0.12.0
pymc3                     3.8

unzip all the zip files under data/ before running any python scripts. 
Some regquired data may not be available. Please double check before running the python scripts.

This folder contains 7 python scripts for creating and validating the T1D lasso regression preditor developed from our study. The scripts should be run by python directly as 7 sequential steps. All the required files are location in data/ directory. The python_scripts_of_T1D_manuscript/ diectory contains all text files for recording the commands for developing the preditor model. They are only used for references. The python_scripts_of_T1D_manuscript/UKBiobank_Test_data/ directory contains the information and scripts for the T1D lasso regression predictor model with the UKBio bank derived test data. For using the the default files, your should start with B_python_script_split_std_80-20_WTCCC_T1D_weight_matrix.py.

You can start with running A_python_script_for_creating_WTCCC_T1D_weighted_matrix.py (optional):
Before running A_python_script_for_creating_WTCCC_T1D_weighted_matrix.py to create a individual tissue specific eQTL effect table, it is required to have your own inidividual genetype data in plink bed format with three files (bed,fam,bin) to be converted to a .raw text file. The pink command is as the following:

plink -bfile your_genotype_bed_file --recode A  --out individual_genotype_table

After this command, individual_genotype_table.raw is created and be placed in data/. Due to the data sharing agreements, we could not provide any indiviual genotype data from Wellcome Trust Case and Control Consortium and UK Biobank. Nevertheless, you can go to https://www.wtccc.org.uk/ and https://www.ukbiobank.ac.uk/ for applying their genotype data. 

Or you can start with skipping A_python_script_for_creating_WTCCC_T1D_weighted_matrix.py and running B_python_script_split_std_80-20_WTCCC_T1D_weight_matrix.py using default files at data/.

P.S.: Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019 is the Suppliementary table 5 in the T1D manuscript.
