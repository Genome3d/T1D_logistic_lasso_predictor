# create a genotype table
#
# Assume your_genotype_bed_file is your plink genotype data 
# your_genotype_bed_file is in plink bed format with three files (bed,fam,bin)
# It needs to be converted to a text table format (.raw)

# plink -bfile your_genotype_bed_file --recode A --out individual_genotype_table



#--------------------------------------------------------------------------------------
# python program for creating the Tissue specific eQTL expression table for new corrected Denis_total_SNPs_sig_eQTL_effects 
#

import numpy as np
import pandas as pd

#Tissue_SNP_Gene_mapping.txt is a mapping file created from the Supplementary Table 4 of the T1D manuscript
#Tissue_SNP_Gene_mapping.txt is a mapping table of tissue specific eQTL-effects. It has four columns (Tissue, SNP, Gene_Name, Effect_Size)
#Tissue_SNP_Gene_mapping.txt is created from the significant tissue specific SNP-gene eQTL data from CodEs3D
TSG_list = pd.read_table('data/Tissue_SNP_Gene_mapping.txt')

#individual_genotype_table.raw is created from your plink bed genotype data (README.md)
eQTL_table = pd.read_table('data/individual_genotype_table.raw',sep = " ")

tmp_eQTL_table = eQTL_table.copy()
samples_header = list (eQTL_table.columns[0:6])
SNP_w_list = list(eQTL_table.columns[6:])
tmp_SNP_w_list = SNP_w_list.copy()
TSG_cols = []
drop_SNPs = []


# check the T1D SNPs if they are in the mapping table. If Yes, mutiple (weight) them with their tissue specific eQTL effect sizes from GTEx. 
for rsSNP_w in tmp_SNP_w_list :
	rsSNP = rsSNP_w.split('_')[0]
	tmp_TSGs = list ( TSG_list.Tissue[TSG_list.SNP == rsSNP] + "--" + rsSNP_w + "--" + TSG_list.Gene_Name[TSG_list.SNP == rsSNP])
	tmp_Effects = list ( TSG_list.Effect_Size[TSG_list.SNP == rsSNP])
	if len(tmp_Effects) != 0 :
		count = 0
		for TSG in tmp_TSGs:
			tmp_eQTL_table[TSG] = tmp_eQTL_table[rsSNP_w] * tmp_Effects[count]
			count = count + 1
			TSG_cols.append(TSG)
		drop_SNPs.append(rsSNP_w)

for rsSNP_w in drop_SNPs:
	SNP_w_list.remove(rsSNP_w)
	
TSG_col_sorted = sorted(TSG_cols)
new_eQTL_table = tmp_eQTL_table[samples_header + TSG_col_sorted + SNP_w_list]
new_eQTL_table.to_csv('data/Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019_v2.txt', sep='\t', index=False)
# Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019.txt is the individual tissue specific eQTL effect table (referring to the Supplementary Table 5 of the T1D manuscript)


#--------------------------------------------------------------------------------------
