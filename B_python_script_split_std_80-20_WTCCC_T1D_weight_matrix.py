# Split the tissue specific eQTL matrix into two files (80% and 20%)
# Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019.txt is the individual tissue specific eQTL effect table matrix with 4893 individual samples (the Supplementary Table 5)


import numpy as np
import pandas as pd

T1D_eQTL_table = pd.read_csv("data/Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019_v2.txt", sep="\t")


T1D_eQTL_table_tmp = T1D_eQTL_table.copy().reindex( np.random.permutation(T1D_eQTL_table.index))

T1D_eQTL_table2 = T1D_eQTL_table_tmp.copy().reindex( np.random.permutation(T1D_eQTL_table_tmp.index))


# The split is appled to the total number of the individual genotype samples. Our Type 1 diabetes data have 4893 individuls. 80% of them is 3914 and 20% is 979. You may need to adjust the number according to the number of you genotype samples.  
# Choose the first 3914 (out of 4893) examples for training.
std_table80 = T1D_eQTL_table2.head(3914)


# Choose the last 979 (out of 4893) examples for validation.
std_table20 = T1D_eQTL_table2.tail(979)

std_table80.to_csv('data/std80_Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019_v2.txt', sep='\t', index=False)
std_table20.to_csv('data/std20_Denis_total_Gwas_cat_Denis_2017_all_combined_eQTL_table02042019_v2.txt', sep='\t', index=False)
