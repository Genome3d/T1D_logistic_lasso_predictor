# To validate the final T1D logistic lasso regression predictor model with 30 UKBiobank test datasets

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn import linear_model
import pandas as pd
from sklearn.model_selection import GridSearchCV
from tsfresh import extract_features, select_features
from  tsfresh.feature_selection.relevance import calculate_relevance_table
import pickle
import sklearn.metrics as metrics

# Selected_by_Fdr0.2_mann_on_Train_MARCH1.txt is the data_fields_selected by FDR0.2 Mann Whitney U test from WTCCC traning dataset 
eQTL_table20_Mann = pd.read_table('Selected_by_Fdr0.2_mann_on_Train_MARCH1.txt')
X = eQTL_table20_Mann[eQTL_table20_Mann.columns[6:]]


#Retrieve the saved final T1D lasso regression predictor model created from the full WTCCC traning dataset
filename = 'finalized_lgl1c1max500ran1T1D_model_2042.sav'
loaded_model = pickle.load(open(filename, 'rb'))

lg_clf = loaded_model



a = pd.read_csv('UKBio_T1D_case.table', sep='\t')
a['PHENOTYPE'] = 1

num = [1,2,3,4,5,6,7,8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

AUCs = []

for x in num:
	con_file =  "UKBio_T1D_control" + str(x) + ".table"
	b = pd.read_csv(con_file, sep='\t')
	b['PHENOTYPE'] = 0
	full_data = pd.concat([a,b],ignore_index=True)
	eQTL_table_test = full_data
	x_test = eQTL_table_test[eQTL_table_test.columns[6:]]
	for fe in X.columns:
		if fe not in eQTL_table_test.columns[6:]:
			x_test[fe] = 0
	y_test = eQTL_table_test['PHENOTYPE'] 
	x_test2 = x_test[X.columns]
	AUCs.append(roc_auc_score(y_test, lg_clf.predict_proba(x_test2)[:,1]))

results = pd.DataFrame(AUCs,columns=['UKBIO T1D AUC'])
results.to_csv("30_UKBio_Controls_AUCs.txt", sep='\t', index=False)





