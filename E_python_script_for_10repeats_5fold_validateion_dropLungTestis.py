
## Use RepeatedKFold to do 10 x 5 fold validation of the T1D logistic lasso regression models with removing 'Lung--rs3087243_A--CTLA4' or 'Testis--rs3087243_A--CTLA4'
# The script was used to create the Supplementary Table 12 in the T1D manuscript

# Inputs: Weighted_eQTL_matrix.txt
# The full tissue specific eQTL table

# Outputs: sk_grid_lg_0.2_onTrain_50model_dropLung + str(count) + weights.txt where count = 1 to 50,  sk_grid_lg_0.2_onTrain_50model_dropTesis + str(count) + weights.txt where count 1 to 50
# 	   AUC_results_50modeldropLung.txt, AUC_results_50modeldropTestis.txt, Man0.2_onTrianBest_fullmax500_50modelsdropLungTestis.txt 
# sk_grid_lg_0.2_onTrain_50model_dropLung + str(count) + weights.txt contains the model components and wights for each randomized predictor without eQTL Lung--rs3087243_A--CTLA4
# sk_grid_lg_0.2_onTrain_50model_dropTesis + str(count) + weights.txt contains the model components and wights for each randomized predictor without eQTL Testis--rs3087243_A--CTLA4
# AUC_results_50modeldropLung.txt contains AUC results from the 50 predictors without 'Lung--rs3087243_A--CTLA4'
# AUC_results_50modeldropTestis.txt contains AUC results from the 50 predictors without 'Testis--rs3087243_A--CTLA4'
# Man0.2_onTrianBest_fullmax500_50modelsdropLungTestis.txt contains the performance results for the randomized predictor models

import numpy as np
import pandas as pd
import statistics as st
from tsfresh import extract_features, select_features
from  tsfresh.feature_selection.relevance import calculate_relevance_table
from tsfresh.transformers import FeatureSelector


import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn import linear_model
from sklearn.model_selection import RepeatedKFold

from sklearn.model_selection import GridSearchCV

eQTL_table = pd.read_table('data/Weighted_eQTL_matrix.txt')

Y_pd = eQTL_table['PHENOTYPE'] - 1

X_pd = eQTL_table[eQTL_table.columns[6:]]

#X = X_pd.values
#Y = Y_pd.values

X = X_pd
Y = Y_pd

AUC_result_dl = list() # This list will store your results
AUC_result_dt = list()
selected_features = list()
classifiers = list()
# Initialize tsfresh feature selector
# Initialize 10 times repeated 10-fold cross validation
rkf = RepeatedKFold(n_splits=5, n_repeats=10, random_state=1)
# This loop will iterate 100 times
# Every row of your data set will be in the test set at least in 10 iterations
# However, here we are assuming a fixed set of hyperparameters
# Please regard the following as pseudo-code
# I have omitted any hyperparameter configuration

clf_dl = LogisticRegression(C=1,max_iter=500, l1_ratio=1, random_state=1, solver='saga',n_jobs=-1, penalty='elasticnet')
clf_dt = LogisticRegression(C=1,max_iter=500, l1_ratio=1, random_state=1, solver='saga',n_jobs=-1, penalty='elasticnet')
select = FeatureSelector(fdr_level=0.2)

count = 0
for train_index, test_index in rkf.split(X):
# print("TRAIN:", train_index, "TEST:", test_index)
	X_train, X_test = X.iloc[train_index], X.iloc[test_index]
	y_train, y_test = Y.iloc[train_index], Y.iloc[test_index]
# 2. Apply Mann_Whitley test on the 90% training dataset
	
	X_train_sel= select.fit_transform(X_train, y_train)
# 3. Fit predictor
	X_train_seldropLung = X_train_sel.drop(columns=['Lung--rs3087243_A--CTLA4'])
	X_train_seldropTestis = X_train_sel.drop(columns=['Testis--rs3087243_A--CTLA4'])
	clf_dl.fit(X_train_seldropLung, y_train)
	clf_dt.fit(X_train_seldropTestis, y_train)
# 4. Validate predictor
# Use .transform() NOT .fit_transform()
	X_test_sel = select.transform(X_test)
	X_test_seldropLung = X_test_sel.drop(columns=['Lung--rs3087243_A--CTLA4'])
	X_test_seldropTestis = X_test_sel.drop(columns=['Testis--rs3087243_A--CTLA4'])
	y_pred_dl = clf_dl.predict_proba(X_test_seldropLung)[:,1]
	y_pred_dt = clf_dt.predict_proba(X_test_seldropTestis)[:,1]
	auc_dl = roc_auc_score(y_test, y_pred_dl)
	auc_dt = roc_auc_score(y_test, y_pred_dt)
# 5. Store result
	AUC_result_dl.append(auc_dl)
	AUC_result_dt.append(auc_dt)
#	selected_features.append(select.relevant_features)
#	classifiers.append(clf)

#	X_header_dl = np.array(X_train_seldropLung.columns)
#	X_header_dt = np.array(X_train_seldropTestis.columns)
#        best_clf = clf
#	best_clf = clf_dl
	count = count + 1
	X_header_dl = np.array(X_train_seldropLung.columns)
	best_clf = clf_dl
	data_array = np.vstack((X_header_dl,best_clf.coef_[0,:]))
	model_weights = pd.DataFrame(data=data_array.T,columns=['Data_feature', 'Weight'])
	m_name = 'sk_grid_lg_0.2_onTrain_50model_dropLung' + str(count) + 'weights.txt'
	model_weights.to_csv(m_name, sep='\t',index=False, line_terminator='\n')
	best_clf = clf_dt
	X_header_dt = np.array(X_train_seldropTestis.columns)
	data_array = np.vstack((X_header_dt,best_clf.coef_[0,:]))
	model_weights = pd.DataFrame(data=data_array.T,columns=['Data_feature', 'Weight'])
	m_name = 'sk_grid_lg_0.2_onTrain_50model_dropTesis' + str(count) + 'weights.txt'
	model_weights.to_csv(m_name, sep='\t',index=False, line_terminator='\n')

# Now do an in-sample evaluation
# 6. Apply Mann-Whitley test on the 100% original eQTL dataset.
X_sel = select.fit_transform(X, Y)
# 7. Fit predictor to statistically significant features (just once!!!)
X_sel_dl = X_sel.drop(columns=['Lung--rs3087243_A--CTLA4'])
clf_dl.fit(X_sel_dl, Y)
y_pred_dl = clf_dl.predict_proba(X_sel_dl)[:,1]
# This in-sample AUC should be better than your the AUCs from your repeated cross-validation
auc_dl = roc_auc_score(Y, y_pred_dl)

X_sel_dt = X_sel.drop(columns=['Testis--rs3087243_A--CTLA4'])
clf_dt.fit(X_sel_dt, Y)
y_pred_dt = clf_dt.predict_proba(X_sel_dt)[:,1]
auc_dt = roc_auc_score(Y, y_pred_dt)





# AUC results from the 50 predictors without 'Lung--rs3087243_A--CTLA4' or 'Testis--rs3087243_A--CTLA4'
AUC_out = pd.DataFrame(AUC_result_dl, columns=['AUC'])
AUC_out.to_csv("AUC_results_50modeldropLung.txt", sep='\t',index=False, line_terminator='\n')
AUC_out = pd.DataFrame(AUC_result_dt, columns=['AUC'])
AUC_out.to_csv("AUC_results_50modeldropTestis.txt", sep='\t',index=False, line_terminator='\n')

AUC_dlstd= st.stdev(AUC_result_dl)
AUC_dtstd= st.stdev(AUC_result_dt)
AUC_dlmean= st.mean(AUC_result_dl)
AUC_dtmean= st.mean(AUC_result_dt)

num_coef_dl = np.sum(clf_dl.coef_[0,:] != 0)
num_coef_dt = np.sum(clf_dt.coef_[0,:] != 0)


#Performance results for the randomized predictors
f= open("Man0.2_onTrianBest_fullmax500_50modelsdropLungTestis.txt","w+")

f.write('Man0.2_onTrianBest_fullmax500_50models_dropLungTestis\n')

f.write('In-Sample AUC_dropLung: ' + str(auc_dl) + '\n')
f.write('MeanCV AUC_dropLung: ' + str(AUC_dlmean) + '\n')
f.write('Standard Deviation CV AUC_dropLung: ' + str(AUC_dlstd) + '\n')

f.write('num_coef drop Lung: ' + str(num_coef_dl) + '\n\n')

f.write('In-Sample AUC_dropTestis: ' + str(auc_dt) + '\n')
f.write('MeanCV AUC_dropTestis: ' + str(AUC_dtmean) + '\n')
f.write('Standard Deviation CV AUC_dropTestis: ' + str(AUC_dtstd) + '\n')
f.write('num_coef drop Tesis: ' + str(num_coef_dt) + '\n\n')



f.close()


