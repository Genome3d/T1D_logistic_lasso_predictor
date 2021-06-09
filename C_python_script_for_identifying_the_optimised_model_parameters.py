
# To apply mann whitney u test from tsfresh to filter out data features > FDR 0.2 on the 80% individual tissue specific eQTL effect table (the Supplementary table 5)
# Use GridSearchCV to search the optimized hyperparameter for the T1D logistic predictor model
# Sklearn Machine Learning algorithms are used to create T1D predictive models from the individual tissue specific eQTL effect table
# Inputs: std80_Weighted_eQTL_matrix.txt, std20_Weighted_eQTL_matrix.txt
# The 80% and 20% individual tissue specific eQTL effect tables

# Outputs: sk_grid_lg_0.2_onTrain_saga_elnet2_20092019.txt, sk_grid_lg_0.2_onTrain_saga_elnet2_20092019weights.txt
# sk_grid_lg_0.2_onTrain_saga_elnet2_20092019.txt : Predictor performance information
# sk_grid_lg_0.2_onTrain_saga_elnet2_20092019weights.txt : Predictor model components and their weight information


import numpy as np
import pandas as pd


from tsfresh import extract_features, select_features
from  tsfresh.feature_selection.relevance import calculate_relevance_table

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn import linear_model

from sklearn.model_selection import GridSearchCV

eQTL_table = pd.read_table('data/std80_Weighted_eQTL_matrix.txt')
x_features = eQTL_table[eQTL_table.columns[6:]]
y_phenotype = eQTL_table['PHENOTYPE'] - 1


X_selected0_1 = select_features(x_features, y_phenotype,test_for_binary_target_real_feature='mann', test_for_real_target_binary_feature='mann', fdr_level=0.2 )

full_columns = list(eQTL_table.columns[:6]) + list( X_selected0_1.columns)
selected_std_full = eQTL_table[full_columns]


eQTL_table80 = selected_std_full
temp20 = pd.read_table('data/std20_Weighted_eQTL_matrix.txt')

eQTL_table20 = temp20[full_columns]

y_train = eQTL_table80['PHENOTYPE'] - 1
y_test = eQTL_table20['PHENOTYPE'] - 1
x_train = eQTL_table80[eQTL_table80.columns[6:]]
x_test = eQTL_table20[eQTL_table20.columns[6:]]


X = x_train
Y = y_train

parameters = {'C':[1e-4,1e-3,1e-2,1e-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1e0,3,5,7,8,10,15,20,25,30,40,50,100],'max_iter':[1,5,70,100,130,150,170, 180, 200, 300,500,1000,1200,1400,1600,1800,2000,2200,2400,2600,3000],'l1_ratio':[1,0.9,0.8,0.7,0.6]}
lg_clf = LogisticRegression(random_state=1, solver='saga',n_jobs=-1, penalty='elasticnet' )
grid_clf = GridSearchCV(lg_clf, parameters, scoring='roc_auc', n_jobs=-1,iid=False, cv=10)
grid_clf.fit(X,Y)

f= open("sk_grid_lg_0.2_onTrain_saga_elnet2_20092019.txt","w+")
#sk_grid_lg_0.2_onTrain_saga_elnet2_20092019.txt contains the predictor performance information

f.write('sk_grid_lg_0.2_onTrain_saga_elnet2_20092019\n')
f.write('grid_clf.best_score_: ' + str(grid_clf.best_score_) + '\n')
f.write('grid_clf.best_estimator_: ' + str(grid_clf.best_estimator_) + '\n')
f.write('grid_clf.best_params_: ' + str(grid_clf.best_params_) + '\n')
f.write('grid_clf.scorer_: ' + str(grid_clf.scorer_) + '\n')
f.write('grid_clf.best_score_: ' + str(grid_clf.best_score_) + '\n')


lg_clf_best_grid = grid_clf.best_estimator_

model_Max = roc_auc_score(Y, lg_clf_best_grid.predict_proba(X)[:,1])
test_score = roc_auc_score(y_test, lg_clf_best_grid.predict_proba(x_test)[:,1])
num_coef = np.sum(lg_clf_best_grid.coef_[0,:] != 0)



cv_results = grid_clf.cv_results_
std_cv = cv_results['std_test_score'][grid_clf.best_index_]

f.write('--------------------------------------------------------------------\n\n')
f.write('In-Sample AUC: ' + str(model_Max) + '\n')
f.write('MeanCV AUC: ' + str(grid_clf.best_score_) + '\n')
f.write('Standard Deviation CV AUC: ' + str(std_cv) + '\n')
f.write('Test sample AUC: ' + str(test_score) + '\n')
f.write('num_coef: ' + str(num_coef) + '\n')


f.close()

X_header = np.array(X.columns)
best_clf =  grid_clf.best_estimator_
data_array = np.vstack((X_header,best_clf.coef_[0,:]))
model_weights = pd.DataFrame(data=data_array.T,columns=['Data_feature', 'Weight'])
model_weights.to_csv('data/sk_grid_lg_0.2_onTrain_saga_elnet2_20092019weights.txt', sep='\t',index=False,line_terminator='\n')
#sk_grid_lg_0.2_onTrain_saga_elnet2_20092019weights.txt contains the components and their weights of the best predictor model from sklearn grid search algorithm


