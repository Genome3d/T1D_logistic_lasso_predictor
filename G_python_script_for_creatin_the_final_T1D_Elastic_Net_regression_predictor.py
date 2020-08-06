
# To create the final T1D logistic lasso regression predictor with the full the tissue specific eQTL effect matrix

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


eQTL_table100 = pd.read_table('data/Weighted_eQTL_matrix.txt')
eQTL_table100_selected = eQTL_table100[selected_std_full.columns]

y_train = eQTL_table100_selected['PHENOTYPE'] - 1
x_train = eQTL_table100_selected[eQTL_table100_selected.columns[6:]]

X = x_train
Y = y_train


parameters = {'C':[1],'max_iter':[500],'l1_ratio':[1]}
lg_clf = LogisticRegression(random_state=1, solver='saga',n_jobs=-1, penalty='elasticnet' )
grid_clf = GridSearchCV(lg_clf, parameters, scoring='roc_auc', n_jobs=-1,iid=False, cv=10)
grid_clf.fit(X,Y)


# The final T1D predictor performance
f= open("sk_grid_lg_0.2_onTrain_saga_c1l1max500_full_modelb_26092019.txt","w+")

f.write('sk_grid_lg_0.2_onTrain_saga_c1l1max500_full_modelb_26092019\n')
f.write('10 fold:\n')
f.write('grid_clf.best_score_: ' + str(grid_clf.best_score_) + '\n')
f.write('grid_clf.best_estimator_: ' + str(grid_clf.best_estimator_) + '\n')
f.write('grid_clf.best_params_: ' + str(grid_clf.best_params_) + '\n')
f.write('grid_clf.scorer_: ' + str(grid_clf.scorer_) + '\n')
f.write('grid_clf.best_score_: ' + str(grid_clf.best_score_) + '\n')
f.write('grid_clf.best_score_: ' + str(grid_clf.best_score_) + '\n')

lg_clf_best_grid = grid_clf.best_estimator_

model_Max = roc_auc_score(Y, lg_clf_best_grid.predict_proba(X)[:,1])

num_coef = np.sum(lg_clf_best_grid.coef_[0,:] != 0)

cv_results = grid_clf.cv_results_
std_cv = cv_results['std_test_score'][grid_clf.best_index_]

f.write('--------------------------------------------------------------------\n\n')
f.write('In-Sample AUC: ' + str(model_Max) + '\n')
f.write('10 fold MeanCV AUC: ' + str(grid_clf.best_score_) + '\n')
f.write('Standard Deviation CV AUC: ' + str(std_cv) + '\n')

f.write('num_coef: ' + str(num_coef) + '\n')
f.write('\n\n\n')


# The selected data fields and their model weights in the lasso regression predictor model 
# The model weights estimated by the final predictor are the contributions of the eQTL or SNP effects to T1D predictions 
X_header = np.array(X.columns)
best_clf =  grid_clf.best_estimator_
data_array = np.vstack((X_header,best_clf.coef_[0,:]))
model_weights = pd.DataFrame(data=data_array.T,columns=['Data_feature', 'Weight'])
model_weights.to_csv('sk_grid_lg_0.2_onTrain_saga_c1l1max500_full_modelb_26092019weights.txt', sep='\t',index=False,line_terminator='\n')




