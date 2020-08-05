# python script for bayesian estimation supersedes the t test from pymc3 
# for AUC difference analysis of 50 predictor models 'Lung--rs3087243_A--CTLA4' or 'Testis--rs3087243_A--CTLA4'


import numpy as np
import pymc3 as pm
import pandas as pd
import matplotlib.pyplot as plt

AUCdropTestis_table = pd.read_table('data/AUC_results_50modeldropTestis.txt')
AUCdropLung_table = pd.read_table('data/AUC_results_50modeldropLung.txt')

y1 = AUCdropTestis_table['AUC'].values
y2 = AUCdropLung_table['AUC'].values

y_diff = y1 - y2
y = y_diff

#y.hist('value', by='group', figsize=(12, 4));

#plt.show()

µ_m = np.mean(y)
µ_s = np.std(y) * 2


with pm.Model() as model:
    group1_mean = pm.Normal('group1_mean', mu=µ_m, sd=µ_s)
   
	
s_low = 0
s_high = 1

with model:
    group1_std = pm.Uniform('group1_std', lower=s_low, upper=s_high)
    

with model:
    v = pm.Exponential('v_minus_one', 1/29.) + 1

	
with model:
    a1 = group1_std**-2
    group1 = pm.StudentT('y_diff', nu=v, mu=group1_mean, lam=a1, observed=y_diff)
    


with model:
    trace = pm.sample(2000,tune=1000)


pm.plot_posterior(trace, var_names=['group1_mean', 'group1_std'], ref_val=0);

plt.show()
pm.summary(trace)



def major_formatter(x, pos):
    return "%.1e" % x

import matplotlib.ticker as ticker

ax, = pm.plot_posterior(trace, var_names=['group1_mean'], ref_val=0);	  

ax.xaxis.set_major_formatter(ticker.FuncFormatter(major_formatter))

plt.show()

ax, = pm.plot_posterior(trace, var_names=['group1_std'], ref_val=0);
ax.xaxis.set_major_formatter(ticker.FuncFormatter(major_formatter))
plt.show()
