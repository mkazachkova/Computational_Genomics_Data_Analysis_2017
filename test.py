import numpy as np
from sklearn import linear_model
from sklearn.model_selection import PredefinedSplit, cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import mean_squared_error
from sklearn.svm import SVR
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.feature_selection import VarianceThreshold
import operator
from sklearn.preprocessing import normalize  
from regression import *


print "LINEAR REGRESSION"
#calculate R^2 wth test sets
print "--fit quality--"
print "R squared value SNP: " + str(lin_fit.score(snp_result_test, y_test))
print "R squared value GENE: " + str(lin2_fit.score(gene_result_test, y_test2)) 

#cross-validation
lins = cross_val_score(lin_fit, snp_result_train, y_train, cv=10)  
lins2 = cross_val_score(lin2_fit, gene_result_train,  y_train, cv = 10) 
print lins
print lins2
print("CV Score: %0.2f (+/- %0.2f)" % (lins.mean(), lins.std() * 2))
print("CV Score: %0.2f (+/- %0.2f)" % (lins2.mean(), lins2.std() * 2))

#calculate MSE
lin_predict = lin_fit.predict(snp_result_test)
lin2_predict = lin2_fit.predict(gene_result_test)

print "--prediction quality--"
print "mean square error of linear regression on SNP is " + str(mean_squared_error(y_test, lin_predict)) #y_true, y_predict
print "mean square error of linear regression on gene expression is " + str(mean_squared_error(y_test2,lin2_predict))



print "RIDGE REGRESSION"
#calculate R^2 with test sets
print "--fit quality--"
print "R squared value SNP: " + str(ridge1_fit.score(snp_result_test, y_test)) 
print "R squared value GENE: " + str(ridge2_fit.score(gene_result_test, y_test2)) 

#cross-validation 
r = cross_val_score(ridge1_fit, snp_result_train,  y_train, cv=10) 
r2 = cross_val_score(ridge2_fit, gene_result_train,  y_train, cv = 10) 
print r
print r2
print("CV Score: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("CV Score: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 

#calculate MSE
ridge_predict = ridge1_fit.predict(snp_result_test)
ridge2_predict = ridge2_fit.predict(gene_result_test)

print "--prediction quality--"
print "mean square error of ridge regression on SNP is " + str(mean_squared_error(y_test, ridge_predict)) #y_true, y_predict
print "mean square error of ridge regression on gene expression is " + str(mean_squared_error(y_test2,ridge2_predict))



print "LASSO MODEL"

#calculate R^2 with test sets
print "--fit quality--"
print "R squared value SNP: " + str(lasso1_fit.score(snp_result_test, y_test)) 
print "R squared value GENE: " + str(lasso2_fit.score(gene_result_test, y_test2)) 

#cross-validation 
r = cross_val_score(lasso1_fit, snp_result_train,  y_train, cv=10) 
r2 = cross_val_score(lasso2_fit, gene_result_train,  y_train, cv = 10) 
print r
print r2
print("CV Score: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("CV Score: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 

#calculate MSE
lasso1_predict = lasso1_fit.predict(snp_result_test)
lasso2_predict = lasso2_fit.predict(gene_result_test)

print "--prediction quality--"
print "mean square error of lasso on SNP is " + str(mean_squared_error(y_test, lasso1_predict)) #y_true, y_predict
print "mean square error of lasso on gene expression is " + str(mean_squared_error(y_test2,lasso2_predict))


print "SUPPORT VECTOR MACHINE REGRESSION"
#calculate R^2 with test sets
print "--fit quality--"
print "R squared value SNP: " + str(svr1_fit.score(snp_result_test, y_test)) 
print "R squared value GENE: " + str(svr2_fit.score(gene_result_test, y_test2)) 

#cross-validation 
r = cross_val_score(svr1_fit, snp_result_train,  y_train, cv=10) 
r2 = cross_val_score(svr2_fit, gene_result_train,  y_train, cv = 10) 
print r
print r2
print("CV Score: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("CV Score: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 

#calculate MSE
svr1_predict = svr1_fit.predict(snp_result_test)
svr2_predict = svr2_fit.predict(gene_result_test)

print "--prediction quality--"
print "mean square error of svr on SNP is " + str(mean_squared_error(y_test, svr1_predict)) #y_true, y_predict
print "mean square error of svr on gene expression is " + str(mean_squared_error(y_test2,svr2_predict))
