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
import warnings

warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

snp = np.load('data_output/final_snp_data.npy') #X1
exp = np.load('data_output/final_gene_data.npy') #X2
ics = np.load('data_output/final_ics.npy')  # y value
cell_lines = np.load('data_output/final_cell_lines.npy')

snp = snp.astype(np.float)
exp = exp.astype(np.float)
ics = ics.astype(np.float)
snp_trans = np.transpose(snp)

# split training and test sets
X_train, X_test, y_train, y_test = train_test_split(snp, ics, test_size = 0.2, random_state=0)
X_train2, X_test2, y_train2, y_test2 = train_test_split(exp, ics, test_size = 0.2, random_state=0)


'''Pre-feature selection'''

'''linear regression'''
print '--linear regression--'
lin = linear_model.LinearRegression()
lin_fit = lin.fit(X_train, y_train)
lin2 = linear_model.LinearRegression()
lin2_fit = lin2.fit(X_train2, y_train2)

#calculate R^2 wth test sets
print 'R^2 for SNP:' + str(lin_fit.score(X_test, y_test)) 
print 'R^2 for gene expression:' + str(lin2_fit.score(X_test2, y_test2)) 

#calculate MSE
lin_predict = lin_fit.predict(X_test)
lin2_predict = lin2_fit.predict(X_test2)
print "mean square error of linear regression on SNP is " + str(mean_squared_error(y_test, lin_predict)) #y_true, y_predict
print "mean square error of linear regression on gene expression is " + str(mean_squared_error(y_test2,lin2_predict))

#cross-validation
lins = cross_val_score(lin, snp, ics, cv=10)  
lins2 = cross_val_score(lin2, exp, ics, cv = 10) 
print("CV for SNP: %0.2f (+/- %0.2f)" % (lins.mean(), lins.std() * 2))
print("CV for gene expression: %0.2f (+/- %0.2f)" % (lins2.mean(), lins2.std() * 2))

'''ridge regression'''
print '--ridge regression--'
ridge1 = linear_model.Ridge()
ridge1_fit = ridge1.fit(X_train, y_train)
ridge2 = linear_model.Ridge()
ridge2_fit = ridge2.fit(X_train2, y_train2)

#calculate R^2 with test sets
print 'R^2 for SNP' + str(ridge1.score(X_test, y_test))  
print 'R^2 for gene expression' + str(ridge2.score(X_test2, y_test2)) 

#calculate MSE
ridge_predict = ridge1_fit.predict(X_test)
ridge2_predict = ridge2_fit.predict(X_test2)

print "mean square error of ridge regression on SNP is " + str(mean_squared_error(y_test, ridge_predict)) #y_true, y_predict
print "mean square error of ridge regression on gene expression is " + str(mean_squared_error(y_test2,ridge2_predict))

#cross-validation 
r = cross_val_score(ridge1, snp, ics, cv=10) 
r2 = cross_val_score(ridge2, exp, ics, cv = 10) 
print("CV for SNP: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("CV for gene expression: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 


'''LASSO (on SNP)'''
print '--LASSO--'
clf = linear_model.Lasso(alpha=1)
clf_fit = clf.fit(X_train, y_train)
clf2 = linear_model.Lasso(alpha=1)
clf2_fit = clf2.fit(X_train2, y_train2)

#calculate R^2 with test sets
print "R^2 for lasso (SNP): " + str(clf_fit.score(X_test, y_test))  
print "R^2 for lasso (gene expression): " + str(clf2_fit.score(X_test2, y_test2)) 

#calculate MSE
clf_predict = clf_fit.predict(X_test)
clf2_predict = clf2_fit.predict(X_test2)

print "mean square error of LASSO on SNP is " + str(mean_squared_error(y_test, clf_predict)) #y_true, y_predict
print "mean square error of LASSO on gene expression is " + str(mean_squared_error(y_test2,clf2_predict))

#cross-validation 
r = cross_val_score(clf, snp, ics, cv=10) 
r2 = cross_val_score(clf2, exp, ics, cv = 10) 
print("CV for SNP: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("CV for gene expression: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 

'''SVR'''
print '--Support Vector Regression--'
svr = SVR()
svr_fit = svr.fit(X_train, y_train)
svr2 = SVR()
svr2_fit = svr2.fit(X_train2, y_train2)

#calculate R^2 with test sets
print "R^2 for SVR (SNP): " + str(svr_fit.score(X_test, y_test))  
print "R^2 for SVR (gene expression): " + str(svr2_fit.score(X_test2, y_test2)) 

#calculate MSE
svr_predict = svr_fit.predict(X_test)
svr2_predict = svr2_fit.predict(X_test2)

print "mean square error of SVR on SNP is " + str(mean_squared_error(y_test, svr_predict)) #y_true, y_predict
print "mean square error of SVR on gene expression is " + str(mean_squared_error(y_test2,svr2_predict))

#cross-validation 
r = cross_val_score(svr, snp, ics, cv=10) 
r2 = cross_val_score(svr2, exp, ics, cv = 10) 
print("CV for SNP: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("CV for gene expression: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 


print 'plotting variance over features on SNP data...'
#feature selection for SNP
selector = VarianceThreshold()
selector.fit(X_train)
all_var = selector.variances_
sorted_all_var = sorted(all_var, reverse=True)
y_pos = np.arange(22419)
plt.bar(y_pos,sorted_all_var,align='center')
plt.ylabel('Variances')
plt.xlabel('Features')
plt.title('Variance over Features on SNP Data')
plt.savefig('data_output/snp_var.png')
plt.close()

print 'plotting variance over features on gene expression data...'
#feature selection for gene expression
selector = VarianceThreshold()
selector.fit(X_train2)
all_var = selector.variances_
sorted_all_var = sorted(all_var, reverse=True)
y_pos = np.arange(18926)
plt.bar(y_pos,sorted_all_var,align='center')
plt.ylabel('Variances')
plt.xlabel('Features')
plt.title('Variance over Features on Gene Expression Data')
plt.savefig('data_output/gene_var.png')