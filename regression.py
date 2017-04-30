import numpy as np
from sklearn import linear_model
from sklearn.model_selection import PredefinedSplit, cross_val_score
from sklearn import svm
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split

snp = np.load('final_snp_data.npy') #X1
exp = np.load('final_gene_data.npy') #X2
ics = np.load('final_ics.npy')  # y value
cell_lines = np.load('final_cell_lines.npy')


snp = snp.astype(np.float)
exp = exp.astype(np.float)
ics = ics.astype(np.float)

# split training and test sets
X_train, X_test, y_train, y_test = train_test_split(snp, ics, test_size = 0.2, random_state=0)
X_train2, X_test2, y_train2, y_test2 = train_test_split(exp, ics, test_size = 0.2, random_state=0)


#linear regression
lin = linear_model.LinearRegression()
lin_fit = lin.fit(X_train, y_train)
lin2 = linear_model.LinearRegression()
lin2_fit = lin2.fit(X_train2, y_train2)

#calculate R^2 wth test sets
print lin_fit.score(X_test, y_test) #-0.904163717163
print lin2_fit.score(X_test2, y_test2) #-0.172005711613

#cross-validation
print cross_val_score(lin, snp, ics)  #[-0.90431656 -0.91727751 -0.75193867]
print cross_val_score(lin2, exp, ics) #[-0.07655055 -0.25684399 -0.31656965]


#ridge regression
ridge1 = linear_model.Ridge()
ridge1_fit = ridge1.fit(X_train, y_train)
ridge2 = linear_model.Ridge()
ridge2_fit = ridge2.fit(X_train2, y_train2)

#calculate R^2 with test sets
print ridge1.score(X_test, y_test)  #-0.894554561435
print ridge2.score(X_test2, y_test2) #-0.171934191368

#cross-validation 
print cross_val_score(ridge1, snp, ics) #[-0.89550766 -0.79916886 -0.75113175]
print cross_val_score(ridge2, exp, ics) #[-0.07631259 -0.25673478 -0.31650685]
  

