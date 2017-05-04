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


snp = np.load('final_snp_data.npy') #X1
exp = np.load('final_gene_data.npy') #X2
ics = np.load('final_ics.npy')  # y value
cell_lines = np.load('final_cell_lines.npy')


snp = snp.astype(np.float)
exp = exp.astype(np.float)
ics = ics.astype(np.float)
snp_trans = np.transpose(snp)
print snp.shape
print snp_trans.shape

print 'snp dimension is ' + str(snp.shape)
print 'exp dimension is ' + str(exp.shape) 

# split training and test sets
X_train, X_test, y_train, y_test = train_test_split(snp, ics, test_size = 0.2, random_state=0)
X_train2, X_test2, y_train2, y_test2 = train_test_split(exp, ics, test_size = 0.2, random_state=0)


#print X_train.shape
#print y_train.shape

# #linear regression
# lin = linear_model.LinearRegression()
# lin_fit = lin.fit(X_train, y_train)
# lin2 = linear_model.LinearRegression()
# lin2_fit = lin2.fit(X_train2, y_train2)

# #calculate R^2 wth test sets
# print lin_fit.score(X_test, y_test) 
# print lin2_fit.score(X_test2, y_test2) 

# #calculate MSE
# lin_predict = lin_fit.predict(X_test)
# lin2_predict = lin2_fit.predict(X_test2)
# print "mean square error of linear regression on SNP is " + str(mean_squared_error(y_test, lin_predict)) #y_true, y_predict
# print "mean square error of linear regression on gene expression is " + str(mean_squared_error(y_test2,lin2_predict))

# #cross-validation
# lins = cross_val_score(lin, snp, ics, cv=10)  
# lins2 = cross_val_score(lin2, exp, ics, cv = 10) 
# print lins
# print lins2
# print("Accuracy: %0.2f (+/- %0.2f)" % (lins.mean(), lins.std() * 2))
# print("Accuracy: %0.2f (+/- %0.2f)" % (lins2.mean(), lins2.std() * 2))





# #ridge regression
# ridge1 = linear_model.Ridge()
# ridge1_fit = ridge1.fit(X_train, y_train)
# ridge2 = linear_model.Ridge()
# ridge2_fit = ridge2.fit(X_train2, y_train2)

# #calculate R^2 with test sets
# print ridge1.score(X_test, y_test)  
# print ridge2.score(X_test2, y_test2) 

# #calculate MSE
# ridge_predict = ridge1_fit.predict(X_test)
# ridge2_predict = ridge2_fit.predict(X_test2)

# print 'length of ridge predict is ' + str(len(ridge_predict))
# print "mean square error of ridge regression on SNP is " + str(mean_squared_error(y_test, ridge_predict)) #y_true, y_predict
# print "mean square error of ridge regression on gene expression is " + str(mean_squared_error(y_test2,ridge2_predict))

#cross-validation 
r = cross_val_score(ridge1, snp, ics, cv=10) 
r2 = cross_val_score(ridge2, exp, ics, cv = 10) 
print r
print r2
print("Accuracy: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("Accuracy: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 


# #LASSO (on SNP)
# clf = linear_model.Lasso(alpha=1)
# clf_fit = clf.fit(X_train, y_train)
# print clf.score(X_test,y_test)

# #LASSO (on Gene Expression)
# clf2 = linear_model.Lasso(alpha=1)
# clf2_fit = clf2.fit(X_train2, y_train2)
# print clf2.score(X_test2,y_test2)

#support vector regression (on SNP)
svr = SVR()
svr_fit = svr.fit(X_train, y_train)
print svr_fit.score(X_test, y_test)

#support vector regression (on Gene Expression)
svr2 = SVR()
svr2_fit = svr2.fit(X_train2, y_train2)
print svr2_fit.score(X_test2, y_test2)





# # #plot linear regression (doesn't work)
# # # plt.scatter(X_train, y_train,  color='black')
# # # plt.plot(X_train, lin.predict(X_train), color='blue',
# # #          linewidth=3)
# # # plt.xticks(())
# # # plt.yticks(())

# # # plt.show()





<<<<<<< HEAD
#feature selection (does not give better result)
# selector = VarianceThreshold()
# selector.fit(X_train)
# all_var = selector.variances_
# print np.amax(all_var)
# print np.amin(all_var)


# var_d = {}
=======
#feature selection 
"""for SNP"""
selector = VarianceThreshold()
selector.fit(X_train)
all_var = selector.variances_
print np.amax(all_var)
print np.amin(all_var)

var_d = {}
>>>>>>> 5313af1c6755b6e9352d0216e4d5f15f249790e8

# for i in range(0, len(all_var)):
# 	var_d[i] = all_var[i]

<<<<<<< HEAD

# sorted_var_d = sorted(var_d.items(), key=operator.itemgetter(1), reverse=True)

# ind = []
# for i in range(50,5000):
# 	ind.append(sorted_var_d[i][0])
=======
sorted_var_d = sorted(var_d.items(), key=operator.itemgetter(1), reverse=True)


"""for GENE"""
selector2 = VarianceThreshold()
selector2.fit(X_train2)
all_var2 = selector2.variances_
print np.amax(all_var2)
print np.amin(all_var2)

var_d2 = {}

for i in range(0, len(all_var2)):
	var_d2[i] = all_var2[i]

sorted_var_d2 = sorted(var_d2.items(), key=operator.itemgetter(1), reverse=True)




"""Do feature selection for snp. Not taking features with very low variance and very high variance"""

ind = []
for i in range(50,5500):
	ind.append(sorted_var_d[i][0])

# X_train_selected = []
# X_test_selected = []

# tempX = []
# tempXX = []
# for i in range(0,len(ind)):
# 	tempX= X_train[0:314,ind[i]]
# 	tempXX = X_test[0:79,ind[i]]
# 	X_train_selected.append(tempX)
# 	X_test_selected.append(tempXX)


# X_train_selected = np.array(X_train_selected)
# X_test_selected = np.array(X_test_selected)

# X_train_selected = np.transpose(X_train_selected)
# X_test_selected = np.transpose(X_test_selected)

# #linear regression after feature selection
# lin = linear_model.LinearRegression()
# lin_fit = lin.fit(X_train_selected, y_train)

# #calculate R^2 wth test sets
# print lin_fit.score(X_test_selected, y_test) 

# #calculate MSE
# lin_predict = lin_fit.predict(X_test_selected)
# print "mean square error of linear regression (after feature selection) on SNP is " + str(mean_squared_error(y_test, lin_predict)) #y_true, y_predict


tempX = []
tempXX = []
for i in range(0,len(ind)):
	tempX= X_train[0:314,ind[i]]
	tempXX = X_test[0:79,ind[i]]
	X_train_selected.append(tempX)
	X_test_selected.append(tempXX)




"""Do feature selection for Gene. Not taking features with very low variance and very high variance"""
ind2 = []
for i in range(0,6000):
	ind2.append(sorted_var_d2[i][0])

X_train_selected2 = []
X_test_selected2 = []


tempX2 = []
tempXX2 = []
for i in range(0,len(ind2)):
	tempX2= X_train2[0:314,ind2[i]]
	tempXX2 = X_test2[0:79,ind2[i]]
	X_train_selected2.append(tempX2)
	X_test_selected2.append(tempXX2)


"""Now we move on to doing pca for gene"""
pca_gene = PCA(30).fit(X_train2)
components_gene = pca_gene.components_
print (pca_gene.components_).shape
gene_result_train = np.dot(X_train2, components_gene.transpose())
print gene_result_train.shape

gene_result_test = np.dot(X_test2, components_gene.transpose())
print gene_result_test.shape



"""Now we move on to doing pca for SNP"""
pca = PCA(30).fit(X_train)
components = pca.components_
print (pca.components_).shape
snp_result_train = np.dot(X_train, components.transpose())
print snp_result_train.shape

snp_result_test = np.dot(X_test, components.transpose())
print snp_result_test.shape





"""

#linear regression
lin = linear_model.LinearRegression()
lin_fit = lin.fit(snp_result_train, y_train)
lin2 = linear_model.LinearRegression()
lin2_fit = lin2.fit(gene_result_train, y_train2)

#calculate R^2 wth test sets
print "R squared value SNP: " + str(lin_fit.score(snp_result_test, y_test))
print "R squared value GENE: " + str(lin2_fit.score(gene_result_test, y_test2)) 

#calculate MSE
lin_predict = lin_fit.predict(snp_result_test)
lin2_predict = lin2_fit.predict(gene_result_test)
print "mean square error of linear regression on SNP is " + str(mean_squared_error(y_test, lin_predict)) #y_true, y_predict
print "mean square error of linear regression on gene expression is " + str(mean_squared_error(y_test2,lin2_predict))

#cross-validation
lins = cross_val_score(lin, snp, ics, cv=10)  
lins2 = cross_val_score(lin2, exp, ics, cv = 10) 
print lins
print lins2
print("Accuracy: %0.2f (+/- %0.2f)" % (lins.mean(), lins.std() * 2))
print("Accuracy: %0.2f (+/- %0.2f)" % (lins2.mean(), lins2.std() * 2))

"""





#ridge regression
ridge1 = linear_model.Ridge()
ridge1_fit = ridge1.fit(snp_result_train, y_train)
ridge2 = linear_model.Ridge()
ridge2_fit = ridge2.fit(gene_result_train, y_train2)

#calculate R^2 with test sets
print "R squared value SNP: " + str(ridge1.score(snp_result_test, y_test)) 
print "R squared value GENE: " + str(ridge2.score(gene_result_test, y_test2)) 

#calculate MSE
ridge_predict = ridge1_fit.predict(snp_result_test)
ridge2_predict = ridge2_fit.predict(gene_result_test)


print 'length of ridge predict is ' + str(len(ridge_predict))
print "mean square error of ridge regression on SNP is " + str(mean_squared_error(y_test, ridge_predict)) #y_true, y_predict
print "mean square error of ridge regression on gene expression is " + str(mean_squared_error(y_test2,ridge2_predict))

#cross-validation 
r = cross_val_score(ridge1, snp, ics, cv=10) 
r2 = cross_val_score(ridge2, exp, ics, cv = 10) 
print r
print r2
print("Accuracy: %0.2f (+/- %0.2f)" % (r.mean(), r.std() * 2))
print("Accuracy: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std() * 2)) 

>>>>>>> 5313af1c6755b6e9352d0216e4d5f15f249790e8







"""

X_train_selected2 = np.array(X_train_selected2)
X_test_selected2 = np.array(X_test_selected2)

X_train_selected2 = np.transpose(X_train_selected2)
X_test_selected2 = np.transpose(X_test_selected2)

#linear regression after feature selection
lin2 = linear_model.LinearRegression()
lin_fit2 = lin2.fit(X_train_selected2, y_train2)

#calculate R^2 wth test sets
print lin_fit2.score(X_test_selected2, y_test2) 

#calculate MSE
lin_predict2 = lin_fit2.predict(X_test_selected2)
print "mean square error of linear regression (after feature selection) on gene is " + str(mean_squared_error(y_test2, lin_predict2)) #y_true, y_predict
"""



# #Construct a dictionary (key: value) = (index: variance)
# #sort the dictionary's values in ascending order 
# #get the first 500 key-val pairs


# #variance over features for SNP
# sorted_all_var = sorted(all_var, reverse=True)
# first500 = sorted_all_var[0:500]
# y_pos = np.arange(22419)
# plt.bar(y_pos,sorted_all_var,align='center')
# # plt.xticks(y_pos, objects)
# plt.ylabel('Variances')
# plt.xlabel('Features')
# plt.title('Variance over Features on SNP Data')
# plt.savefig('snp_var.png')



#variance over features for Gene Expression

#feature selection (does not give better result)
# selector = VarianceThreshold()
# selector.fit(X_train2)
# all_var = selector.variances_
# sorted_all_var = sorted(all_var, reverse=True)
# y_pos = np.arange(18926)
# plt.bar(y_pos,sorted_all_var,align='center')
# # plt.xticks(y_pos, objects)
# plt.ylabel('Variances')
# plt.xlabel('Features')
# plt.title('Variance over Features on Gene Expression Data')
# plt.savefig('gene_var.png')


















# #support vector machine (Doesn't work bc this project is not a classication problem)
# # svm1 = SVC(kernel="linear")
# # svm1_fit = svm1.fit(snp, ics)
# svm2 = SVC(kernel='linear')
# svm2_fit = svm2.fit(X_train2, y_train2)

# #mean accuracy on the given test data and labels
# # print svm1_fit.score(X_test, y_test)
# # print svm1.score(X_test, y_test)

# print sv2.score(X_test2, y_test2)
# print 'end of file'


# print 'dimension of x train is '
# print X_train.shape
# pca = PCA().fit(X_train)  
# print (pca.components_).shape


