from sklearn import datasets
from sklearn import feature_extraction
from sklearn import feature_selection
from sklearn.model_selection import train_test_split
from sklearn import tree
from sklearn import metrics
from sklearn import model_selection
import pandas as pd
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

# def one_hot_dataframe(data, cols, replace=False):
#     vec = feature_extraction.DictVectorizer()
#     mkdict = lambda row: dict((col, row[col]) for col in cols)
#     vecData = pd.DataFrame(vec.fit_transform(
#     data[cols].apply(mkdict, axis=1)).toarray())
#     vecData.columns = vec.get_feature_names()
#     vecData.index = data.index
#     if replace:
#         data = data.drop(cols, axis=1)
#         data = data.join(vecData)
#     return (data, vecData)


cd_path = './collapsed_data2.tsv'
cd = pd.read_csv(cd_path, sep='\t')
cd = cd.T
cd.columns = cd.iloc[0]
cd = cd.drop(cd.index[0])
cd['Cohort'] = ""
for index, row in cd.iterrows():
    stat = str(index).split('_')[0]
    if stat == 'IBD' or stat == 'IR':
        coh = 0
    elif stat == 'NDAR':
        coh = 1
    # else:
    #     coh = 2
    row['Cohort'] = coh


# # print(cd.describe())
# cd, cd_n = one_hot_dataframe(cd, ['Cohort'], replace=True)
# print(cd)

cd.fillna(0, inplace=True)
cd_target = cd['Cohort']
cd_data = cd.drop(['Cohort'], axis=1)
X_train, X_test, y_train, y_test =train_test_split(cd_data, cd_target, test_size=0.25, random_state=33)
# print(X_train)

dt = tree.DecisionTreeClassifier(criterion='entropy')
dt = dt.fit(X_train, y_train)
# y_pred = dt.predict(X_test)
# print("Accuracy:{0:.3f}".format(metrics.accuracy_score(y_test, y_pred)), "\n")
#
fs = feature_selection.SelectPercentile(feature_selection.chi2, percentile=30)
x_train_fs = fs.fit_transform(X_train, y_train)
dt = dt.fit(x_train_fs, y_train)
X_test_fs = fs.transform(X_test)
y_pred_fs = dt.predict(X_test_fs)
print("Accuracy:{0:.3f}".format(metrics.accuracy_score(y_test, y_pred_fs)), "\n")

# percentiles = range(1, 100, 5)
# results = []
# acclist = []
# for i in range(1, 100, 5):
#     fs = feature_selection.SelectPercentile(feature_selection.chi2, percentile=i)
#     X_train_fs = fs.fit_transform(X_train, y_train)
#     X_test_fs = fs.transform(X_test)
#     scores = model_selection.cross_val_score(dt, X_train_fs, y_train, cv=5)
#     results = np.append(results, scores.mean())
#     dta = dt.fit(X_train_fs, y_train)
#     y_pred_fs = dta.predict(X_test_fs)
#     # print(y_pred_fs)
#     acc = metrics.accuracy_score(y_test, y_pred_fs)
#     acclist = np.append(acclist, acc)
#     # print(score, acc)

# optimal_percentil = np.where(results == results.max())[0]
# print("Optimal number of features:{0}".format(percentiles[optimal_percentil]), "\n")


# plt.figure(1)
# plt.plot(percentiles, results)
# plt.savefig('./feat_vs_acc_plot_trainset.png')
#
# plt.figure(2)
# plt.plot(percentiles, acclist)
# plt.savefig('./feat_vs_acc_plot_testset.png')


# plottrain = pl.figure()
# plottrain.xlabel("Number of features selected")
# plottrain.ylabel("Cross-validation accuracy)")
# plottrain.plot(percentiles, results)
# plottrain.savefig('./feat_vs_acc_plot_trainset.png')
#
# plottest = pl.figure()
# plottest.xlabel("Number of features selected")
# plottest.ylabel("Cross-validation accuracy)")
# plottest.plot(percentiles, acclist)
# plottest.savefig('./feat_vs_acc_plot_testset.png')