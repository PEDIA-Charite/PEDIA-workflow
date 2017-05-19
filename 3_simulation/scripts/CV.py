import sys, getopt
import numpy as np
import csv
from sklearn.preprocessing import Imputer
from sklearn import preprocessing
from sklearn.model_selection import LeaveOneGroupOut
from sklearn import svm, datasets, ensemble
from scipy import interp
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import gzip

np.random.seed(42)

argv = sys.argv[1:]

inputfile = ''
probabilityfile= ''
try:
	opts, args = getopt.getopt(argv,"hi:p:",["ifile=","pfile="])
except getopt.GetoptError:
	print('jsonToTable.py -i <input-csv-file> -p <output probability file>')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('jsonToTable.py -i input-csv-file> -p <output probability file>')
		sys.exit()
	elif opt in ("-i", "--ifolder"):
		inputfile = arg
	elif opt in ("-p", "--pfile"):
		probabilityfile = arg
	else:
		print('jsonToTable.py -i <input-csv-file> -p <output probability file>')
		sys.exit(2)
print('Input csv file is ', inputfile)
print('Probability file is', probabilityfile)


#random_state = np.random.RandomState(0)

data = []
y = []
grouping = {}
groups = []

with open(inputfile) as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		data.append([row["feature_score"], row["cadd_phred_score"], row["combined_score"], row["cadd_raw_score"], row["gestalt_score"], row["boqa_score"], row["pheno_score"]])
		label = int(row["label"])
		y.append(label)
		groups.append(row["case"])
		if (label==1):
			grouping[row["case"]] = row["gene_id"]

y = np.array(y)
data = np.array(data)

'''
	for row in data:
		clean = [float(i) for i in row]
		clean_data.append(clean)
'''

for i, case in enumerate(groups):
	if case in grouping:
		groups[i] = grouping[case]
	else:
		groups[i] = "nan"
#classifier = svm.SVC(kernel='poly', probability=True, class_weight='balanced')
classifier = ensemble.RandomForestClassifier(class_weight='balanced')

logo = LeaveOneGroupOut()

probabilities = []
labels = []

r = 0
for train, test in logo.split(data, y, groups=groups):
	r+=1
	print("Train/Test round "+ str(r) + "/"+ str(len(set(groups))))
	X_train, X_test, y_train, y_test = data[train], data[test], y[train], y[test]

	for i in range(X_train.shape[1]):
		m = min(X_train[:,i])
		X_train[X_train[:,i]=='nan',i]=m
		X_test[X_test[:,i]=='nan',i]=m

	X_train = X_train.astype(float)
	X_test = X_test.astype(float)
'''
	#### impute median
	imp = Imputer(missing_values='NaN', strategy='median', axis=0)
	imp.fit(X_train)
	X_train_transformed = imp.transform(X_train)
'''
	X_train_transformed = X_train
	####
	#### normalize
	normalizer = preprocessing.Normalizer()
	normalizer.fit(X_train_transformed)
	X_train_normalized = normalizer.transform(X_train_transformed)
	'''
	X_train_normalized = X_train

	X_test_transformed  = imp.transform(X_test)
	'''
	X_test_transformed = X_test
	X_test_normalized = normalizer.transform(X_test_transformed)

	X_test_normalized = X_test
	probabilities.extend(classifier.fit(X_train_normalized, y_train).predict_proba(X_test_normalized)[:, 1])
	labels.extend(y_test)



fpr, tpr, thresholds = roc_curve(labels, probabilities)
roc_auc = auc(fpr, tpr)
print("AUROC: "+ str(roc_auc))

plt.plot(fpr, tpr, linewidth=2, color='darkorange', label='ROC fold (area = %0.2f)' % ( roc_auc))

precision, recall, thresholds = precision_recall_curve(labels, probabilities)
prc_auc = auc(recall,precision)
print("AUPRC: "+ str(prc_auc))

if probabilityfile != '':
	with gzip.open(probabilityfile, 'wb') as f:
		for i, prob in enumerate(probabilities):
			if i == 0:
				toWrite = str(labels[i]) + "\t" + str(prob)
			else:
				toWrite = "\n" + str(labels[i]) + "\t" + str(prob)
			f.write(toWrite.encode())
