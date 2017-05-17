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


argv = sys.argv[1:]
	
trainfile = ''
testfile = ''
try:
	opts, args = getopt.getopt(argv,"h::",["train=","test="])
except getopt.GetoptError:
	print('jsonToTable.py --train <train-csv-file> --test <test-csv-file>')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('jsonToTable.py --train <train-csv-file> --test <test-csv-file>')
		sys.exit()
	elif opt in ("--train"):
		trainfile = arg
	elif opt in ("--test"):
		testfile = arg
	else:
		print('jsonToTable.py --train <train-csv-file> --test <test-csv-file>')
		sys.exit(2)
print('Train is ',trainfile)
print('Test is ',testfile)


def loadData(file):

	random_state = np.random.RandomState(0)

	X = []
	y = []

	with open(trainfile) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			X.append([row["feature_score"], row["cadd_phred_score"], row["combined_score"], row["cadd_raw_score"], row["gestalt_score"], row["boqa_score"], row["pheno_score"]])
			y.append(int(row["label"]))

	y = np.array(y)
	X = np.array(X)

	return X,y

X_train, y_train = loadData(trainfile)
X_test, y_test = loadData(testfile)

#classifier = svm.SVC(kernel='poly', probability=True, class_weight='balanced', random_state=random_state)
classifier = ensemble.RandomForestClassifier(class_weight='balanced')


imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
imp.fit(X_train)
X_train_transformed = imp.transform(X_train)
	
normalizer = preprocessing.Normalizer()
normalizer.fit(X_train_transformed)
X_train_normalized = normalizer.transform(X_train_transformed)

classifier.fit(X_train_normalized, y_train)

probabilities = classifier.predict_proba(normalizer.transform(imp.transform(X_test)))[:, 1]

fpr, tpr, thresholds = roc_curve(y_test, probabilities)
roc_auc = auc(fpr, tpr)
print("AUROC: "+ str(roc_auc))

plt.plot(fpr, tpr, linewidth=2, color='darkorange', label='ROC fold (area = %0.2f)' % ( roc_auc))

precision, recall, thresholds = precision_recall_curve(y_test, probabilities)
prc_auc = auc(recall,precision)
print("AUPRC: "+ str(prc_auc))



