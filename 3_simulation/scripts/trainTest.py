import sys, getopt
import numpy as np
import csv
from sklearn.preprocessing import Imputer
from sklearn import preprocessing
from sklearn.model_selection import LeaveOneGroupOut
from sklearn import svm, datasets, ensemble
from sklearn.utils import shuffle
from scipy import interp
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import gzip
import random

random.seed(21561234)



argv = sys.argv[1:]

trainfile = ''
testfile = ''
probabilityfile= ''
repetitions=1
try:
	opts, args = getopt.getopt(argv,"h::",["train=","test=","prediction=","repetitions="])
except getopt.GetoptError:
	print('jsonToTable.py --train <train-csv-file> --test <test-csv-file> --prediction <p file> --repetitions 100')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('jsonToTable.py --train <train-csv-file> --test <test-csv-file> --prediction <p file> --repetitions 100')
		sys.exit()
	elif opt in ("--train"):
		trainfile = arg
	elif opt in ("--test"):
		testfile = arg
	elif opt in ("--prediction"):
		probabilityfile = arg
	elif opt in ("--repetitions"):
			repetitions = int(arg)
	else:
		print('jsonToTable.py --train <train-csv-file> --test <test-csv-file> --prediction <p file> --repetitions 100')
		sys.exit(2)
print('Train is ',trainfile)
print('Test is ',testfile)
print('Probability file is', probabilityfile)
print('Repetitions are ', repetitions)


def loadData(file, remove=set()):

	gene_ids = set()

	X = []
	y = []

	with open(file) as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			label = int(row["label"])
			if (label == 1):
				gene_id = int(row['gene_id'])
				gene_ids.add(gene_id)
				if gene_id in remove:
					continue
				else:
					i = i+1
			X.append([row["feature_score"], row["cadd_phred_score"], row["combined_score"], row["cadd_raw_score"], row["gestalt_score"], row["boqa_score"], row["pheno_score"]])
			y.append(label)

		print(i)

	y = np.array(y)
	X = np.array(X)

	return X,y,gene_ids

X_test, y_test, remove = loadData(testfile)
X_train, y_train, gene_ids = loadData(trainfile, remove)

for i in range(X_train.shape[1]):
	m = min(X_train[:,i])
	X_train[X_train[:,i]=='nan',i]=m
	X_test[X_test[:,i]=='nan',i]=m

X_train = X_train.astype(float)
X_test = X_test.astype(float)

probabilities = []
y = []

for i in range(repetitions):
	np.random.seed(random.randint(0, 4294967295))

	X_train_rnd, y_train_rnd = shuffle(X_train, y_train, random_state=random.randint(0, 4294967295))
	X_test_rnd, y_test_rnd = shuffle(X_test, y_test, random_state=random.randint(0, 4294967295))

	#classifier = svm.SVC(kernel='poly', probability=True, class_weight='balanced')
	classifier = ensemble.RandomForestClassifier(n_estimators = 100,max_features=3,n_jobs=2)

	'''
	imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
	imp.fit(X_train)
	X_train_transformed = imp.transform(X_train)
	'''
	X_train_transformed = X_train_rnd
	normalizer = preprocessing.Normalizer()
	normalizer.fit(X_train_transformed)
	X_train_normalized = normalizer.transform(X_train_transformed)

	classifier.fit(X_train_normalized, y_train_rnd)

	#X_test_transformed = imp.transform(X_test)
	X_test_transformed = X_test_rnd

	probabilities.extend(classifier.predict_proba(normalizer.transform(X_test_transformed))[:, 1])
	y.extend(y_test_rnd)

if probabilityfile != '':
	with gzip.open(probabilityfile, 'wb') as f:
		for i, prob in enumerate(probabilities):
			if i == 0:
				toWrite = str(y[i]) + "\t" + str(prob)
			else:
				toWrite = "\n" + str(y[i]) + "\t" + str(prob)
			f.write(toWrite.encode())

fpr, tpr, thresholds = roc_curve(y, probabilities)
roc_auc = auc(fpr, tpr)
print("AUROC: "+ str(roc_auc))

plt.plot(fpr, tpr, linewidth=2, color='darkorange', label='ROC fold (area = %0.2f)' % ( roc_auc))

precision, recall, thresholds = precision_recall_curve(y, probabilities)
prc_auc = auc(recall,precision)
print("AUPRC: "+ str(prc_auc))
