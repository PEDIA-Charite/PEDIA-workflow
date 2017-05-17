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
from sklearn.metrics import roc_curve, auc


argv = sys.argv[1:]
	
inputfile = ''
try:
	opts, args = getopt.getopt(argv,"h:i:",["ifile="])
except getopt.GetoptError:
	print('jsonToTable.py -i <input-csv-file>')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('jsonToTable.py -i input-csv-file>')
		sys.exit()
	elif opt in ("-i", "--ifolder"):
		inputfile = arg
	else:
		print('jsonToTable.py -i <input-csv-file>')
		sys.exit(2)
print('Input folder is ',inputfile)


random_state = np.random.RandomState(0)

data = []
y = []
groups = []

with open(inputfile) as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		data.append([row["feature_score"], row["cadd_phred_score"], row["combined_score"], row["cadd_raw_score"], row["gestalt_score"], row["boqa_score"], row["pheno_score"]])
		y.append(int(row["label"]))
		groups.append(row["case"])

y = np.array(y)
data = np.array(data)

#classifier = svm.SVC(kernel='poly', probability=True, class_weight='balanced', random_state=random_state)
classifier = ensemble.RandomForestClassifier(class_weight='balanced')


logo = LeaveOneGroupOut()

probabilities = []
labels = []

i = 0
for train, test in logo.split(data, y, groups=groups):
	i+=1
	print("Train/Test round "+ str(i) + "/"+ str(len(set(groups))))
	X_train, X_test, y_train, y_test = data[train], data[test], y[train], y[test]

	imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
	imp.fit(X_train)
	X_train_transformed = imp.transform(X_train)
	
	normalizer = preprocessing.Normalizer()
	normalizer.fit(X_train_transformed)
	X_train_normalized = normalizer.transform(X_train_transformed)

	probabilities.extend(classifier.fit(X_train_normalized, y_train).predict_proba(normalizer.transform(imp.transform(X_test)))[:, 1])
	labels.extend(y_test)

fpr, tpr, thresholds = roc_curve(labels, probabilities)
roc_auc = auc(fpr, tpr)
print(roc_auc)
plt.plot(fpr, tpr, linewidth=2, color='darkorange', label='ROC fold (area = %0.2f)' % ( roc_auc))

