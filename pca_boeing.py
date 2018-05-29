from sklearn.decomposition import PCA
import numpy as np
import csv

path = "/dfs/scratch0/dataset/20161030-Boeing/data/WD673_rt1.txt"
#path="boeing.out"

arr = []
with open(path, 'r') as f:
	r = csv.reader(f, delimiter='\t', skipinitialspace=False)
	for i, row in enumerate(r):
		if i==0: continue
		if i % 1000 == 0: print(i)
		arr += [[float(v if v!= '' else 0) for v in row]]

#inputdata = np.loadtxt(path, skiprows=1)
arr = np.array(arr)
dims = 13
pca = PCA(n_components=dims)
inputdata = pca.fit_transform(arr)
np.savetxt("boeing/data.out", inputdata, delimiter=",")
