import numpy as np 
from sklearn.cluster import KMeans
from hmmlearn.hmm import GaussianHMM
from sklearn.mixture import GaussianMixture
import sys

def getData(data_file):
    data = np.loadtxt(data_file, delimiter=",")
    return data

def performBaseline(data_file, out_file, baselineType="KMEANS", K=10):
    data = np.loadtxt(data_file, delimiter=",")
    labels=None
    if baselineType == "KMEANS":
        kmeans = KMeans(n_clusters=K, random_state=0, n_jobs=30)
        kmeans.fit(data)
        labels = kmeans.labels_
    elif baselineType == "GMM":
        gmm = GaussianMixture(n_components=K, covariance_type="full")
        gmm.fit(data)
        labels = gmm.predict(data)
    elif baselineType == "HMM":
        hmm = GaussianHMM(n_components=K, covariance_type="full", n_iter=100)
        hmm.fit(data)
        labels = hmm.predict(data)
    else:
        assert False

    labels = np.array(labels)
    np.savetxt(out_file, labels, delimiter=",", fmt='%d')

if __name__ == "__main__":
    assert len(sys.argv) > 3
    baseline = sys.argv[1]
    infile = sys.argv[2]
    outfile = sys.argv[3]
    performBaseline(infile, outfile, baseline)
