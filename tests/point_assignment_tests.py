import sys 
sys.path.append('..')
from ticc.motif import *
from ticc.test_utils import *
import matplotlib.pyplot as plt

def test():
    breakpoints = [100, 10, 10, 10, 10] # A, B, C
    numClusters = 4 # A, B, C, D
    ll, labels = GenerateFakeData(breakpoints, 0.75, 0.05, numClusters)
    motif_hmm = MotifHMM(ll, np.array([0,1,2]),0)
    null_hmm = NullHMM(ll, 0)
    scores = []
    for i in range(len(labels)):
        motif_hmm.UpdateStep(i)
        null_hmm.UpdateStep(i)
        # score generation
        score = motif_hmm.getEndingScore(i)/i
        scores.append(score)
    plotSequence(scores, labels, breakpoints)

def plotSequence(scores, labels, breakpoints):
    x = np.arange(len(labels))
    currLine = 0
    for i in breakpoints:
        currLine += i
        plt.axvline(x=currLine, linestyle=':')
    plt.xticks(x, labels)
    plt.plot(x, scores, '-o')
    plt.show()

    
test()