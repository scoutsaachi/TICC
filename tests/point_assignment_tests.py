import sys 
sys.path.append('..')
from ticc.motif import *
from ticc.test_utils import *
import matplotlib.pyplot as plt

def scoreNormalizedScore(motifHMM, nullHMM, t):
    return motifHMM.getEndingScore(t)/t

def scoreNormalizedSqrt(motifHMM, nullHMM, t):
    return motifHMM.getEndingScore(t)/np.sqrt(t)

def scoreMotifNull(motifHMM, nullHMM, t):
    return motifHMM.getEndingScore(t)/nullHMM.getMostLikelyScore(t)

def scoreNormMotifNull(motifHMM, nullHMM, t):
    return scoreMotifNull(motifHMM, nullHMM, t)/t

scoringDict = [
    ("LL/n", scoreNormalizedScore),
    ("LL/sqrt(n)", scoreNormalizedSqrt),
    ("LL_Motif/LL_Null", scoreMotifNull)
]

def test(breakpoints, confidenceMean, confidenceVariance, numClusters, T):
    ll = GenerateFakeData(breakpoints, confidenceMean, confidenceVariance, numClusters)
    motif_hmm = MotifHMM(ll, np.array([0,1,2]),0)
    null_hmm = NullHMM(ll, 0)
    scores = []
    for i in range(len(scoringDict)):
        scores.append([])
    for i in range(T):
        motif_hmm.UpdateStep(i)
        null_hmm.UpdateStep(i)
        for j in range(len(scoringDict)):
            _, fn = scoringDict[j]
            score = fn(motif_hmm, null_hmm, i)
            scores[j].append(score)
    return scores

def plotSequence(scores, labels, breakpoints, title, idx):
    '''
    scores is a list of lists (which are the lines to plot)
    '''
    x = np.arange(len(labels))
    currLine = 0
    for i in breakpoints:
        currLine += i
        plt.axvline(x=currLine, linestyle=':')
    plt.title(title)
    plt.xticks(x, labels)
    for s in scores:
        plt.plot(x, s, '-o')
    plt.show()

def GenerateLabels(breakpoints):
    '''
    Create a labels list
    breakpoints - a list of the length of each cluster, with the first and last values being random
    '''
    assert len(breakpoints) > 2
    labels = ['_']*breakpoints[0]
    for i in range(1, len(breakpoints) - 1):
        identifier = i-1
        label = chr(ord('A') + identifier)
        labels += [label]*breakpoints[i]
    labels += ['_']*breakpoints[1]
    return labels
    
def testMultiple(iterations):
    # input params
    breakpoints = [50, 10, 10, 10, 10] # A, B, C
    confidenceMean = 0.5
    confidenceVariance = 0.1
    numClusters = 4 # A, B, C, D

    bpStr = ','.join([str(i) for i in breakpoints])
    toStr = "[%s], conf mean: %s, conf var: %s" % (bpStr, confidenceMean, confidenceVariance)
    labels = GenerateLabels(breakpoints)

    results = []
    for i in range(len(scoringDict)):
        results.append([])
    for i in range(iterations):
        scores = test(breakpoints, confidenceMean, confidenceVariance, numClusters, len(labels))
        for j, s in enumerate(scores):
            results[j].append(s)

    for i in range(len(scoringDict)):
        fnName = scoringDict[i][0]
        title = "%s: %s" % (fnName, toStr)
        plotSequence(results[i], labels, breakpoints, title, i+1)
    plt.show()
testMultiple(5)