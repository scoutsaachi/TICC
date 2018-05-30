import sys
import numpy as np
from sklearn.metrics import f1_score, accuracy_score
from generateDatasets.constants import NUM_CLUSTERS, GARBAGE_CLUSTERS, NUM_SEQS, NUM_GARBAGE, LEN_SEGMENT, CLUSTER_SEQUENCE
import pickle


ONLY_PICK_BEST = False

def analyzeMotifs(motifFname):
    true_result = []
    for _ in range(NUM_SEQS):
        true_result += [0 for _ in range(NUM_GARBAGE*LEN_SEGMENT)]
        true_result += [1 for _ in range(len(CLUSTER_SEQUENCE)*LEN_SEGMENT)]
    # true_result is the true motif assign vals
    
    f = open(motifFname, "rb")
    motifs = pickle.load(f) 
    motifResult = [0 for _ in range(len(true_result))]
    if ONLY_PICK_BEST:
        lengths = sorted([(len(v),k) for k,v in motifs.items()], reverse=True)
        top = lengths[0][1]
        print(lengths[0])
        motifs  = {top:motifs[top]}
    for k,v in motifs.items():
        for s, e in v:
            motifResult[s:e+1] = [1 for _ in range(s, e+1)]
    allResult = [1 for _ in range(len(true_result))]
    scores = [f1_score(true_result, motifResult), f1_score(true_result, allResult), accuracy_score(true_result, motifResult), accuracy_score(true_result, allResult)]
    print("F1 Score: %s, Naive Score: %s, Accuracy Score: %s, Naive Accuracy score: %s" % tuple(scores))
    return scores

EPSILONS = np.arange(0, 0.8, 0.1).tolist()

def getScores(directory):
    epsilons = EPSILONS
    epsilons = ["%.1f" % e for e in epsilons]
    all_scores = [[],[],[]]
    for e in epsilons:
        print(e)
        s = analyzeMotifs("%s/%s/motifs.pkl" % (directory, e))
        for i in range(len(s)): all_scores[i].append(s[i])
    return all_scores

if __name__ == "__main__":
    fileDir = sys.argv[1]
    analyzeMotifs(fileDir)
