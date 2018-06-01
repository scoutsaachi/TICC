import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from analyze_synthetic import getValidMappings 
print("Started")
EPSILONS = np.arange(0, 0.8, 0.1).tolist()

mapping = {
    'casc06': ('0.6/assign.out', "CASC $\gamma$ = 0.6", "-ro"),
    'casc99': ('0.99/assign.out', "CASC, $\gamma$ = 0.99", "-ms"),
    'ticc': ('old/assign.out', "TICC", "-.bd"),
    'kmeans': ('kmeans.out', "KMEANS", ":k*"),
    'gmm': ('gmm.out', "GMM", ":g1"),
    'hmm': ('hmm.out', "HMM", ":c2"),
}

def getScores(directory):
    epsilons = EPSILONS
    epsilons = ["%.1f" % e for e in epsilons]
    all_scores = {}
    correctname = "%s/correct.out" % directory
    for k, v in mapping.items():
        print(k)
        scores = [[],[],[]]
        for e in epsilons:
            print(e)
            s = getValidMappings(correctname,  "%s/%s/%s" % (directory, e, v[0]))
            for i in range(len(s)): scores[i].append(s[i])
        all_scores[k] = scores
    return all_scores

def saveScores(directory, output):
    all_scores = getScores(directory)
    for k,v in mapping.items():
        np.savetxt("%s/%s" % (output, k), np.array(all_scores[k]), delimiter=',')

def loadAndPlotScores(directory):
    scoreValues = {}
    for k,v in mapping.items():
        scoreValues[k] = np.loadtxt("%s/%s" % (directory, k), delimiter=",")
    labels = EPSILONS
    titles = ["Weighted Macro F1 Score on Motif Clusters", 
        "Weighted Macro F1 Score on All Clusters",
        "Accuracy Score on Motif Segments"]
    for i in range(3):
        plt.figure(i, figsize=(6, 7))
        ax = plt.subplot(111)
        markerCount = 0
        for k,v in scoreValues.items():
            plt.plot(labels, v[i][:len(labels)], mapping[k][2], label=mapping[k][1])
            markerCount += 1
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
          fancybox=True, shadow=True, ncol=3)
        plt.xlabel("Fraction of noisy segments")
        plt.ylim(ymin=0)
        plt.title(titles[i])
        if i == 2:
            plt.ylabel("Accuracy score on motif segments vs Noise")
        else:
            plt.ylabel("Weighted Macro F1 Score")
        # plt.legend()
    plt.show()

saveScores("ordered_synthetic_retry2", "final_scores_retry2")
#loadAndPlotScores("final_scores_retry2")


