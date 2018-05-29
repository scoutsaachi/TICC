import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from analyze_synthetic import getValidMappings 
print("Started")
markers = ["*", "o", "d", "1", "s", "p", "D"]
EPSILONS = np.arange(0, 0.8, 0.1).tolist()

mapping = {
    'ticc': 'old/assign.out',
    'casc06': '0.6/assign.out',
    'casc08': '0.8/assign.out',
    'casc99': '0.99/assign.out',
    'kmeans': 'kmeans.out',
    'gmm': 'gmm.out',
    'hmm': 'hmm.out'
}

def getScores(directory):
    epsilons = EPSILONS
    epsilons = ["%.1f" % e for e in epsilons]
    TICC = [[],[],[]]
    CASC0_6 = [[],[],[]]
    CASC0_8 = [[],[],[]]
    KMEANS = [[],[],[]]
    GMM = [[],[],[]]
    HMM = [[],[],[]]

    all_scores = {}

    correctname = "%s/correct.out" % directory
    for k, v in mapping.items():
        scores = [[],[],[]]
        for e in epsilons:
            s = getValidMappings(correctname,  "%s/%s/%s" % (directory, e, v))
            for i in range(len(s)): scores[i].append(s[i])
        all_scores[k] = scores
    return all_scores

def saveScores(directory, output):
    all_scores = getScores(directory)
    for k,v in mapping.items():
        np.savetxt("%s/%s" % (output, k), np.array(all_scores[k]), delimiter=',')

def loadAndPlotScores(directory):
    TICC = np.loadtxt("%s/ticc" % directory, delimiter=",")
    CASC1 = np.loadtxt("%s/casc0_6" % directory, delimiter=",")
    CASC2 = np.loadtxt("%s/casc0_8" % directory, delimiter=",")
    KMEANS = np.loadtxt("%s/kmeans"% directory, delimiter=",")
    GMM = np.loadtxt("%s/gmm" % directory, delimiter=",")
    HMM = np.loadtxt("%s/hmm" % directory, delimiter=",")

    mapping = {
        "TICC": TICC,
        "CASC, gamma=0.6": CASC1,
        "CASC, gamma=0.8": CASC2,
        "kmeans": KMEANS,
        "gmm": GMM,
        # "hmm": HMM
    }
    endPoint = 0.7
    labels = np.arange(0.1, endPoint, 0.1).tolist()
    print("read")
    titles = ["Weighted Macro F1 Score on Motif Clusters", 
        "Weighted Macro F1 Score on All Clusters",
        "Accuracy Score on Motif Segments"]
    for k,v in mapping.items():
        stringVals = ["%.4f" % val for val in v[0][:len(labels)]]
        print("%s, %s" % (k, ' & '.join(stringVals)))
    for i in range(3):
        plt.figure(i, figsize=(6, 7))
        ax = plt.subplot(111)
        markerCount = 0
        for k,v in mapping.items():
            plt.plot(labels, v[i][:len(labels)], label=k, marker=markers[markerCount])
            markerCount += 1
        box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
        #         fancybox=True, shadow=True, ncol=1)
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


saveScores("ordered_synthetic", "final_scores_newset")
# loadAndPlotScores("final_scores_perturbed_0.7_single_seq/macro_final_scores")
# loadAndPlotScores("final_scores_perturbed_0.7_single_seq/micro_final_scores")
#loadAndPlotScores("final_scores_perturbed_0.7_single_seq/weighted_final_scores")
# loadAndPlotScores("final_scores")
# loadAndPlotScores("old_final_scores")


