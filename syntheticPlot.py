import numpy as np
import matplotlib.pyplot as plt
from analyze_synthetic import getValidMappings 
print("Started")
markers = ["*", "o", "d", "1", "s", "p", "D"]

def getScores(directory):
    epsilons = np.arange(0.1, 1, 0.1).tolist()
    epsilons = ["%.1f" % e for e in epsilons]
    TICC = [[],[],[]]
    CASC0_6 = [[],[],[]]
    CASC0_8 = [[],[],[]]
    KMEANS = [[],[],[]]
    GMM = [[],[],[]]
    HMM = [[],[],[]]

    correctname = "%s/correct.out" % directory
    for e in epsilons:
        print(e)
        ticcscores = getValidMappings(correctname,  "%s/%s/old/assign.out" % (directory, e) )
        casc06scores =  getValidMappings(correctname,  "%s/%s/0.6/assign.out" % (directory, e) )
        casc08scores =  getValidMappings(correctname,  "%s/%s/0.8/assign.out" % (directory, e) )
        kmeans =  getValidMappings(correctname,  "%s/%s/kmeans.out" % (directory, e) )
        gmm =  getValidMappings(correctname,  "%s/%s/gmm.out" % (directory, e) )
        hmm =  getValidMappings(correctname,  "%s/%s/hmm.out" % (directory, e) )

        for i in range(len(ticcscores)):
            TICC[i].append(ticcscores[i])
            CASC0_6[i].append(casc06scores[i])
            CASC0_8[i].append(casc08scores[i])
            KMEANS[i].append(kmeans[i])
            GMM[i].append(gmm[i])
            HMM[i].append(hmm[i])
    return TICC, CASC0_6, CASC0_8, KMEANS, GMM, HMM

def saveScores(directory, output):
    TICC, CASC0_6, CASC0_8, KMEANS, GMM, HMM = getScores(directory)
    np.savetxt("%s/ticc" % output, np.array(TICC), delimiter=",")
    np.savetxt("%s/casc0_6" % output, np.array(CASC0_6), delimiter=",")
    np.savetxt("%s/casc0_8" % output, np.array(CASC0_8), delimiter=",")
    np.savetxt("%s/kmeans"% output, np.array(KMEANS), delimiter=",")
    np.savetxt("%s/gmm" % output, np.array(GMM), delimiter=",")
    np.savetxt("%s/hmm" % output, np.array(HMM), delimiter=",")

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


# saveScores("ordered_synthetic", "final_scores")
# loadAndPlotScores("final_scores_perturbed_0.7_single_seq/macro_final_scores")
# loadAndPlotScores("final_scores_perturbed_0.7_single_seq/micro_final_scores")
loadAndPlotScores("final_scores_perturbed_0.7_single_seq/weighted_final_scores")
# loadAndPlotScores("final_scores")
# loadAndPlotScores("old_final_scores")


