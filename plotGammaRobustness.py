import numpy as np
import matplotlib
# matplotlib.use('agg')
import matplotlib.pyplot as plt
from analyze_synthetic import getValidMappings 
print("Started")
GAMMAS = ["0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "0.99"]

def getScores(directory):
    correctname = "%s/correct.out" % directory
    scores = [[],[],[]]
    for g in GAMMAS:
        print(g)
        s = getValidMappings(correctname,  "%s/%s" % (directory, g))
        for i in range(len(s)): scores[i].append(s[i])
    return scores

def saveScores(directory, output):
    all_scores = getScores(directory)
    np.savetxt("%s/%s" % (output, "gammaRobust.out"), np.array(all_scores), delimiter=',')

def loadAndPlotScores(directory):
    scoreValues = np.loadtxt("%s/%s" % (directory, "gammaRobust.out"), delimiter=',')
    labels = [str(v) for v in GAMMAS]
    titles = ["Weighted Macro F1 Score on Motif Clusters", 
        "Weighted Macro F1 Score on All Clusters",
        "Accuracy Score on Motif Segments"]
    for i in range(3):
        plt.figure(i, figsize=(6, 7))
        ax = plt.subplot(111)
        plt.plot(labels, scoreValues[i])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),
          fancybox=True, shadow=True, ncol=3)
        plt.xlabel("Gamma")
        plt.ylim(ymin=0)
        plt.title(titles[i])
        if i == 2:
            plt.ylabel("Accuracy score on motif segments vs Noise")
        else:
            plt.ylabel("Weighted Macro F1 Score")
        # plt.legend()
    plt.show()


# saveScores("ordered_synthetic_retry2", "final_scores_retry2")
loadAndPlotScores("final_scores_retry2")


