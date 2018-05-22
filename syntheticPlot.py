import numpy as np
from analyze_synthetic import getValidMappings 

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

saveScores("ordered_synthetic", "final_scores")
