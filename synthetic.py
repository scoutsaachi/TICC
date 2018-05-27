from TICC_solver import TICCSolver
import numpy as np
import sys
import pickle
import os

# def dataset2():
# 	beta = 70
# 	inputName = "synthetic/dataset2/input.csv"
# 	outputName = "synthetic/dataset2/old_assign.out"
# 	number_of_clusters = 10
# 	#runNonMotifTICC(inputName, outputName, number_of_clusters, beta, outputName)
# 	runHyperParameterTests(inputName, "synthetic/dataset2/new_assign", number_of_clusters,beta,outputName)


def dataset(mode, input_name, output_dir):
    beta = 20 # used 20 earlier
    number_of_clusters = 10
    if mode == 1:
        outputName = "%s/old/assign.out" % output_dir
    else: outputName = runNonMotifTICC(input_name, output_dir, number_of_clusters, beta, None)
    runHyperParameterTests(input_name, output_dir, number_of_clusters, beta, outputName)

def runHyperParameterTests(inputName, outputDir, clusters, beta, oldAssignmentsName):
    gammas = [0.6, 0.8]
    motifReqs = 10
    for g in gammas:
        gammaDir = "%s/%s/" % (outputDir, g)
        makeDir(gammaDir)
        runTest(1, inputName, gammaDir, clusters,
                beta, g, motifReqs, oldAssignmentsName, 10)

def runNonMotifTICC(inputName, outputDir, clusters, beta, oldAssignmentsName):
    oldDir = "%s/old/" % outputDir
    makeDir(oldDir)
    return runTest(0, inputName, oldDir, clusters, beta, 1, 1, oldAssignmentsName, 15)

def pickleObject(fname, data):
    f = open(fname, "wb")
    pickle.dump(data, f)
    f.close()

def makeDir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def runTest(mode, inputName, outputDir, clusters, beta, gamma, motifReq, oldAssignmentsName, maxIters):
    print("TESTING %s" % (gamma))
    #maxIters used to be 30
    solver = TICCSolver(window_size=1, number_of_clusters=clusters, lambda_parameter=1e-3, beta=beta, threshold=2e-5,
                        gamma=gamma, input_file=inputName, num_proc=10, maxMotifs=50, motifReq=motifReq, maxIters=maxIters)
    old_assign = None
    usemotif = False
    if mode == 1:
        old_assign = np.loadtxt(oldAssignmentsName, dtype=int)
        usemotif = True
    (cluster_assignment, cluster_MRFs, motifs, motifRanked) = solver.PerformFullTICC(
        initialClusteredPoints=old_assign, useMotif=usemotif)
    solver.CleanUp()
    if usemotif:
        # save the motifs and motifsRanked
        motifFile = "%smotifs.pkl" % outputDir
        pickleObject(motifFile, motifs)
        motifRanked = "%smotifRanked.pkl" % outputDir
        pickleObject(motifRanked, motifs)

    outputName = "%sassign.out" % outputDir
    np.savetxt(outputName, cluster_assignment, fmt='%d')
    return outputName

if __name__ == "__main__":
    # mode of 1 to skip old assign
    assert len(sys.argv) == 4
    mode, input_fname, output_fdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    dataset(mode, input_fname, output_fdir)

