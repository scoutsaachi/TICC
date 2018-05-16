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


def dataset1(input_name, output_dir):
    beta = 70
    number_of_clusters = 10
    outputName = runNonMotifTICC(input_name, output_dir, number_of_clusters, beta, None)
    runHyperParameterTests(input_name, output_dir, number_of_clusters, beta, outputName)

def runHyperParameterTests(inputName, outputDir, clusters, beta, oldAssignmentsName):
    gammas = [0.7, 0.8, 0.9, 0.99]
    motifReqs = 5
    for g in gammas:
        gammaDir = "%s/%s/" % (outputDir, g)
        makeDir(gammaDir)
        runTest(1, inputName, gammaDir, clusters,
                beta, g, motifReqs, oldAssignmentsName)

def runNonMotifTICC(inputName, outputDir, clusters, beta, oldAssignmentsName):
    oldDir = "%s/old/" % outputDir
    makeDir(oldDir)
    return runTest(0, inputName, oldDir, clusters, beta, 1, 1, oldAssignmentsName)

def pickleObject(fname, data):
    f = open(fname, "wb")
    pickle.dump(data, f)
    f.close()

def makeDir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def runTest(mode, inputName, outputDir, clusters, beta, gamma, motifReq, oldAssignmentsName):
    print("TESTING %s" % (gamma))
    solver = TICCSolver(window_size=1, number_of_clusters=clusters, lambda_parameter=1e-3, beta=beta, threshold=2e-5,
                        gamma=gamma, input_file=inputName, num_proc=30, maxMotifs=50, motifReq=motifReq, maxIters=20)
    old_assign = None
    usemotif = False
    if mode == 1:
        print("using motif")
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
    assert len(sys.argv) == 3
    mode, input_fname, output_fdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    if mode == 1:
        dataset1(input_fname, output_fdir)

