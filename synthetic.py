from TICC_solver import TICCSolver
import numpy as np
import sys

def dataset2():
	beta = 70
	inputName = "synthetic/dataset2/input.csv"
	outputName = "synthetic/dataset2/old_assign.out"
	number_of_clusters = 10
	#runNonMotifTICC(inputName, outputName, number_of_clusters, beta, outputName)
	runHyperParameterTests(inputName, "synthetic/dataset2/new_assign", number_of_clusters,beta,outputName)

def dataset1():
	beta = 70
	inputName = "synthetic/dataset1/input.csv"
	outputName = "synthetic/dataset1/old_assign.out"
	number_of_clusters = 10
	#runNonMotifTICC(inputName, outputName, number_of_clusters, beta, outputName)
	runHyperParameterTests(inputName, "synthetic/dataset1/new_assign", number_of_clusters,beta,outputName)

def runHyperParameterTests(inputName, outputName, clusters, beta, oldAssignmentsName):
	gammas = [0.7, 0.8, 0.9, 0.99]
	motifReqs = [2]
	motifDict = {}
	for g in gammas:
		for m in motifReqs:
			motifs = runTest(1, inputName, outputName, clusters, beta, g, m, oldAssignmentsName)
			motifDict[(g, m)] = motifs
	print(motifDict)

def runNonMotifTICC(inputName, outputName, clusters, beta, oldAssignmentsName):
	runTest(0, inputName, outputName, clusters, beta, 1, 1, oldAssignmentsName)

def runTest(mode, inputName, outputName, clusters, beta, gamma, motifReq, oldAssignmentsName):
	print("TESTING %s %s" % (beta, gamma))
	solver = TICCSolver(window_size=1, number_of_clusters=clusters, lambda_parameter=1e-3, beta=beta, threshold=2e-5, gamma=gamma, input_file=inputName, num_proc=30, maxMotifs=50, motifReq=motifReq, maxIters=20)
	old_assign=None
	usemotif = False
	if mode == 1:
		print("using motif")
		old_assign = np.loadtxt(oldAssignmentsName, dtype=int)
		usemotif=True

	(cluster_assignment, cluster_MRFs, motifs, motifRanked) = solver.PerformFullTICC(initialClusteredPoints=old_assign,useMotif=usemotif)
	solver.CleanUp()
	print(motifs)
	if mode  == 1:
		fname = "%s_gamma%s_req%s.out" % (outputName, gamma, motifReq)
	else:
		fname = outputName
	print ("saving to ", fname)
	np.savetxt(fname, cluster_assignment, fmt='%d')
	return motifs, motifRanked

dataset1()
