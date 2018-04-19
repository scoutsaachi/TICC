from TICC_solver import TICCSolver
import numpy as np
import sys

fname = "ticcsensor.csv"
mode = int(sys.argv[1])
solver = TICCSolver(window_size=1, number_of_clusters=10, lambda_parameter=5e-3, beta=100, threshold=2e-5, gamma=0.8, input_file=fname, num_proc=10, maxMotifs=30)
old_assign=None
usemotif = False
if mode == 1:
	print("using motif")
	old_assign = np.loadtxt("old_assignments.out", dtype=int)
	usemotif=True

(cluster_assignment, cluster_MRFs, motifs) = solver.PerformFullTICC(initialClusteredPoints=old_assign,useMotif=usemotif)
solver.CleanUp()
print(motifs)

if mode  == 1:
	fname = "assignments.out"
else:
	fname = "old_assignments.out"

print ("saving to ", fname)
np.savetxt(fname, cluster_assignment, fmt='%d')
