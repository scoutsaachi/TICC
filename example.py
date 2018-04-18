from TICC_solver import TICCSolver
import numpy as np
import sys

fname = "example_data.txt"
solver = TICCSolver(window_size=1, number_of_clusters=8, lambda_parameter=11e-2, beta=600, maxIters=5, threshold=2e-5, gamma=0.9999, input_file=fname, num_proc=1)

old_assign = np.loadtxt("old_assignments.out", dtype=int)
(cluster_assignment, cluster_MRFs, motifs) = solver.PerformFullTICC(initialClusteredPoints=old_assign,useMotif=True)
solver.CleanUp()

np.savetxt("assignments.out", cluster_assignment, fmt='%d')
