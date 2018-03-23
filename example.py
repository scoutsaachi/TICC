import TICC_solver as TICC
import numpy as np
import sys

fname = "out.csv"
(cluster_assignment, cluster_MRFs) = TICC.solve(window_size = 5,number_of_clusters = 4, lambda_parameter = 11e-2, beta = 0.2, maxIters = 100, threshold = 2e-5, write_out_file = False, input_file = fname, prefix_string = "output_folder/", num_proc=1)

print(cluster_assignment)
np.savetxt('Results.txt', cluster_assignment, fmt='%d', delimiter=',')
