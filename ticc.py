from TICC_solver import solve

BIC_ITER_NUMBER = 4 # The max number of iterations to run for BIC

def RunTicc(input_filename, output_filename, cluster_number=None, process_pool_size=1,
            window_size=1, lambda_param=11e-2, beta=None,
            maxIters=1000, threshold=2e-5, covariance_filename=None,
            input_format='matrix'):
    '''
    Required Parameters:
    -- input_filename: the path to the data file. see input_format below
    -- output_filename: the output file name to write the cluster assignments

    Optional Parameters:
    -- cluster_number: The number of clusters to classify. If not specified, then will
       perform on [3, 5, 10] and use BIC to choose
    -- process_pool_size: the number of processes to spin off for optimization. Default 1
    -- window_size: The size of the window for each cluster. Default 1
    -- lambda_param: sparsity penalty. Default 11e-2
    -- beta: the switching penalty. If not specified, will perform on 
       [50, 100, 200, 400] and then use BIC to choose
    -- maxIters: the maximum number of iterations to allow TICC to run. Default 1000
    -- threshold: the convergence threshold. Default 2e-5
    -- covariance_filename: if not None, write the covariance into this file
    -- file_type is the type of data file. the data file must 
       be a comma separated CSV. the options are:
       -- "matrix": a numpy matrix where each column is a feature and each
          row is a time step
       -- "graph": an adjacency list with each row having the form:
          <start label>, <end label>, value
    '''
    (cluster_assignment, cluster_MRFs) = TICC.solve(
        window_size=window_size,number_of_clusters=cluster_number, lambda_parameter=lambda_param,
        beta=beta, maxIters=maxIters, threshold=threshold,
        input_file=input_filename, num_proc=process_pool_size)

