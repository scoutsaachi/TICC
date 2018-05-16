import numpy as np
from snap import GenRndGnm, PNGraph

SAVE_COVS = False

'''
block_matrices = {}  # Stores all the block matrices
w = window_size
n = number_of_sensors
sparsity_inv_matrix = sparsity_inv_matrix
block_matrices = {}  # Stores all the block matrices
seg_ids = cluster_ids
'''


def genInvCov(size, low=0.3, upper=0.6, portion=0.2):
    ''' Generate off diagonal blocks'''
    portion = portion/2
    S = np.zeros((size, size))
    G = GenRndGnm(PNGraph, size, int((size*(size-1))*portion))
    for EI in G.Edges():
        value = (np.random.randint(2) - 0.5)*2 * \
            (low + (upper - low)*np.random.rand(1)[0])
        S[EI.GetSrcNId(), EI.GetDstNId()] = value
    S = S + S.T
    return np.matrix(S)


def genRandInv(size, low=0.3, upper=0.6, portion=0.2):
    ''' Generate diagonal blocks '''
    S = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if np.random.rand() < portion:
                value = (np.random.randint(2) - 0.5)*2 * \
                    (low + (upper - low)*np.random.rand(1)[0])
                S[i, j] = value
    return np.matrix(S)


def generate_inverse(rand_seed, n, sparsity, w):
    '''
    n: number of sensors
    w: window size
    sparsity: sparasity penalty (lambda)
    '''
    N = n*w
    block_matrices = {}
    np.random.seed(rand_seed)
    # Generate all the blocks
    # create diagonal blocks
    block_matrices[0] = genInvCov(size=n, portion=sparsity)
    for block in range(1, w):
        block_matrices[block] = genRandInv(size=n, portion=sparsity)

    inv_matrix = np.zeros([N, N])
    # go through all the blocks
    for i in range(w):
        for j in range(w):
            block_result = block_matrices[np.abs(i - j)]
            if i <= j:
                block_result = block_result.T
            inv_matrix[i*n:(i+1)*n, j*n:(j+1)*n] = block_result

    # Make the matrix positive definite
    lambda_min = min(np.linalg.eig(inv_matrix))

    inv_matrix = inv_matrix + (0.1 + abs(lambda_min)) * np.identity(N)
    eigs, _ = np.linalg.eig(inv_matrix)
    print "Modified Eigenvalues are:", np.sort(eigs)
    return inv_matrix


def generate_data(K, n, w, sparsity, assignments):
    '''
    K: number of clusters
    n: number of sensors
    w: window size
    
    '''
    

    # GENERATE POINTS
    N = n*w
    cluster_mean = np.zeros([n, 1])
    cluster_mean_stacked = np.zeros([N, 1])

    # Generate two inverse matrices
    cluster_invs = {}
    cluster_covs = {}
    for c in range(K):
        cluster_invs[c] = generate_inverse(c, n, sparsity, w)
        cluster_covs[c] = np.linalg.inv(cluster_invs[c])
        if SAVE_COVS:
            fn = "inv_cov_cluster_%s.csv" % (c)
            np.savetxt(fn, cluster_invs[c], delimiter=",", fmt='%1.6f')
            fn = "cov_cluster_%s.csv" % (c)
            np.savetxt(fn, cluster_covs[c], delimiter=",", fmt='%1.6f')

    # Data matrix
    T = len(assignments)
    Data = np.zeros((T, n))
    for counter in range(len(break_points)):
        break_pt = break_points[counter]
        cluster = seg_ids[counter]
        old_break_pt = 0 if counter == 0 else break_points[counter - 1]
        # construct the first point
        cov_matrix = cluster_covs[cluster][0:n, 0:n]
        new_mean = cluster_mean_stacked[n * (w-1):N]
        new_row = np.random.multivariate_normal(
            new_mean.reshape(n), cov_matrix)
        Data[0, :] = new_row
        for breakIndex in range(old_break_pt+1, break_pt):
            # generate the following points
            num = breakIndex if breakIndex < w else w-1
            cov_matrix = cluster_covs[cluster][0:(num+1)*n, 0:(num+1)*n]
            Sig22 = cov_matrix[(num)*n:(num+1)*n, (num)*n:(num+1)*n]
            Sig11 = cov_matrix[0:(num)*n, 0:(num)*n]
            Sig21 = cov_matrix[(num)*n:(num+1)*n, 0:(num)*n]
            Sig12 = np.transpose(Sig21)
            cov_mat_tom = Sig22 - \
                np.dot(np.dot(Sig21, np.linalg.inv(Sig11)),
                        Sig12)  # sigma2|1
            a = np.zeros((num*n, 1))
            for idx in range(num):
                a[idx*n:(idx+1)*n, 0] = Data[idx, :].reshape((n))
            new_mean = cluster_mean + np.dot(np.dot(Sig21, np.linalg.inv(Sig11)),
                                                (a - cluster_mean_stacked[0:(num)*n, :]))
            new_row = np.random.multivariate_normal(
                new_mean.reshape(n), cov_mat_tom)
            Data[breakIndex, :] = new_row

    print "length of generated Data is:", Data.shape[0]

    # save the generated matrix
    np.savetxt(out_file_name, Data, delimiter=",", fmt='%1.4f')
