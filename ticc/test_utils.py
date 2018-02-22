import numpy as np


def GenerateFakeData(breakpoints, confidenceMean, confidenceVariance, numClusters):
    '''
    Create a probability matrix given the above breakpoints
    @ params
    breakpoints - a list of the length of each cluster, with the first and last values being random
    '''
    assert len(breakpoints) > 2
    r = [(-1, -1, -1)]*breakpoints[0]
    labels = ['_']*breakpoints[0]
    for i in range(1, len(breakpoints) - 1):
        identifier = i-1
        r += [(identifier, confidenceMean, confidenceVariance)] * breakpoints[i]
        label = chr(ord('A') + identifier)
        labels += [label]*breakpoints[i]
    r += [(-1, -1, -1)]*breakpoints[-1]
    labels += ['_']*breakpoints[1]
    likelihood = GenerateLikelihoods(r, numClusters)
    return likelihood, labels


def GenerateLikelihoods(assignments, numClusters):
    '''
    assignments: a list of (id, probability, variance) tuples. If -1 then garbage

    returns a negative log likelihood matrix
    '''
    results = []
    for c, p, v in assignments:
        if c == -1:
            # this is garbage, so just create a random probability
            a = np.random.rand(numClusters)
            a = a/np.sum(a)
            results.append(a)
        else:
            a = np.random.rand(numClusters)
            a[c] = 0
            val = min(1, np.random.normal(p, v))
            rest = 1-val
            a = rest*a/np.sum(a)
            a[c] = val
            results.append(a)
    results = -1*np.log(np.array(results))
    return results
