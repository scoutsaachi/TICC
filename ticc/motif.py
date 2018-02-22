import numpy as np 
from test_utils import *

class HMM:
    def __init__(self, adjacencyMatrix, llMatrix, initialCosts=None):
        '''
        @params
        adjacencyMatrix - m x m matrix, where A[i,j] is the cost of going from state
            i to state j (infty if not possible)
        llMatrix: n x m matrix, where A[i,j] is the neg ll of data i being classified as j
        '''
        self.adjacencyMatrix = adjacencyMatrix
        self.llMatrix = llMatrix
        self.n = np.shape(llMatrix)[0] # of time steps
        self.m = np.shape(llMatrix)[1] # of clusters
        self.initialCosts = np.full(self.m, 0) if initialCosts is None else initialCosts # initial costs
        self.viterbiGrid = np.zeros((self.n,self.m)) # contains probabilities that will be filled by row
        self.backPointers = np.zeros((self.n,self.m)) # contains the back pointers for the viterbi path
        self.backPointers[0].fill(-1)

    def UpdateStep(self, ts):
        if ts == 0:
            self.viterbiGrid[0] = self.initialCosts + self.llMatrix[0]
            return
        prevCosts = np.reshape(self.viterbiGrid[ts - 1], (self.m,1)) # get the costs from the timestamp before
        withSwitchingCosts = prevCosts + self.adjacencyMatrix
        # we want the min value for each row 
        indices = np.argmin(withSwitchingCosts, axis=0) # backpointers
        maxvals =  withSwitchingCosts[indices, np.arange(self.m)]
        newRow = maxvals + self.llMatrix[ts] # the new ll
        self.viterbiGrid[ts] = newRow
        self.backPointers[ts] = indices

class MotifHMM(HMM):
    def __init__(self, llMatrix, motif, beta):
        '''
        This is the HMM for the motifs
        @params
        motif - a list of IDs that represent the motif, i.e [0,1,0,2]
        '''
        self.motif = motif
        numStates = len(motif)
        adjacencyMatrix = self.createAdjacencyMatrix(motif, numStates, beta)
        initCosts = np.full(numStates, np.infty)
        initCosts[0] = 0
        ll = np.take(llMatrix, motif, axis=1) # grab only the relevant likelihoods
        HMM.__init__(self, adjacencyMatrix, ll, initCosts)
    
    def createAdjacencyMatrix(self, motif, numStates, beta):
        r = np.full((numStates, numStates), np.infty)
        np.fill_diagonal(r, 0) # costs nothing to go to same cluster
        np.fill_diagonal(r[:,1:], beta ) # allow transitions to the next state
        return r
    
    def getEndingScore(self, ts):
        ''' return log likelihood of ending at a specific ts '''
        return -1 * self.viterbiGrid[ts, -1]

class NullHMM(HMM):
    def __init__(self, llMatrix, beta):
        numStates = np.shape(llMatrix)[1]
        adjacencyMatrix = self.createAdjacencyMatrix(numStates, beta)
        HMM.__init__(self, adjacencyMatrix, llMatrix, None)
    
    def createAdjacencyMatrix(self, numStates, beta):
        r = np.full((numStates, numStates), beta)
        np.fill_diagonal(r, 0) # no cost to go to same cluster
        return r
    
    def getMostLikelyScore(self, ts):
        return np.min(self.viterbiGrid[ts])


