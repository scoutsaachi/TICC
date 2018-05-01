from .rstr_suffix.rstr_max import GetMotifs
from .hmm import MotifHMM
import heapq
import numpy as np
from bitarray import bitarray
from collections import Counter, namedtuple,defaultdict
from scipy.stats import poisson,binom
# https://code.google.com/archive/p/py-rstr-max/

def PerformAssignment(sequence, negLLMatrix, solver):
    '''
    Perform the motif guided sequence assignment

    Parameters
    -------
    sequence: the sequence of assignments outputted by TICC
    negLLMatrix: the negative log likelihood without any motifs
    sovler: the TICC_solver calling this method

    Returns
    -------
    result: the new assignment guided by the motifs
    motifs: the motifs found, as a dict of {motif: [(start1,end1), (start2,end2)...]}
    '''
    sequence = [i.astype(int) for i in sequence]
    logFreqProbs = getFrequencyProbs(sequence)
    _, collapsed = collapse(sequence)
    totLength = len(collapsed)
    # find common motifs with scores
    motifs = find_motifs(sequence, solver.maxMotifs)
    nMotifsFound = len(motifs)
    instanceList = []  # list of (score, motif, indices)
    garbageCol, betaGarbage = getGarbageCol(
        sequence, negLLMatrix, solver.beta, solver.gamma)
    futures = [None]*nMotifsFound
    for i, motifTuple in enumerate(motifs):
        # motifTuple is motif lengths, motif
        futures[i] = solver.pool.apply_async(motifWorker,
            (totLength, motifTuple, solver.beta, solver.gamma, negLLMatrix,
            garbageCol, betaGarbage, logFreqProbs))
    instanceList = []
    for i in range(nMotifsFound):
        worker_result = futures[i].get()
        instanceList += worker_result
    instanceList.sort()    
    final_assignment, motif_result = greedy_assignv2(sequence, instanceList, solver.motifReq)
    print (motif_result)
    if np.all(final_assignment == sequence):
        print ("assignment not changed")
    return final_assignment, motif_result, computeFinalMotifScores(final_assignment, motif_result)

def computeFinalMotifScores(final_assignment, motif_result):
    logFreqProbs = getFrequencyProbs(final_assignment)
    motif_scores = []
    for m, incidences in motif_result.items():
        score = MotifScore(len(final_assignment), logFreqProbs, m, len(incidences))
        motif_scores.append((m, score))
    motif_scores.sort(reverse=True, key=lambda m:m[1])
    return motif_scores
    

def motifWorker(totLength, motifTuple, beta, gamma, negLLMatrix, garbageCol, betaGarbage, logFreqProbs):
    m, motifIncidenceLengths = motifTuple
    instanceList = []
    motif_hmm = MotifHMM(negLLMatrix, m, beta, gamma,
                         motifIncidenceLengths, garbageCol, betaGarbage)
    _, motifInstances = motif_hmm.SolveAndReturn()  # (changepoints)
    score = MotifScore(totLength, logFreqProbs, m, len(motifInstances))
    for motifIndices, neg_likelihood in motifInstances:
        logodds = computeLogOdds(m, motifIncidenceLengths, motifIndices, garbageCol, negLLMatrix)
        motifScore = logodds + score 
        instanceList.append((-1*motifScore, tuple(m), motifIndices))
    return instanceList


def getGarbageCol(sequence, negLLMatrix, beta, gamma):
    '''
    return the garbage column and the times where the garbage 
    column should be added by beta
    '''
    n = len(sequence)
    origVals = negLLMatrix[range(n), sequence] - np.log(gamma)
    betaGarbage = set()
    # add in where the original subsumed switching costs
    currValue = sequence[0]
    for i in range(1, n):
        newValue = sequence[i]
        if currValue != newValue:
            betaGarbage.add(i)
        currValue = newValue
    return origVals, betaGarbage


def greedy_assignv2(sequence, instanceList, motifReq):
    T = len(sequence)
    locked = bitarray(T) # locked corresponds to indices that have been irrecovocably taken
    locked.setall(False)
    tentative = [] # tentative corresponds to motifs that have made a bid on an index
    for i in range(T):
        tentative.append(set())
    result = sequence.copy() # the final result
    motifResult = defaultdict(set) # motifs to instance ids
    locked_vals = []

    def lock(motifIDX): # lock down an instance. assume that we're already in motif set, and already tentative
        nonlocal instanceList, locked, result, tentative, motifResult
        assert not isLocked(motifIDX)
        locked_vals.append(motifIDX)
        
        _, motif, indices = instanceList[motifIDX]
        oldSet = motifResult[motif].copy()
        gap = (indices[0], indices[-1])
        for i in range(gap[0], gap[1]+1):
            losers = tentative[i].copy()
            for loser in losers: # for each loser who lost out on this spot
                if loser == motifIDX: continue # not counting ourselves as losers
                _, loser_m, loserIndices = instanceList[loser]
                loserGap = (loserIndices[0], loserIndices[-1])
                for j in range(loserGap[0], loserGap[1]+1): # remove this loser from tentative
                    tentative[j].discard(loser)
                motifResult[loser_m].discard(loser) # remove this loser from the motifResult
            assert len(tentative[i]) == 1
            tentative[i] = None
        # lock it down
        locked[gap[0]:gap[1]+1] = 1
        result[gap[0]:gap[1]+1] = generateExpandedMotif(motif, indices)

    def makeTentative(motifIDX):
        nonlocal instanceList, tentative, motifResult
        _, motif, indices = instanceList[motifIDX]
        gap = (indices[0], indices[-1])
        # add ourselves to the result set
        motifResult[motif].add(motifIDX)
        for i in range(gap[0], gap[1]+1): # make ourselves tentative
            tentative[i].add(motifIDX)
    
    def isLocked(idx):
        nonlocal instanceList, locked
        _, _, indices = instanceList[idx]
        gap = (indices[0], indices[-1])
        return locked[gap[0]:gap[1]+1].any()

    for motifIndex, instance in enumerate(instanceList):
        score, motif, indices = instance
        assert len(motif) == len(indices) - 1
        if isLocked(motifIndex):
            # spot is not feasible
            continue
        # make tentative
        makeTentative(motifIndex) # are now tentative and have added to motifResult
        if len(motifResult[motif]) > motifReq:
            # was already a locked motif, so just update this index
            lock(motifIndex)
        if len(motifResult[motif]) == motifReq:
            # now a locked motif so lock everything in this set
            winners = motifResult[motif].copy()
            for newlyLockedIDX in winners:
                lock(newlyLockedIDX)
    finalResult = {}
    for m,instanceIdxSet in motifResult.items():
        if len(instanceIdxSet) < motifReq:
            continue
        gaps = []
        for idx in instanceIdxSet:
            assert isLocked(idx)
            _, _, indices = instanceList[idx]
            gap = (indices[0], indices[-1])
            gaps.append(gap)
        finalResult[m] = gaps
    return result, finalResult

def generateExpandedMotif(motif, motifIndices):
    '''
    given the motif and the indices corresponding to each
    stage of the motif, expand out the motif
    '''
    motifStart = motifIndices[0]
    motifLength = motifIndices[-1] - motifIndices[0] + 1
    result = np.zeros(motifLength, dtype='int')
    for i in range(len(motif)):  # fill this motif
        val = motif[i]
        start = motifIndices[i] - motifStart
        end = motifIndices[i+1] + 1 - motifStart
        result[start:end] = val
    return result

def checkPeriodic(S):
    S = "".join(map(str, S))
    doubleS = S + S
    curtailedDoubleS = doubleS[1:-1]
    return S in curtailedDoubleS

def replaceRedundancy(motifStr, candidateStrings):
    '''candidateList sorted from smallest to biggest'''
    ''' motifStr as string that is comma sep'''
    m = motifStr
    for i,cs in reversed(list(enumerate(candidateStrings))):
        replacementStr = str(-1*(i+1))+"."
        m = m.replace(cs, replacementStr)
    return m


def find_motifs(sequence, maxMotifs=None):
    '''
    Get the maximal motifs in the sequence along with their scores

    Returns
    -------
    A list of [(motif, motifIncidenceLengths, logscore)]
    '''
    orig_indices, collapsed = collapse(sequence)
    logFreqProbs = getFrequencyProbs(collapsed)
    print(logFreqProbs)
    totLength = len(collapsed)
    motif_results = GetMotifs(collapsed)  # [(motif length), [<start_indices>]]
    motif_results.sort(key=lambda mr: mr[0]) # sort so that the shortest is first
    print(totLength)
    candidateStrings = [] # list of candidate motifs in string form sep by commas
    candidates = [] # motifs, incidences
    alpha = 0.01/len(motif_results)

    for length, incidences in motif_results:
        numIncidences = filterOverlapping(incidences, length)
        if numIncidences == 1: continue
        # dynamic candidate replacement
        motif = collapsed[incidences[0]:incidences[0]+length]
        motifStr = ",".join([str(m) for m in motif])
        motifStrReplaced = replaceRedundancy(motifStr, candidateStrings)
        motifReconstructed = [int(float(m)) for m in motifStrReplaced.split(",")]
        log_prob_ind = getMotifIndepProb(motifReconstructed,  logFreqProbs)
        # calculate pscore
        pscore = 1-binom.cdf(numIncidences, totLength, np.exp(log_prob_ind))
        if pscore < alpha: # significant
            candidateStrings.append(motifStr)
            idx = -1*len(candidateStrings)
            logFreqProbs[idx] = np.log(float(numIncidences)/totLength) # update probs
            motifIncidenceLengths = inflateMotifLengths(incidences, orig_indices, len(motif))
            candidates.append((motif, motifIncidenceLengths))
    print(candidateStrings)
    return candidates


def filterOverlapping(incidences, length):
    count = 1
    currVal = incidences[0]
    for i in range(1, len(incidences)):
        newVal = incidences[i]
        if newVal - currVal < length:
            # this new value overlaps so don't count it
            continue
        else:
            currVal = newVal
            count += 1
    return count

def GetPoissonMotifLambda(totLength, logFreqProbs, motif):
    ''' return lambda for poisson estimating non-overlapping
        motif is an array
    '''
    logscore_indep = getMotifIndepProb(motif, logFreqProbs)
    currOverlapSum = 0
    currProb = logscore_indep
    motifStrArr = [str(m) for m in motif]
    motifStr = ",".join(motifStrArr)
    currMotifSuffix = motifStrArr[:]
    for i in range(len(motif)-1):
        deletedLetter = int(currMotifSuffix[-1])
        currMotifSuffix = currMotifSuffix[:-1]
        currProb += logFreqProbs[deletedLetter]
        # check if overlap is true
        currStr = ",".join(currMotifSuffix)
        motifStr = ",".join(motifStrArr)
        if motifStr.endswith(currStr):
            currOverlapSum += np.exp(currProb)
    lamb = totLength*(np.exp(logscore_indep)/(1+currOverlapSum))
    return lamb

def PoissonMotifScore(totLength, logFreqProbs, motif, numIncidences):
    lamb = GetPoissonMotifLambda(totLength, logFreqProbs, motif)
    prob = 1 - poisson.cdf(numIncidences, lamb)
    return prob


def MotifScore(totLength, logFreqProbs, motif, numIncidences):
    '''
        perform motif score: G-test
        return a score. 
    '''
    E = GetPoissonMotifLambda(totLength, logFreqProbs, motif)
    log_E = np.log(E)
    N = numIncidences
    score = 2*N*(np.log(N)-log_E)
    return score


def collapse(sequence):
    '''
    Given a sequence, collapse consecutive instances within the sequence

    Example: collapse([1,1,1,2,2,2,3,3,3,1,1,1,1]) ->
        ([(0, 2), (3, 5), (6, 8), (9, 12)], [1, 2, 3, 1])
    Returns
    -------
    start_end_indices: A list of the start and end indices for each segment
    collapsed_string: The list with all consecutive values collapsed
    '''
    start_end_indices = []
    collapsedString = []
    oldC = -1
    currSeg = []
    for i, c in enumerate(sequence):
        if c != oldC:
            if len(currSeg) == 1:
                # this is the end of this segment
                currSeg.append(i-1)  # the last one was the end index
                start_end_indices.append(tuple(currSeg))
                currSeg = []
            currSeg.append(i)
            collapsedString.append(c)
            oldC = c
    assert len(currSeg) == 1
    currSeg.append(len(sequence) - 1)
    start_end_indices.append(tuple(currSeg))
    return start_end_indices, collapsedString


def getFrequencyProbs(arr):
    '''
    Return the log probabilities of each state from it's appearance in the array
    '''
    freqs = Counter(arr)
    tot = len(arr)
    result = {}
    for cluster, num in freqs.items():
        result[cluster] = np.log(float(num)/tot)
    return result


def getMotifIndepProb(motif, logFreqProbs):
    # assumes that each state is independent, gets probability of motif ocurring
    logscore_indep = 0
    for val in motif:
        logscore_indep += logFreqProbs[val]
    return logscore_indep


def computeLogOdds(motif, motifIncidenceLengths, motifIndices, garbageCol, negLLMatrix):
    negLLSubset = negLLMatrix[motifIndices[0]:motifIndices[-1]+1, :]
    n = negLLSubset.shape[0]
    expanded_seq = generateExpandedMotif(motif, motifIndices)
    assert len(expanded_seq) == n
    likelihood = negLLSubset[range(n), expanded_seq.tolist()]
    likelihood = -1*np.sum(likelihood)

    # indiv_probs are just the garbage columns
    garbage_likelihoods = garbageCol[motifIndices[0]:motifIndices[-1]+1]
    indiv_prob = -1*np.sum(garbage_likelihoods)
    '''
    #distance probs:
    distanceProbs = 0
    means, stdevs = motifIncidenceLengths
    for indexPair in motifIndices:
        dist = indexPair[1]-indexPair[0]
        prob = scipy.stats.norm.pdf(dist, means, x)
        distanceProbs += np.log(prob)
    distanceProbs = distanceProbs/n
    '''
    return 2*(likelihood - indiv_prob)

def inflateMotifLengths(collapsedStartIndices, orig_indices, length):
    '''
    Parameters
    -----------
    collapsedIndices: the start indices in the collapsed version of each incidence
    origIndices: an array such that origIndices[i]=(origstart, origend) corresponds
        to the segment in the original sequence that got mapped to index i in the collapsed sequence
    length: length of the collapsed motif

    Returns
    ----------
    --Edit if collapsed motif is length K, return two lists: [means, stdevs], each
    of length K

    If the collapsed motif is length K, and there are N incidences of the motif,
    return a N x K matrix A where A[i,j] is the length of the jth motif segment 
    in incidence i
    '''
    N = len(collapsedStartIndices)
    K = length
    result = np.zeros((N, K))
    for i, collapsedStart in enumerate(collapsedStartIndices):
        for j in range(K):
            index = collapsedStart + j
            segment = orig_indices[index]
            segLength = segment[1] - segment[0] + 1
            result[i, j] = segLength
    means = np.mean(result, axis=0)
    stdevs = np.std(result, axis=0)
    assert len(means) == K
    return (means, stdevs)
