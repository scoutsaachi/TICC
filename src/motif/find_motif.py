from rstr_suffix.rstr_max import GetMotifs
from hmm import MotifHMM
import heapq
import numpy as np
from bitarray import bitarray
from collections import Counter
# https://code.google.com/archive/p/py-rstr-max/
# TODO: cite


def PerformAssignment(sequence, negLLMatrix, gamma, MaxMotifs=None):
    '''
    Perform the motif guided sequence assignment

    Parameters
    -------
    sequence: the sequence of assignments outputted by TICC
    negLLMatrix: the negative log likelihood without any motifs
    gamma: aggressiveness of motif finding
    MaxMotifs: the max number of motifs allowed

    Returns
    -------
    result: the new assignment guided by the motifs
    motifs: the motifs found, as a dict of {motif: [(start1,end1), (start2,end2)...]}
    '''
    logFreqProbs = getFrequencyProbs(sequence)
    motifs = find_motifs(sequence, MaxMotifs)  # find common motifs with scores
    instanceList = []  # list of (score, motif, indices)
    # TODO: perform this in parallel
    for m, score in motifs:
        logscore = np.log(score)
        motif_hmm = MotifHMM(negLLMatrix, m, gamma)
        _, motifInstances = motif_hmm.SolveAndReturn()  # (changepoints, score)
        for motifIndices, neg_likelihood in motifInstances:
            logodds = computeLogOdds(
                neg_likelihood, m, motifIndices, logFreqProbs)
            motifScore = logodds + logscore
            instanceList.append((motifScore, m, motifIndices))
    heapq.heapify(instanceList)
    final_assignment, motif_result = greedy_assign(sequence, instanceList)
    return final_assignment, motif_result


def greedy_assign(sequence, instanceList):
    '''
    Given the motif assignments in instanceList, greedily assign
    the sequence to motifs

    Parameters
    ---------
    sequence: the list of current assignments without motif guiding
    instanceList: a heapified list of (motifScore, motif, [motifIndices])

    Returns
    -------
    value1: the new assignments
    value2: dictionary of motifs: {motif: [indices of instances]}
    '''
    result = sequence.copy()
    taken = bitarray(len(sequence))
    taken.setall(False)
    motifs = {}  # motif to instances
    while len(instanceList) != 0:
        _, motif, indices = instanceList.heappop()
        assert len(motif) == len(indices) - 1
        gap = (indices[0], indices[-1])
        if taken[gap[0]:gap[1]+1].any():
            # TODO: allow motifs to intrude on the ends of other motifs
            continue  # this spot has already been taken
        taken[gap[0]:gap[1]+1] = 1  # reserve this spot
        # fill motif
        result[gap[0]:gap[1]+1] = generateExpandedMotif(motif, indices)
        key = tuple(motif)
        if key not in motifs:
            motifs[key] = []
        motifs[key].append(gap)
    return result, motifs


def generateExpandedMotif(motif, motifIndices):
    '''
    given the motif and the indices corresponding to each
    stage of the motif, expand out the motif
    '''
    motifLength = motifIndices[-1] - motifIndices[0] + 1
    result = np.zeros(motifLength)
    for i in range(len(motif)):  # fill this motif
        val = motif[i]
        start = motifIndices[i]
        end = motifIndices[i+1] + 1
        result[start:end] = val
    return result


def find_motifs(sequence, maxMotifs=None):
    '''
    Get the maximal motifs in the sequence along with their scores

    Returns
    -------
    A list of [(motif), logscore]
    '''
    orig_indices, collapsed = collapse(sequence)
    logFreqProbs = getFrequencyProbs(collapsed)
    totLength = len(collapsed)
    motif_results = GetMotifs(collapsed)  # [(motif length), [<start_indices>]]
    processed_motif_list = []  # score, motif
    for length, incidences in motif_results:
        newCount = filterOverlapping(incidences, length)
        if newCount == 1: continue
        motif = collapsed[incidences[0]:incidences[0]+length]
        # TODO get original distances and incorporate them into motif description
        score = MotifScore(totLength, logFreqProbs, motif, newCount)
        processed_motif_list.append((score, motif))
    processed_motif_list.sort(reverse=True)  # sort by score
    if maxMotifs:
        processed_motif_list = processed_motif_list[:maxMotifs]
    scores = np.array([s for s, _ in processed_motif_list])
    scores = scores/np.linalg.norm(scores)
    result = [(processed_motif_list[i][1], scores[i]) for i in range(len(scores))]
    print result 
    return result


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


def MotifScore(totLength, logFreqProbs, motif, numIncidences):
    '''
        perform motif score: G-test
        return a score. do not log yet
    '''
    motifLength = len(motif)
    database_size = totLength - (motifLength - 1)

    logscore_indep = getMotifIndepProb(motif, logFreqProbs)
    log_E = logscore_indep + np.log(database_size)  # p_ind

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

def computeLogOdds(neg_likelihood, motif, motifIndices, logFreqProbs):
    likelihood = -1*neg_likelihood
    expanded_seq = generateExpandedMotif(motif, motifIndices)
    indiv_prob = getMotifIndepProb(expanded_seq, logFreqProbs)
    return likelihood - indiv_prob


find_motifs([1, 2, 3, 5, 1, 2, 3, 6, 1, 2, 3, 4, 5, 3, 4,
             4, 6, 5, 1, 2, 3, 4, 5, 4, 5, 3, 1, 2, 3, 1, 2])
