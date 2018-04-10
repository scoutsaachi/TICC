from rstr_suffix.rstr_max import GetMotifs
from hmm import MotifHMM
import heapq
import numpy as np
from bitarray import bitarray
# https://code.google.com/archive/p/py-rstr-max/
# TODO: cite


def PerformAssignment(sequence, negLLMatrix, gamma):
    '''
    Perform the motif guided sequence assignment

    Parameters
    -------
    sequence: the sequence of assignments outputted by TICC
    negLLMatrix: the negative log likelihood without any motifs
    gamma: aggressiveness of motif finding

    Returns
    -------
    result: the new assignment guided by the motifs
    motifs: the motifs found, as a dict of {motif: [(start1,end1), (start2,end2)...]}
    '''
    motifs = find_motifs(sequence)  # find common motifs with scores
    instanceList = []  # list of (score, motif, indices)
    # TODO: perform this in parallel
    for m, logscore in motifs:
        motif_hmm = MotifHMM(negLLMatrix, m, gamma)
        _, motifInstances = motif_hmm.SolveAndReturn() # (changepoints, score)
        for motifIndices, likelihood in motifInstances:
            motifLength = motifIndices[-1] - motifIndices[0] + 1
            motifScore = (likelihood/motifLength) + logscore
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
        for i in range(len(motif)):  # fill this motif
            val = motif[i]
            start = indices[i]
            end = indices[i+1] + 1
            result[start:end] = val
        key = tuple(motif)
        if key not in motifs:
            motifs[key] = []
        motifs[key].append(gap)
    return result, motifs


def find_motifs(sequence):
    '''
    Get the maximal motifs in the sequence along with their scores

    Returns
    -------
    A list of [(motif), score]
    '''
    orig_indices, collapsed = collapse(sequence)
    motif_results = GetMotifs(collapsed) # [(motif length), [<start_indices>]]
    processed_motif_list = []
    for length, incidences in motif_results:
        motif = collapsed[incidences[0]:incidences[0]+length]
        # TODO get original distances and incorporate them into motif description
        score = MotifScore(sequence, length, incidences)
        processed_motif_list.append(motif, score)
    return processed_motif_list


def MotifScore(sequence, length, incidences):
    '''
        perform motif score
        return a log score
        TODO implement
    '''
    return 0


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


find_motifs([1, 1, 2, 3, 2, 1, 1, 2, 3, 1, 1, 1, 2, 3, 2])
