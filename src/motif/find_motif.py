from rstr_suffix.rstr_max import GetMotifs
from hmm import MotifHMM
import heapq
import numpy as np
from bitarray import bitarray
# https://code.google.com/archive/p/py-rstr-max/
# TODO: cite

def PerformAssignment(sequence, negLLMatrix, gamma):
    '''
    given the negLLMatrix for the null hypothesis
    '''
    motifs = FindMotifs(sequence) # find common motifs with scores
    instanceList = [] # list of (score, motif, indices)
    for m, logscore in motifs:
        motif_hmm = MotifHMM(negLLMatrix, m, 0.9)
        _, motifInstances = motif_hmm.SolveAndReturn()
        for motifIndices, likelihood in motifInstances:
            motifLength = motifIndices[-1] - motifIndices[0] + 1
            motifScore = (likelihood/motifLength) + logscore # should I be adding this?
            instanceList.append((motifScore, m, motifIndices))        
    heapq.heapify(instanceList)

def GreedyAssign(sequence, instanceList):
    result = sequence.copy()
    taken = bitarray(len(sequence))
    taken.setall(False)
    motifs = {} # motif to instances
    while len(instanceList) != 0:
        _, motif, indices = instanceList.heappop()
        assert len(motif) == len(indices) - 1
        gap = (indices[0], indices[-1])
        if taken[gap[0]:gap[1]+1].any():
            continue # this spot has already been taken
        taken[gap[0]:gap[1]+1] = 1 # reserve this spot
        for i in range(len(motif)): # fill this motif
            val = motif[i]
            start = indices[i]
            end = indices[i+1] + 1
            result[start:end] = val
        if tuple(motif) not in motifs:
            motifs[tuple(motif)] = []
        motifs[tuple(motif)].append(gap)
    return result, motifs

def FindMotifs(sequence):
    '''
    Return a list of [(motif, motif_score)]
    '''
    orig_indices, collapsed = Collapse(sequence)
    motif_results = GetMotifs(collapsed)
    processed_motif_list = []
    print orig_indices, collapsed, motif_results
    for length, incidences in motif_results:
        if length == 1:
            continue
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

def Collapse(sequence):
    start_end_indices = []
    collapsedString = []
    oldC = -1
    currSeg = []
    for i,c in enumerate(sequence):
        if c != oldC:
            if len(currSeg) == 1:
                # this is the end of this segment
                currSeg.append(i-1) # the last one was the end index
                start_end_indices.append(tuple(currSeg))
                currSeg = []
            currSeg.append(i)
            collapsedString.append(c)
            oldC = c
    assert len(currSeg) == 1
    currSeg.append(len(sequence) - 1)
    start_end_indices.append(tuple(currSeg))
    return start_end_indices, collapsedString

FindMotifs([1,1,2,3, 2,1,1,2,3, 1,1, 1,2,3,2])
