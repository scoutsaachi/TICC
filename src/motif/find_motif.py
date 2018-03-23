from rstr_suffix.rstr_max import Rstr_max
# https://code.google.com/archive/p/py-rstr-max/
# TODO: cite

def FindMotifs(sequence):
    pass

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
    print start_end_indices
    print collapsedString
