from generate_synthetic_util import generate_data
import numpy as np
import sys

WINDOW_SIZE = 1
NUM_SENSORS = 5
SPARSITY = 0.2
RAND_SEED = 10


NUM_CLUSTERS = 10
GARBAGE_CLUSTERS = 10
CLUSTER_SEQUENCE = [6,7,8,9]
NUM_SEQS = 500 # number of macro segs
NUM_GARBAGE = 10 # number of garbage segs
LEN_SEGMENT = 25 # length of each segment

def createCorrect(outputFilename):
    assigns = []
    for _ in range(NUM_SEQS):
        assigns += np.random.choice(GARBAGE_CLUSTERS, NUM_GARBAGE).tolist()
        assigns += CLUSTER_SEQUENCE
    finalAssign = []
    for seg in assigns:
        finalAssign += [seg for _ in range(LEN_SEGMENT)]
    np.savetxt(outputFilename, np.array(finalAssign), delimiter=",", fmt='%d')

def createSegments(assigns, lenSegment):
    assignment = []
    correctAssignment = []
    for cluster, othercluster in assigns:
        correctAssignment += np.full(lenSegment, cluster, dtype=int).tolist()
        assignment.append(((cluster, othercluster), lenSegment))
    return assignment, correctAssignment


def createDataset2(outputFilename, outputCorrectName, noiseFraction, eps):
    num_clusters = 10
    garbage_clusters = 10 # 6
    clustersequence = [6, 7, 8, 9]
    numSeqs = 500
    numGarbage = 10  # number of garbage segments
    lenSegment = 25  # length of each segment
    # with probability epsilon, add noise

    # create assigns
    assigns = []
    for _ in range(numSeqs):
        garbage = np.random.choice(garbage_clusters, numGarbage).tolist()
        garbage = zip(garbage, [None for _ in range(len(garbage))])
        assigns += garbage
        noneVals = None
        if np.random.random() < eps:
            # pick a random one to edit
            noneVals = [None for _ in range(len(clustersequence))]
            chosenIdx = np.random.choice(len(clustersequence))
            chosenCluster = np.random.choice(garbage_clusters)
            noneVals[chosenIdx] = chosenCluster
        else:  # don't add noise
            noneVals = [None for _ in range(len(clustersequence))]
        assigns += list(zip(clustersequence, noneVals))

    #alts = [None for _ in range(len(assigns))]
    #assigns = list(zip(assigns, alts))
    assignment, correctAssignment = createSegments(assigns, lenSegment)
    generate_data(num_clusters, NUM_SENSORS, WINDOW_SIZE,
                  SPARSITY, assignment, outputFilename, RAND_SEED, noiseWeight=noiseFraction)

    correctAssignment = np.array(correctAssignment)
    np.savetxt(outputCorrectName, correctAssignment, delimiter=",", fmt='%d')


if __name__ == "__main__":
    createCorrect("blah.out")
    assert False
    assert len(sys.argv) >= 4
    mode = int(sys.argv[1])
    mode, output_file, correct_file = int(
        sys.argv[1]), sys.argv[2], sys.argv[3]
    if mode == 1:
        createDataset1(output_file, correct_file)
    else:
        assert len(sys.argv) == 6
        noisyFrac, epsVal = float(sys.argv[4]), float(sys.argv[5])
        createDataset2(output_file, correct_file, noisyFrac, epsVal)
