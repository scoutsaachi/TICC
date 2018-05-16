from generate_synthetic_util import generate_data
import numpy as np
import sys

WINDOW_SIZE = 1
NUM_SENSORS = 5
SPARSITY = 0.2
RAND_SEED = 10


def createGarbage(numGarbageSegments, lenSegment, clusters):
    assignment = []
    correctAssignment = []
    clusterAssignments = np.random.choice(clusters, numGarbageSegments)
    for cluster in clusterAssignments:
        correctAssignment += np.full(lenSegment, cluster, dtype=int).tolist()
        assignment.append((cluster, lenSegment))
    return assignment, correctAssignment


def createSegments(assigns, lenSegment):
    assignment = []
    correctAssignment = []
    for cluster in assigns:
        correctAssignment += np.full(lenSegment, cluster, dtype=int).tolist()
        assignment.append((cluster, lenSegment))
    return assignment, correctAssignment


def createDataset1(outputFilename, outputCorrectName):
    num_clusters = 10
    garbage_clusters = 6
    clustersequence = [6, 7, 8, 9]
    numSeqs = 500
    numGarbage = 10  # number of garbage segments
    lenSegment = 25  # length of each segment

    # create assigns
    assigns = []
    for _ in range(numSeqs):
        assigns += np.random.choice(garbage_clusters, numGarbage).tolist()
        assigns += clustersequence

    assignment, correctAssignment = createSegments(assigns, lenSegment)
    generate_data(num_clusters, NUM_SENSORS, WINDOW_SIZE,
                  SPARSITY, assignment, outputFilename, RAND_SEED)
    correctAssignment = np.array(correctAssignment)
    np.savetxt(outputCorrectName, correctAssignment, delimiter=",", fmt='%d')

if __name__ == "__main__":
    assert len(sys.argv) == 4
    mode = int(sys.argv[1])
    mode, output_file, correct_file = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    if mode == 1:
        createDataset1(output_file, correct_file)

