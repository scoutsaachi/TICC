from collections import Counter
import itertools
import heapq
import numpy as np
from sklearn.metrics import f1_score, confusion_matrix
import matplotlib.pyplot as plt
segLength = 200
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


def plot_confusion_matrix(cm, classes,output_name,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues,):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt), fontsize=16,
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.show()

def performMapping(arr, perm):
    return [perm[i] for i in arr]

def getAssigns(fname):
    with open(fname, 'r') as instream:
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines]
    return assigns

def getMappingDists(correctAssigns, testAssigns):
    T = len(correctAssigns)
    K = 10
    results = [np.zeros(K) for i in range(K)]
    for i in range(T):
        correct = correctAssigns[i]
        test = testAssigns[i]
        results[correct][test] += 1.0
    x = []
    indexes = [i for i in range(K)]
    for i in range(K):
        final_result = results[i]
        final_result *= -1
        # final_result /= -1*np.sum(results[i]) # invert order
        final_result = final_result.tolist()
        indices_result = list(zip(final_result, indexes[:]))
        indices_result.sort()
        x.append((indices_result,i))
        # resultStr = ["%s:%.2f" % (j, results[i][j]) for j in range(K)]
        # print("%s: %s" % (i, resultStr))
    heapq.heapify(x)
    taken = [False for _ in range(K)]
    mapping = [None for _ in range(K)]

    while len(x) != 0:
        r, idx = heapq.heappop(x)
        best_score, best_index = r[0]
        if taken[best_index]:
            r = r[1:]
            assert len(r) != 0
            heapq.heappush(x, (r, idx))
        else:
            mapping[idx] = best_index
            taken[best_index] = True
    return mapping

def getValidMappings(correctFname, assignFname):
    correctAssigns = getAssigns(correctFname)
    testAssigns = getAssigns(assignFname)
    mapping = getMappingDists(correctAssigns, testAssigns)
    mapped_correct = performMapping(correctAssigns, mapping)
    cared_about_values = mapping[6:]
    score = f1_score(mapped_correct, testAssigns, average='macro', labels=cared_about_values)
    score2 = f1_score(mapped_correct, testAssigns, average='macro')
    # cf = confusion_matrix(mapped_correct, testAssigns)
    # plot_confusion_matrix(cf, list(range(10)),"blah", title='Confusion matrix specific', normalize=True)
    print(score, score2)
    


getValidMappings("synthetic_testing/raw/correct.out", "synthetic_testing/old/assign.out")
getValidMappings("synthetic_testing/raw/correct.out", "synthetic_testing/0.8/assign.out")



# def run(ground_truth, input_name, input_map):
#     caredAboutValues = [7,8,9]
#     mapping = range(10)
#     for i,k in input_map.items():
#         mapping[i] = k
#     ground_truth = [g if g in caredAboutValues else -1 for g in ground_truth]
#     truth = []
#     for i in ground_truth: truth += [i]*segLength
#     with open(input_name, 'r') as instream:
#         lines = instream.readlines()
#         assigns = [int(val.strip()) for val in lines]
#     mappedAssigns = performMapping(assigns, mapping)
#     subsetOutAssigns = [g if g in caredAboutValues else -1 for g in mappedAssigns]
#     score = f1_score(truth, subsetOutAssigns, average='macro', labels=caredAboutValues)
#     return score

# def runMotif(motifTruth, m, size, output_name):
#     truth = np.full(size, 0)
#     for start, end in motifTruth:
#         truth[start:end+1] = 1
#     motifVals = np.full(size, 0)
#     for _, arr in m.items():
#         for start,end in arr:
#             motifVals[start:end+1] = 1
#     cf = confusion_matrix(truth, motifVals)
#     plot_confusion_matrix(cf, ["no motif", "motif"],output_name)
#     print(cf)

# def printMotifTable(input_map, motifs, ranks):
#     mapping = range(10)
#     for i,k in input_map.items():
#         mapping[i] = k
#     for m, rank in ranks:
#         num = len(motifs[m])
#         m_mapped = performMapping(m, mapping)
#         print "%s & %s & %s\\\\" % (m_mapped, num, rank)

    # plt.savefig(output_name)
    # plt.clf()

# def dataset1():
#     ground_truth = [1, 2, 4, 6, 5, 0, 5, 1, 6, 2, 7, 8, 7, 9, 0, 5, 1, 3, 6, 6, 4, 3, 6, 3, 7, 8, 7, 9,6, 5, 0, 6, 3, 2, 5, 5, 2, 3, 7, 8, 7, 9, 0, 6, 2, 2, 5, 1, 6, 3, 6, 4, 7, 8, 7, 9, 6, 1, 1, 4, 0, 0, 4, 4, 0, 3, 7, 8, 7, 9, 0, 3, 3, 2, 6, 1, 4, 6, 0, 2, 7, 8, 7, 9, 4, 4, 4, 4, 5, 0, 1, 6, 5, 5, 7, 8, 7, 9, 2, 1, 2, 1, 4, 1, 4, 4, 0, 0, 7, 8, 7, 9, 1,4, 0, 4, 3, 2, 5, 1, 0, 5, 7, 8, 7, 9, 0, 5, 4, 3, 1, 5, 2, 1, 6, 2, 7, 8, 7, 9]
#     dataset1Map = {9:3, 3:9}
#     score = run(ground_truth, "synthetic/dataset1/old_assign.out", dataset1Map)
#     print("truth f1: %s" % score)
#     gamma = [0.99, 0.9, 0.8]
#     for g in gamma:
#         fname = "synthetic/dataset1/new_assign_gamma%s_req3.out" % g
#         score = run(ground_truth, fname, dataset1Map)
#         print("gamma %s: %s" % (g, score))
#     motifTruth = [(2000, 2800), (4800, 5600), (7600, 8400), (10400, 11200), (13200, 14000), (16000, 16800), (18800, 19600), (21600, 22400), (24400, 25200), (27200, 28000)]
#     motifMap = {(0.8, 3): ({(5, 2, 6): [(9201, 9799), (4200, 4600), (18000, 18399), (9800, 10199), (23600, 24004), (1200, 1799)], (2, 6, 9, 7, 8, 7, 3): [(26599, 27999), (4601, 5596), (15799, 16800)], (8, 7): [(20800, 21200), (16801, 17601), (12400, 12800)], (7, 8, 7, 3, 0, 5): [(2000, 3205), (24401, 25600), (13198, 14205), (7600, 8602), (18798, 19600), (21599, 22408), (10402, 11199)]}, [((7, 8, 7, 3, 0, 5), 61.621988949751), ((2, 6, 9, 7, 8, 7, 3), 31.94092356664973), ((5, 2, 6), -26.828034119646517), ((8, 7), -33.31691138121139)]), (0.9, 3): ({(9, 7, 8, 7,3, 0): [(7200, 8602), (13000, 14203), (1800, 2999), (24383, 25398), (27000, 27999), (15800, 16800)], (7, 8, 7, 3): [(18795, 19600), (21599, 22399), (4601, 5596), (10402, 11199)]}, [((9, 7, 8, 7, 3, 0), 33.15969263523333), ((7, 8, 7, 3), -17.986767600037574)]), (0.99, 3): ({(7, 8, 7, 3): [(7399, 8399), (13001, 13999), (4601, 5596), (24383, 25200), (18795, 19600), (21599, 22399), (1999, 2798), (27200, 27999), (16002, 16800), (10402, 11199)]}, [((7, 8, 7, 3), -26.670261388931653)])}
#     for g in gamma:
#         runMotif(motifTruth, motifMap[(g,3)][0], len(ground_truth)*segLength, "synthetic/pics/1_cf%s.png" % g)
#         printMotifTable(dataset1Map,motifMap[(g,3)][0], motifMap[(g,3)][1])
# def dataset2():
#     ground_truth = [6, 3, 1, 4, 6, 7, 8, 3, 9, 1, 6, 2, 1, 0, 7, 8, 7, 9, 6, 2, 3, 1, 2, 7, 8, 7, 9, 1, 6, 0, 6, 3, 7, 8, 7, 9, 6, 3, 3, 0, 6, 7, 8, 7, 9, 4, 0, 4, 1, 3, 7, 8, 7, 9, 3, 0, 3, 4, 5, 7, 8, 7, 9, 1, 3, 4, 1, 6, 7, 8, 7, 9, 0, 6, 6, 0, 2, 7, 8, 7, 9, 5, 4, 2, 0, 6, 7, 8, 7, 9, 2, 2, 2, 4, 5, 7, 8, 0, 9, 4, 0, 0, 3, 0, 7, 8, 7, 9, 1, 2, 2, 6, 4, 7, 8, 7, 9, 4, 5, 1, 3, 2, 7, 8, 7, 9, 3, 3, 0, 1, 1, 7, 8, 7, 9, 3, 1, 1, 4, 3, 7, 8, 7, 9, 3, 1, 1, 0, 0, 7, 8, 7, 9, 0, 3, 4, 2, 2, 7, 8, 7, 9, 1, 0, 3, 5, 1, 7, 8, 7, 9, 6, 0, 2, 2, 0, 7, 8, 7, 9, 0, 5, 4, 6, 0, 7, 8, 6, 9, 5, 3, 3, 1, 4, 7, 8, 7, 9, 4, 5, 1, 0, 3, 7, 8, 7, 9, 0, 6, 2, 2, 1, 7, 8, 7, 9, 3, 6, 4, 5, 6, 7, 8, 7, 9, 4, 3, 5, 2, 1, 7, 8, 7, 9, 0, 5, 6, 5, 3, 7, 8, 7, 9, 4, 6, 1, 0, 1, 7, 8, 7, 9, 6, 5, 0, 6, 0, 7, 8, 7, 9, 5, 5, 6, 3, 2, 7, 8, 7, 9, 3, 2, 6, 5, 4, 7, 8, 5, 9, 6, 6, 3, 3, 0, 7, 8, 7, 9, 2, 2, 3, 2, 1, 7, 8, 7, 9, 5, 1, 2, 0, 4, 7, 8, 7, 9, 4, 0, 3, 4, 3, 7, 8, 7, 9, 5, 4, 4, 2, 1, 7, 8, 7, 9, 3, 4, 6, 0, 1, 7, 8, 7, 9, 3, 4, 3, 0, 2, 7, 8, 7, 9, 2, 2, 0, 3, 2, 7, 8, 7, 9, 0, 0, 3, 3, 5, 7, 8, 7, 9, 6, 0, 0, 5, 6, 7, 8, 3, 9, 1, 5, 6, 5, 3, 7, 8, 7, 9, 5, 5, 6, 5, 1, 7, 8, 7, 9, 3, 3, 3, 1, 5, 7, 8, 7, 9, 0, 1, 6, 4, 1, 7, 8, 7, 9, 0, 6, 0, 6, 4, 7, 8, 7, 9, 2, 6, 4, 5, 6, 7, 8, 7, 9, 6, 2, 0, 4, 4, 7, 8, 7, 9, 6, 5, 3, 4, 4, 7, 8, 7, 9, 2, 2, 2, 3, 3, 7, 8, 7, 9]
#     # 3,6,3,4 -> 7,8,7,9
#     # dataset2Map = {3:7, 6:8, 4:9, 7:3, 8:6, 9:4}
#     # 8,6,8,4 -> 7,8,7,9
#     dataset2Map = {8:7, 6:8, 4:9, 7:6, 9:4}
#     score = run(ground_truth, "synthetic/dataset2/old_assign.out", dataset2Map)
#     print("truth f1: %s" % score)
#     gamma = [0.99, 0.7, 0.6]
#     dataset2GammaMap = {
#         0.99: {8:7, 6:8, 4:9, 7:6, 9:4},
#         0.7: {3:7, 6:8, 4:9, 7:3, 8:6, 9:4},
#         0.6: {3:7, 6:8, 4:9, 7:3, 8:6, 9:4}
#     }
#     for g in gamma:
#         fname = "synthetic/dataset2/new_assign_gamma%s_req3.out" % g
#         score = run(ground_truth, fname, dataset2GammaMap[g])
#         print("gamma %s: %s" % (g, score))
#     motifTruth = [(1000, 1800), (2800, 3600), (4600, 5400), (6400, 7200), (8200, 9000), (10000, 10800), (11800, 12600), (13600, 14400), (15400, 16200), (17200, 18000), (19000, 19800), (20800, 21600), (22600, 23400), (24400, 25200), (26200, 27000), (28000, 28800), (29800, 30600), (31600, 32400), (33400, 34200), (35200, 36000), (37000, 37800), (38800, 39600), (40600, 41400), (42400, 43200), (44200, 45000), (46000, 46800), (47800, 48600), (49600, 50400), (51400, 52200), (53200, 54000), (55000, 55800), (56800, 57600), (58600, 59400), (60400, 61200), (62200, 63000), (64000, 64800), (65800, 66600), (67600, 68400), (69400, 70200), (71200, 72000), (73000, 73800), (74800, 75600), (76600, 77400), (78400, 79200), (80200, 81000), (82000, 82800), (83800, 84600), (85600, 86400), (87400, 88200), (89200, 90000)]
#     motifMap = {(0.99, 3): ({(8, 6): [(73000, 73400), (19000, 19399), (999, 1398), (36999, 37399), (55000, 55399)], (6, 3, 4): [(28200, 28798), (89402, 89999), (78600, 79200)], (8, 6, 8, 4): [(53195, 53999), (4597, 5400), (51400, 52203), (24398, 25200), (42400, 43202), (17201, 18002), (2799, 3600), (15400, 16201), (6400, 7200), (31600, 32400), (33399, 34199), (40599, 41399), (46000, 46800), (87399, 88199), (11800, 12599), (35201, 36000), (49600, 50399), (82000, 82799), (83800, 84599), (8200, 8999), (22600, 23399), (26200, 26999), (38800, 39599), (65800, 66599), (13601, 14399), (56801, 57599), (20800, 21598), (47800, 48598), (63997, 64795), (69404, 70202), (60402, 61199), (74802, 75599), (71203, 71999), (29806, 30600), (9999, 10800), (58603, 59401), (67598, 68399), (44200, 45000), (76597, 77398), (62198, 62999), (85600, 86399), (80204, 80999)]}, [((8, 6, 8, 4), -3.103308115613693), ((6, 3, 4), -9.601666714133371), ((8, 6), -59.80152641398994)]), (0.6, 3): ({(9, 3, 6, 3, 4, 1, 9, 0): [(55000, 56198), (45799, 47399), (999, 2199), (4596, 5799), (74400, 76198), (11600, 12802), (58398, 59802), (79993, 81399), (51399, 52801), (13597, 14998), (72999, 74399), (42201, 43599), (31599, 32599), (63800, 64806)], (4, 9, 0): [(83400, 83799), (76200, 76597), (79400, 79799)], (1, 7, 9, 3, 6, 3, 4, 0, 1): [(70204, 72599), (77396, 79399), (32600, 34599), (40200, 41599), (25600, 27002), (49200, 50603), (29400, 30799), (44195, 45000), (24393, 25200), (8197, 8999), (36799, 37601), (65400, 66599), (18798, 19600), (9599, 10804), (2600, 3800), (34991, 36200), (6199, 7401), (88800, 89999), (56598, 57608), (83800, 84804), (27799, 28801), (85599, 86600), (47400, 48598), (53194, 54001), (69402, 70203), (82000, 82804)], (3, 6, 3, 4, 9, 5, 2): [(15400, 16821), (62198, 63799), (67596, 68799), (17199, 18599), (20800, 22201)], (5, 3, 6, 3, 4, 5, 9): [(38600, 40199), (60201, 61402), (87000, 88199), (22399, 24199)], (4, 0, 7, 1): [(0, 599), (66998, 67399), (19996, 20400), (20401, 20799)]}, [((1, 7, 9, 3, 6, 3, 4, 0, 1), 565.996155842344), ((9, 3, 6, 3, 4, 1, 9, 0), 207.3598056635494), ((3, 6, 3, 4, 9, 5, 2), 48.58586585042973), ((5, 3, 6, 3, 4, 5, 9), 37.29175225248196), ((4, 0, 7, 1), -4.573148433743643), ((4, 9, 0), -22.31472318986021)]), (0.7, 3): ({(3, 6, 3, 4, 9, 5, 2): [(62198, 63799), (15400, 16799), (87399, 88800), (58603, 60004), (17199, 18599), (28000, 29406), (20800, 22198), (47800, 48802), (4597, 5599), (63997, 65199), (51399, 52601), (67596, 68799), (56801, 57999), (31600, 32602), (74802, 76001)], (1, 7, 9, 3, 6): [(20000, 20602), (36799, 37399), (27000, 27601)], (9, 0, 3, 6): [(73806, 74401), (1801, 2199), (47000, 47401)], (5, 9, 0, 3, 6, 8, 4): [(43595, 45000), (12999, 14399), (83202, 84599), (79798, 80999), (32988, 34199), (76192, 77398), (77998, 79200), (25798, 26999), (11400, 12599), (85200, 86399), (81795, 82799), (38600, 39599), (9400, 10800), (8002, 8999), (22399, 23399), (65599, 66599), (70999, 71999), (45799, 46800), (60201, 61199), (42200, 43202), (49401, 50401), (40587, 41399), (72598, 73405), (6394, 7200), (89199, 89999), (53195, 53999), (69404, 70202), (2798, 3600), (599, 1407), (29804, 30600), (35200, 36000), (24398, 25200), (18600, 19403)], (4, 0, 7, 1): [(14598, 15199), (66998, 67399), (55598, 56800)]}, [((5, 9, 0, 3, 6, 8, 4), 533.959932691621), ((3, 6, 3, 4, 9, 5, 2), 201.5734701818262), ((1, 7, 9, 3, 6), 3.165125283155381), ((4, 0, 7, 1), -5.495244866747187), ((9, 0, 3, 6), -10.689224355288015)]), (0.8, 3): ({(1, 0, 8, 6, 3, 4): [(16800, 18002), (7800, 8999), (29400, 30600), (56598, 57599), (2600, 3600), (34991, 36000), (51200, 52203), (13404, 14399), (83600, 84599), (43998, 45000), (20601, 21598), (47795, 48598), (62197, 62999), (38798, 39599), (82000, 82799), (89200, 89999), (74800, 75599), (6400, 7200), (27999, 28798), (87398, 88199), (22600, 23399), (60398, 61199), (85597, 86399), (9998, 10800)], (2, 9, 3, 6, 8, 4): [(69200, 70202), (41800, 43202), (15200, 16201), (49400, 50401), (70999, 71999), (58200, 59401), (25799, 26999), (33001, 34199), (63599, 64795), (45600, 46800), (77992, 79200), (31200, 32400), (76191, 77398), (67400, 68399), (65599, 66599), (4399, 5400), (11600, 12599), (24200, 25200), (40585, 41399), (79993, 80999), (52997, 53999)], (4, 9, 0): [(1600, 2199), (55598, 56198), (79400, 79798), (72600, 72999), (46999, 47399), (37599, 38001), (73600, 74399)], (3, 6): [(73000, 73400), (999, 1399), (36999, 37399), (19000, 19399), (55000, 55399)]}, [((1, 0, 8, 6, 3, 4), 264.1494832911775), ((2, 9, 3, 6, 8, 4), 213.9264484851046), ((4, 9, 0), -39.81118144279296), ((3, 6), -54.031775170440916)])}
#     for g in gamma:
#         runMotif(motifTruth, motifMap[(g,3)][0], len(ground_truth)*segLength, "synthetic/pics/2_cf%s.png" % g)
#         printMotifTable(dataset2GammaMap[g],motifMap[(g,3)][0], motifMap[(g,3)][1])
# dataset1()
# print("------")
# dataset2()