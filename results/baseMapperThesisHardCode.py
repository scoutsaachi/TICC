from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import csv 
import matplotlib.colors as mcolors
from collections import defaultdict
import sys
import os

colormap = ['blue', 'brown', 'green','orange','purple', 'red', 'gray', 'springgreen', 'springgreen', 'beige', 'yellow', 'olivedrab', 'khaki', 'midnightblue', 'salmon', 'papayawhip']

SAVE = True
SHOW = False

def run():

    # percTuples = [(0,1,None),(0, 0.03, None), (0.07, 0.1, None)]
    # driver2Motifs = {(0.99, 2): {(1.0, 6.0, 3.0, 5.0): [(1383, 1964), (996, 1382)], (5.0, 4.0, 3.0, 4.0, 3.0): [(2961, 3222), (2775, 2902)], (1.0, 5.0, 1.0, 3.0, 5.0): [(460, 705), (2139, 2350)]}, (0.7, 2): {(1.0, 6.0, 3.0, 5.0): [(1383, 1964), (996, 1382)], (4.0, 5.0, 4.0, 3.0): [(2738, 2831), (2955, 3162), (140, 387)], (1.0, 5.0, 1.0, 3.0, 5.0): [(460, 705), (2139, 2350)]}, (0.6, 2): {(1.0, 6.0, 3.0, 5.0): [(1383, 1964), (996, 1382)], (5.0, 7.0, 1.0): [(820, 962), (2566, 2631)], (5.0, 4.0, 3.0, 4.0, 3.0): [(2775, 2902), (157, 387), (2961, 3222)], (1.0, 5.0, 1.0, 3.0, 5.0): [(460, 705), (2139, 2350)]}}
    # for startPerc, endPerc, motifSubset in percTuples:
    #     gamma = 0.6
    #     motifValues = driver2Motifs[(gamma, 2)]
    #     motifColorKeys = motifValues.keys()
    #     motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
    #     create3Plots(startPerc, endPerc, 2, gamma, motifValues, motifColorMap, motifSubset=motifSubset)

    percTuples = [(0,1,[(0.0,4.0,5.0,3.0), (3.0,7.0,6.0), (3.0,6.0,4.0,0.0), (0.0,3.0,6.0), (6.0,4.0)]), (0.477, 0.586, None), (0.01, 0.042, None), (0.375, 0.399, [(3.0, 7.0, 6.0)])]
    driver3Motifs = {(0.99, 2): {(6.0, 4.0): [(21363, 21410), (197, 243), (15003, 15042), (2485, 2628), (5207, 5322), (742, 834), (20869, 20947), (2417, 2484), (17675, 17740), (22384, 22448), (14402, 14452), (14273, 14322)], (3.0, 2.0): [(14764, 14817), (5718, 5749), (13367, 13552), (9220, 13340), (5352, 5517), (9055, 9219), (14818, 14923), (5620, 5717), (5058, 5151), (5750, 5842), (2331, 2416), (2900, 2971)], (0.0, 3.0, 6.0): [(17741, 18052), (928, 1023), (1024, 1113), (22055, 22181)], (4.0, 5.0, 3.0): [(1634, 2236), (3721, 4244), (15694, 16115), (16904, 17322), (16286, 16668), (21662, 22018), (21005, 21351), (3373, 3668), (4328, 4529), (21460, 21661), (15093, 15265), (22182, 22311), (2729, 2840)], (4.0, 5.0): [(3252, 3372), (16211, 16285)], (3.0, 7.0): [(14168, 14272), (6540, 6831), (1271, 1349), (13760, 14015), (6832, 6898), (677, 741), (5843, 6084), (5518, 5619), (18151, 20736), (13595, 13736), (20737, 20868), (17460, 17584), (7969, 9054), (6986, 7522), (2997, 3119), (14016, 14135), (18053, 18150), (6446, 6539), (6899, 6985)], (4.0, 0.0): [(16128, 16210), (601, 676), (5152, 5206), (275, 398), (3131, 3251)]}, (0.7, 2): {(6.0, 4.0): [(2417, 2484), (184, 243), (21363, 21410), (17675, 17740), (2485, 2519)], (3.0, 2.0): [(6084, 6445), (5060, 5151), (9132, 9219), (13488, 13552), (14763, 14817), (5718, 5749), (14818, 14923)], (4.0, 5.0, 3.0, 7.0, 6.0): [(15093, 15693), (15694, 16160), (16286, 16828)], (3.0, 7.0, 6.0): [(676, 830), (13761, 14175), (14176, 14318), (1370, 1513), (17772, 18071), (18072, 18177), (20757, 20914), (1271, 1369), (13605, 13760)], (3.0, 6.0, 4.0, 0.0): [(22355, 22620), (14953, 15092), (492, 675), (14373, 14762), (244, 398)], (4.0, 5.0, 3.0): [(22182, 22312), (4328, 4529)], (6.0, 4.0, 0.0): [(4648, 5015), (5207, 5350)], (3.0, 7.0, 1.0): [(22621, 35999), (6446, 7915)], (7.0, 2.0, 7.0): [(9220, 13340), (0, 183), (7916, 9131), (18195, 20737), (5400, 5646)], (0.0, 4.0, 5.0, 3.0): [(16829, 17315), (1567, 2236), (21002, 21351), (21455, 21650), (3709, 4244), (21651, 22018)], (0.0, 3.0, 6.0): [(914, 1022), (1023, 1113)], (4.0, 0.0, 4.0, 5.0, 4.0, 5.0, 3.0): [(2578, 2941), (3131, 3668)], (4.0, 0.0): [(14319, 14372), (831, 879), (5152, 5206)]}, (0.6, 2): {(6.0, 4.0): [(21363, 21410), (20869, 20946), (2417, 2485), (17675, 17740), (22384, 22448), (163, 243)], (3.0, 2.0): [(5060, 5151), (5644, 5717), (13488, 13552), (5718, 5749), (14818, 14923)], (0.0, 4.0, 5.0, 3.0): [(4326, 4530), (21651, 22018), (21001, 21351), (21454, 21650)], (6.0, 4.0, 0.0): [(272, 398), (14402, 14762), (742, 959), (4648, 5015)], (3.0, 7.0, 1.0): [(22620, 35999), (6446, 7915)], (3.0, 7.0, 6.0, 3.0, 7.0): [(13605, 14023), (1271, 1442), (14024, 14277), (17772, 18161)], (4.0, 5.0, 3.0, 0.0, 3.0): [(22179, 22383), (15694, 16115)], (7.0, 2.0, 7.0): [(9220, 13340), (15266, 15693), (18195, 20737), (7916, 9131)], (6.0, 4.0, 0.0, 3.0, 7.0): [(5207, 5444), (600, 741)], (4.0, 0.0): [(5152, 5206), (14319, 14372)], (4.0, 0.0, 4.0, 5.0, 4.0, 5.0, 3.0): [(16128, 16668), (3131, 3668), (1567, 2236), (2578, 2941), (15028, 15265), (16827, 17321), (3707, 4244)], (3.0, 2.0, 7.0, 2.0): [(5750, 6445), (17460, 17674)]}}
    for startPerc, endPerc, motifSubset in percTuples:
        motifValues = driver3Motifs[(0.7, 2)]
        motifColorKeys = motifValues.keys()
        motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
        create3Plots(startPerc, endPerc, 3, 0.7, motifValues, motifColorMap, motifSubset=motifSubset)


def create3Plots(startPerc, endPerc, driver_id, gamma, motifValues, motifColorMap, motifSubset=None):
    if motifSubset != None:
        tempValues = {}
        for m in motifSubset:
            tempValues[m] = motifValues[m]
        motifValues = tempValues
    output_folder = "final_theses_pics/%s/%s/%s_%s" % (driver_id, gamma, startPerc, endPerc)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    clusterColorMap = colormap[:8]
    lats, longs = getNavData("new_results/raw_data/%s/latlong.csv" % driver_id)
    startVal = int(startPerc * len(lats))
    endVal = int(endPerc * len(lats))

    long_min = min(longs[startVal:endVal])
    long_max = max(longs[startVal:endVal])
    long_center = (long_min + long_max)/2
    lat_min = min(lats[startVal:endVal])
    lat_max = max(lats[startVal:endVal]) + 0.0003
    print(long_min, long_max, lat_min, lat_max)
    lat_center = (lat_min + lat_max)/2
    m = Basemap(llcrnrlon=long_min,
            llcrnrlat=lat_min,
            urcrnrlon=long_max,
            urcrnrlat=lat_max,
            lat_0 = lat_center,
            lon_0 = long_center,
            projection='merc',
            resolution='h')
    
    #plot old assignments 
    old_assignment_name = "new_results/old_assignments/%s__clust8_beta100.0.out" % driver_id
    assignMapLats, assignMapLongs = getClusterData(old_assignment_name, lats, longs, startVal, endVal)
    plotClusterValues(m, assignMapLats, assignMapLongs, clusterColorMap)
    if SHOW: plt.show()
    if SAVE: plt.savefig("%s/old_assignment.png"% output_folder)
    plt.clf()

    # plot new assignments: 0.6
    new_assignment_name = "new_results/new_assignments/%s/%s__clust8_beta100.0_gamma%s_req2.out" % (driver_id, driver_id, gamma)
    assignMapLats, assignMapLongs = getClusterData(new_assignment_name, lats, longs, startVal, endVal)
    plotClusterValues(m, assignMapLats, assignMapLongs, clusterColorMap)
    if SHOW: plt.show()
    if SAVE: plt.savefig("%s/new_assignment.png"% output_folder)
    plt.clf()

    # plot motifs    
    motifLats, motifLongs, endPointLats, endPointLongs=getMotifData(lats, longs, startVal, endVal, motifValues)    
    lgd = plotMotifValues(m, motifLats, motifLongs, endPointLats, endPointLongs, motifColorMap)
    if SAVE:
        art = []
        art.append(lgd)
        plt.savefig("%s/new_motif.png"% output_folder, additional_artists=art, bbox_inches="tight")
    if SHOW: plt.show()
    plt.clf()

# get navigational data
def getNavData(latlongname):
    lats = []
    longs = []
    with open(latlongname, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            lats.append(float(row[0]))
            longs.append(float(row[1]))
    return lats, longs

def getClusterData(input_name, lats, longs, startVal, endVal):
    assignMapLongs = defaultdict(list)
    assignMapLats = defaultdict(list)
    with open(input_name, 'r') as instream:
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines]
        for i,a in enumerate(assigns):
            if i >= startVal and i <= endVal:
                assignMapLongs[a].append(longs[i])
                assignMapLats[a].append(lats[i])
    return assignMapLats, assignMapLongs

def plotClusterValues(m, assignMapLats, assignMapLongs, clusterColorMap):
    keys = assignMapLongs.keys()
    plt.figure(1)
    for k in keys:
        x, y = m(assignMapLongs[k], assignMapLats[k])
        m.scatter(x,y,20,marker='.',color=clusterColorMap[k], label="%s" % k)
    lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,mode="expand", borderaxespad=0., ncol=8)
    plt.axes().set_aspect('equal', 'datalim')
    return lgd
    # plt.show()
    # art = []
    # art.append(lgd)
    # plt.savefig(output_name + "_assign.png", additional_artists=art, bbox_inches="tight")

def getMotifData(lats, longs, startVal, endVal, motifs):
    motifLats = {} # motif to lats
    motifLongs = {} # motif to longs
    endPointLats = {}
    endPointLongs = {}
    if motifs is not None:
        for k, incidents in motifs.items():
            motifResultLats = []
            motifResultLongs = []
            endResultPointLats = []
            endResultPointLongs = []
            for start,end in incidents:
                for idx in range(start, end+1):
                    if idx < startVal or idx > endVal:
                        continue
                    motifResultLats.append(lats[idx])
                    motifResultLongs.append(longs[idx])
                if start > startVal and start < endVal:
                    endResultPointLats.append(lats[start])
                    endResultPointLongs.append(longs[start])
                if end > startVal and end < endVal:
                    endResultPointLats.append(lats[end])
                    endResultPointLongs.append(longs[end])
            if len(motifResultLats) != 0:
                motifLats[k] = motifResultLats
                motifLongs[k] = motifResultLongs
                endPointLats[k] = endResultPointLats
                endPointLongs[k] = endResultPointLongs
    return motifLats, motifLongs, endPointLats, endPointLongs

def plotMotifValues(m, motifLats, motifLongs, endPointLats, endPointLongs, motifColorMap):
    keys = motifLats.keys()
    for i,k in enumerate(keys):
        motifString = "%s" % list(k)
        x, y = m(motifLongs[k], motifLats[k])
        endX, endY = m(endPointLongs[k], endPointLats[k])
        m.scatter(x,y,20,marker='.',color=motifColorMap[k], label=motifString)
        # m.scatter(endX, endY, 60, marker='|',color='black')
    lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,mode="expand", borderaxespad=0., ncol=2)
    plt.axes().set_aspect('equal', 'datalim')
    return lgd

run()