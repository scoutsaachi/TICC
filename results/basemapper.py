from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import csv 
import matplotlib.colors as mcolors
from collections import defaultdict
import sys

# python2.7 basemapper.py old_assignments/concat_5c_500b.out  raw_data/concat/latlongLarge.csv pics/concat_5c_500b_orig.out
motifs = {(0.99, 2): {(1.0, 6.0, 3.0, 5.0): [(1383, 1964), (996, 1382)], (5.0, 4.0, 3.0, 4.0, 3.0): [(2961, 3222), (2775, 2902)], (1.0, 5.0, 1.0, 3.0, 5.0): [(460, 705), (2139, 2350)]}, (0.7, 2): {(1.0, 6.0, 3.0, 5.0): [(1383, 1964), (996, 1382)], (4.0, 5.0, 4.0, 3.0): [(2738, 2831), (2955, 3162), (140, 387)], (1.0, 5.0, 1.0, 3.0, 5.0): [(460, 705), (2139, 2350)]}, (0.6, 2): {(1.0, 6.0, 3.0, 5.0): [(1383, 1964), (996, 1382)], (5.0, 7.0, 1.0): [(820, 962), (2566, 2631)], (5.0, 4.0, 3.0, 4.0, 3.0): [(2775, 2902), (157, 387), (2961, 3222)], (1.0, 5.0, 1.0, 3.0, 5.0): [(460, 705), (2139, 2350)]}}
# motifs = {(0.99, 2): {(6.0, 4.0): [(21363, 21410), (197, 243), (15003, 15042), (2485, 2628), (5207, 5322), (742, 834), (20869, 20947), (2417, 2484), (17675, 17740), (22384, 22448), (14402, 14452), (14273, 14322)], (3.0, 2.0): [(14764, 14817), (5718, 5749), (13367, 13552), (9220, 13340), (5352, 5517), (9055, 9219), (14818, 14923), (5620, 5717), (5058, 5151), (5750, 5842), (2331, 2416), (2900, 2971)], (0.0, 3.0, 6.0): [(17741, 18052), (928, 1023), (1024, 1113), (22055, 22181)], (4.0, 5.0, 3.0): [(1634, 2236), (3721, 4244), (15694, 16115), (16904, 17322), (16286, 16668), (21662, 22018), (21005, 21351), (3373, 3668), (4328, 4529), (21460, 21661), (15093, 15265), (22182, 22311), (2729, 2840)], (4.0, 5.0): [(3252, 3372), (16211, 16285)], (3.0, 7.0): [(14168, 14272), (6540, 6831), (1271, 1349), (13760, 14015), (6832, 6898), (677, 741), (5843, 6084), (5518, 5619), (18151, 20736), (13595, 13736), (20737, 20868), (17460, 17584), (7969, 9054), (6986, 7522), (2997, 3119), (14016, 14135), (18053, 18150), (6446, 6539), (6899, 6985)], (4.0, 0.0): [(16128, 16210), (601, 676), (5152, 5206), (275, 398), (3131, 3251)]}, (0.7, 2): {(6.0, 4.0): [(2417, 2484), (184, 243), (21363, 21410), (17675, 17740), (2485, 2519)], (3.0, 2.0): [(6084, 6445), (5060, 5151), (9132, 9219), (13488, 13552), (14763, 14817), (5718, 5749), (14818, 14923)], (4.0, 5.0, 3.0, 7.0, 6.0): [(15093, 15693), (15694, 16160), (16286, 16828)], (3.0, 7.0, 6.0): [(676, 830), (13761, 14175), (14176, 14318), (1370, 1513), (17772, 18071), (18072, 18177), (20757, 20914), (1271, 1369), (13605, 13760)], (3.0, 6.0, 4.0, 0.0): [(22355, 22620), (14953, 15092), (492, 675), (14373, 14762), (244, 398)], (4.0, 5.0, 3.0): [(22182, 22312), (4328, 4529)], (6.0, 4.0, 0.0): [(4648, 5015), (5207, 5350)], (3.0, 7.0, 1.0): [(22621, 35999), (6446, 7915)], (7.0, 2.0, 7.0): [(9220, 13340), (0, 183), (7916, 9131), (18195, 20737), (5400, 5646)], (0.0, 4.0, 5.0, 3.0): [(16829, 17315), (1567, 2236), (21002, 21351), (21455, 21650), (3709, 4244), (21651, 22018)], (0.0, 3.0, 6.0): [(914, 1022), (1023, 1113)], (4.0, 0.0, 4.0, 5.0, 4.0, 5.0, 3.0): [(2578, 2941), (3131, 3668)], (4.0, 0.0): [(14319, 14372), (831, 879), (5152, 5206)]}, (0.6, 2): {(6.0, 4.0): [(21363, 21410), (20869, 20946), (2417, 2485), (17675, 17740), (22384, 22448), (163, 243)], (3.0, 2.0): [(5060, 5151), (5644, 5717), (13488, 13552), (5718, 5749), (14818, 14923)], (0.0, 4.0, 5.0, 3.0): [(4326, 4530), (21651, 22018), (21001, 21351), (21454, 21650)], (6.0, 4.0, 0.0): [(272, 398), (14402, 14762), (742, 959), (4648, 5015)], (3.0, 7.0, 1.0): [(22620, 35999), (6446, 7915)], (3.0, 7.0, 6.0, 3.0, 7.0): [(13605, 14023), (1271, 1442), (14024, 14277), (17772, 18161)], (4.0, 5.0, 3.0, 0.0, 3.0): [(22179, 22383), (15694, 16115)], (7.0, 2.0, 7.0): [(9220, 13340), (15266, 15693), (18195, 20737), (7916, 9131)], (6.0, 4.0, 0.0, 3.0, 7.0): [(5207, 5444), (600, 741)], (4.0, 0.0): [(5152, 5206), (14319, 14372)], (4.0, 0.0, 4.0, 5.0, 4.0, 5.0, 3.0): [(16128, 16668), (3131, 3668), (1567, 2236), (2578, 2941), (15028, 15265), (16827, 17321), (3707, 4244)], (3.0, 2.0, 7.0, 2.0): [(5750, 6445), (17460, 17674)]}}
startPerc = 0
endPerc = 1

colormap = ['blue', 'red', 'green','orange','purple', 'yellow', 'gray', 'pink', 'springgreen', 'beige', 'brown', 'olivedrab', 'khaki', 'midnightblue', 'salmon', 'papayawhip']
# colormap = mcolors.CSS4_COLORS.keys()
motifValue = None
input_name, latlongname, output_name = sys.argv[1], sys.argv[2], sys.argv[3]
if len(sys.argv) == 6:
    motifDictIndex = (float(sys.argv[4]), int(sys.argv[5]))
    motifValue = motifs[motifDictIndex]
# get navigational data
def getNavData():
    lats = []
    longs = []
    with open(latlongname, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            lats.append(float(row[0]))
            longs.append(float(row[1]))
    return lats, longs

def getClusterData(lats, longs, motifs=None):
    assignMapLongs = defaultdict(list)
    assignMapLats = defaultdict(list)
    startVal = int(startPerc * len(lats))
    endVal = int(endPerc * len(lats))
    with open(input_name, 'r') as instream:
        print("Start, end: ", startVal, endVal)
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines]
        for i,a in enumerate(assigns):
            if i >= startVal and i <= endVal:
                assignMapLongs[a].append(longs[i])
                assignMapLats[a].append(lats[i])
    motifLats = {} # motif to lats
    motifLongs = {} # motif to longs
    endPointLats = {}
    endPointLongs = {}
    if motifs is not None:
        for k, incidents in motifs.items():
            # if tuple(k) != tuple([0,4,2,3,0,2,3]) and tuple(k)!=tuple([4,2]): continue
            motifLats[k] = []
            motifLongs[k] = []
            endPointLats[k] = []
            endPointLongs[k] = []
            for start,end in incidents:
                if start < startVal and end>endVal:
                    continue
                start = max(startVal, start)
                end = min(end, endVal)
                print(start,end)
                motifLats[k] += lats[start:end+1]
                motifLongs[k] += longs[start:end+1]
                endPointLats[k].append(lats[start])
                endPointLongs[k].append(longs[start])
                endPointLats[k].append(lats[end])
                endPointLongs[k].append(longs[end])
    return (assignMapLats, assignMapLongs), (motifLats, motifLongs, endPointLats, endPointLongs)

def plotValues(motifs=None):
    lats, longs = getNavData()
    assignData, motifData = getClusterData(lats, longs,motifs=motifs)
    assignMapLats, assignMapLongs = assignData
    startVal = int(startPerc * len(lats))
    endVal = int(endPerc * len(lats))
    long_min = min(longs)
    long_max = max(longs)
    long_center = (long_min + long_max)/2
    lat_min = min(lats)
    lat_max = max(lats)
    lat_center = (lat_min + lat_max)/2
    print long_min, long_max, lat_min, lat_max
    # plt.figure(figsize=(50,50))
    m = Basemap(llcrnrlon=long_min,
                llcrnrlat=lat_min,
                urcrnrlon=long_max,
                urcrnrlat=lat_max,
                lat_0 = lat_center,
                lon_0 = long_center,
                projection='merc',
                resolution='h')
    assignMapLats, assignMapLongs = assignData
    keys = assignMapLongs.keys()
    plt.figure(1)
    for k in keys:
        x, y = m(assignMapLongs[k], assignMapLats[k])
        m.scatter(x,y,10,marker='*',color=colormap[k], label="%s" % k)
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    art = []
    art.append(lgd)
    plt.savefig(output_name + "_assign.png", additional_artists=art, bbox_inches="tight")
    plt.figure(2)
    if motifs is not None:
        print("showing motifs")
        motifLats, motifLongs, endpointx, endpointy = motifData
        keys = motifLats.keys()
        for i,k in enumerate(keys):
            motifString = "%s" % list(k)
            x, y = m(motifLongs[k], motifLats[k])
            endX, endY = m(endpointy[k], endpointx[k])
            m.scatter(x,y,10,marker='*',color=colormap[i], label=motifString)
            m.scatter(endX, endY, 60, marker='|',color='black')
        art = []
        # lgd = plt.legend(loc=0, borderaxespad=0.)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        art.append(lgd)
        plt.savefig(output_name + "_motifs.png", additional_artists=art, bbox_inches="tight")
    plt.show()

plotValues(motifValue)