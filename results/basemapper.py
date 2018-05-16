from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import csv 
import matplotlib.colors as mcolors
from collections import defaultdict
import sys

# python2.7 basemapper.py old_assignments/concat_5c_500b.out  raw_data/concat/latlongLarge.csv pics/concat_5c_500b_orig.out
motifs = {(0.99, 2): ({(1, 6, 3, 5): [(1383, 1964), (996, 1382)], (5, 4, 3, 4, 3): [(2961, 3222), (2775, 2902)], (1, 5, 1, 3, 5): [(460, 705), (2139, 2350)]}, [((1, 5, 1, 3, 5), 41.034167780174215), ((5, 4, 3, 4, 3), 37.61599983078587), ((1, 6, 3, 5), 25.452823433601697)]), (0.6, 2): ({(1, 6, 3, 5): [(1383, 1964), (996, 1382)], (5, 4, 3, 4, 3): [(157, 387), (2775, 2902), (2961, 3222)], (1, 5, 1, 3, 5): [(460, 705), (2139, 2350)]}, [((5, 4, 3, 4, 3), 58.782753078646444), ((1, 5, 1, 3, 5), 40.91124561581006), ((1, 6, 3, 5), 25.329901269237528)]), (0.7, 2): ({(1, 6, 3, 5): [(1383, 1964), (996, 1382)], (5, 4, 3, 4, 3): [(2775, 2902), (2961, 3222)], (1, 5, 1, 3, 5): [(460, 705), (2139, 2350)]}, [((1, 5, 1, 3, 5), 41.034167780174215), ((5, 4, 3, 4, 3), 37.61599983078587), ((1, 6, 3, 5), 25.452823433601697)])}
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