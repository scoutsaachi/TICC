from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import csv 
import matplotlib.colors as mcolors
from collections import defaultdict
import sys
colormap = mcolors.CSS4_COLORS.keys()

input_name, output_name = sys.argv[1], sys.argv[2]

# input_name = "old_assignments.out"
# output_name = "ticc.png"

# get navigational data
def getNavData(startPerc, endPerc):
    lats = []
    longs = []
    with open('latlong.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            lats.append(float(row[0]))
            longs.append(float(row[1]))
    start = int(startPerc * len(lats))
    end = int(endPerc * len(lats))
    return lats[start:end], longs[start:end]

def getClusterData(startPerc, endPerc, lats, longs, motifs=None):
    assignMapLongs = defaultdict(list)
    assignMapLats = defaultdict(list)
    start = int(startPerc * len(lats))
    end = int(endPerc * len(lats))
    with open(input_name, 'r') as instream:
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines][start:end]
        for i,a in enumerate(assigns):
            assignMapLongs[a].append(longs[i])
            assignMapLats[a].append(lats[i])
    motifLats = {}
    motifLongs = {}
    if motifs is not None:
        for k, incidents in motifs.items():
            motifLats[k] = []
            motifLongs[k] = []
            for start,end in incidents:
                motifLats[k] += lats[start:end+1]
                motifLongs[k] += longs[start:end+1]
    return assignMapLats, assignMapLongs, motifLats, motifLongs

def plotValues(motifs=None):
    startPerc = 0
    endPerc = 1
    lats, longs = getNavData(startPerc, endPerc)
    long_min = min(longs)
    long_max = max(longs)
    long_center = (long_min + long_max)/2
    lat_min = min(lats)
    lat_max = max(lats)
    lat_center = (lat_min + lat_max)/2
    print long_min, long_max, lat_min, lat_max
    m = Basemap(llcrnrlon=long_min,
                llcrnrlat=lat_min,
                urcrnrlon=long_max,
                urcrnrlat=lat_max,
                lat_0 = lat_center,
                lon_0 = long_center,
                projection='merc',
                resolution='h')
    assignMapLats, assignMapLongs, motifLats, motifLongs = getClusterData(startPerc, endPerc, lats, longs,motifs=motifs)
    if motifs is None:
        keys = assignMapLongs.keys()
        for k in keys:
            x, y = m(assignMapLongs[k], assignMapLats[k])
            m.scatter(x,y,3,marker='o',color=colormap[k])
        plt.savefig(output_name)
        plt.show()
    else:
        keys = motifLats.keys()
        for i,k in enumerate(keys):
            x, y = m(motifLongs[k], motifLats[k])
            m.scatter(x,y,3,marker='o',color=colormap[i])
        plt.savefig(output_name)
        plt.show()

plotValues({(7.0, 9.0, 6.0, 9.0): [(982, 1296), (15450, 15680), (0, 575)], (3.0, 2.0, 8.0, 2.0, 4.0): [(16557, 16950), (22820, 23726)], (4.0, 1.0, 4.0): [(8199, 9588), (4210, 5382), (5398, 6357), (24938, 25731), (3064, 3619), (1939, 2746), (3761, 4054), (9937, 10551)], (8.0, 0.0, 8.0, 0.0): [(21291, 21828), (15347, 15449), (20635, 21290), (7178, 8198), (22284, 22475)], (9.0, 8.0, 0.0, 9.0): [(18519, 18832), (18061, 18369), (16017, 16137), (21858, 22137), (18869, 19332)], (0.0, 3.0, 2.0): [(24082, 24283), (17482, 17703), (17704, 17946), (33904, 34135), (20230, 20513)], (4.0, 1.0): [(33718, 33903), (15039, 15280), (14219, 14494), (17134, 17481), (13927, 14051), (13836, 13926), (27465, 28055), (19885, 20097)], (1.0, 8.0, 6.0): [(31702, 32862), (28199, 30222), (11039, 11390), (33091, 33542), (34232, 35999), (11440, 12853), (25732, 26957)], (9.0, 8.0, 7.0, 9.0): [(23827, 24081), (16464, 16556), (10606, 10846)]})