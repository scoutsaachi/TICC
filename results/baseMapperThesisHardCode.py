from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import csv 
import matplotlib.colors as mcolors
from collections import defaultdict
import sys
import os

colormap = ['blue', 'brown', 'green','orange','purple', 'red', 'gray', 'springgreen', 'fuchsia', 'hotpink', 'yellow', 'olivedrab', 'navy', 'midnightblue', 'salmon', 'turquoise', 'black', 'mediumaquamarine','darkslategray', 'turquoise', 'lightcyan','crimson','beige']

SAVE = False
SHOW = True
resultsDir = "conf_paper"
outputDir = "placeholder"

def run():

    driver6Motifs = {(0.99, 3): ({(1, 3, 0, 3): {(6529, 6992), (19926, 20194), (14013, 14144), (18503, 19050)}, (5, 7, 6): {(16006, 16371), (20406, 20626), (4872, 5245), (1338, 2058), (29512, 29734), (19336, 19546), (22635, 24400), (31804, 31973), (25904, 26000), (6993, 7317), (13940, 14012), (9691, 9824), (33528, 33662), (22101, 22260), (4509, 4871), (14145, 14365), (2781, 3386), (9025, 9439), (31597, 31779), (9914, 10240), (26974, 27330), (24681, 25065), (13513, 13629), (11711, 11797)}, (6, 4, 2): {(382, 539), (27781, 28007), (12482, 12576), (17011, 17100), (7867, 7966), (18397, 18478), (34918, 35017), (17892, 17948), (29148, 29373), (33934, 34119), (868, 1274), (12902, 12985), (29778, 29948), (34243, 34343), (17650, 17716), (606, 685), (6270, 6349), (27511, 27626), (5477, 5547), (16639, 16745), (8050, 8120), (33753, 33828), (17397, 17608), (10596, 10701), (26140, 26260), (8768, 8901), (21768, 21852), (21903, 21999), (34541, 34657)}, (3, 5, 7): {(15810, 15934), (15353, 15451), (15531, 15641), (21276, 21767), (25473, 25903)}, (2, 0, 1): {(2074, 2560), (14795, 14964), (31348, 31478), (7525, 7664)}}, [((6, 4, 2), 209.4910155639129), ((5, 7, 6), 133.70355141175907), ((1, 3, 0, 3), 34.09275748463016), ((3, 5, 7), 19.347145296481184), ((2, 0, 1), 13.333705460657582)]), (0.9, 3): ({(5, 7, 6): {(4113, 4871), (4872, 5257), (20406, 20626), (1338, 2058), (24704, 25065), (29512, 29734), (19336, 19546), (31804, 31973), (25904, 26000), (22100, 22260), (6993, 7317), (13940, 14012), (9691, 9824), (33528, 33662), (26973, 27330), (14145, 14365), (2781, 3386), (31597, 31779), (9914, 10240), (13513, 13629), (11800, 12035), (30108, 30650)}, (4, 2, 0, 1): {(32226, 32431), (28261, 28669), (14778, 14964), (12510, 12901), (412, 605)}, (2, 0, 1): {(2074, 2560), (31348, 31478), (7525, 7664)}, (1, 3, 0, 3): {(11052, 11362), (19587, 20194), (11416, 11647), (14013, 14144), (6529, 6992)}, (6, 4, 2, 0, 1): {(17650, 17891), (18397, 18804), (606, 822), (34541, 34737), (10596, 10874), (7867, 8049), (33753, 33902), (5477, 5627), (6270, 6515), (17397, 17649), (26140, 26343), (33934, 34242), (17011, 17396)}, (6, 4, 2): {(34243, 34344), (8050, 8120), (34918, 35017), (27781, 28007), (21768, 21852)}, (1, 6, 4, 2): {(27331, 27626), (21853, 21999), (16393, 16745)}, (3, 5, 7): {(25472, 25903), (15642, 15809), (26577, 26890), (11675, 11768), (15810, 15934), (18887, 19250), (21275, 21767), (13087, 13243), (15353, 15451), (22634, 24277), (15531, 15641), (15935, 16287)}, (6, 4, 2, 0): {(12902, 13020), (17892, 18030), (9332, 9628), (868, 1337), (29778, 30056), (8767, 9024), (29148, 29478)}}, [((6, 4, 2, 0, 1), 178.97359795275588), ((5, 7, 6), 121.78753566121462), ((3, 5, 7), 66.2836678002839), ((6, 4, 2, 0), 59.55555261337467), ((1, 3, 0, 3), 42.39146334695937), ((4, 2, 0, 1), 41.97080669410424), ((1, 6, 4, 2), 20.039872853559594), ((6, 4, 2), 18.402117636514475), ((2, 0, 1), 7.996305752858756)]), (0.6, 3): ({(5, 7, 6): {(33528, 33662), (20406, 20626), (24672, 25074), (20727, 21043)}, (3, 0, 3, 5, 7): {(25402, 25903), (21169, 21767), (18002, 18126), (8962, 9332), (1317, 1944), (19251, 19502), (13087, 13243), (9603, 9791), (26491, 26890), (15295, 15451), (22510, 24277), (11580, 11674), (2759, 3325), (20152, 20405), (18805, 19250)}, (3, 0, 3): {(11249, 11362), (4040, 4102), (24605, 24667), (26891, 26973)}, (7, 6, 1, 6, 4, 2, 0, 1): {(25931, 26343), (16807, 17397), (31872, 32431), (27176, 27628), (12036, 12901), (28072, 28669), (18180, 18804), (7665, 8049)}, (3, 5, 7): {(13314, 13447), (13921, 13987), (8656, 8767), (15642, 15809), (15810, 15934), (29041, 29147), (15531, 15641)}, (6, 4, 2, 0): {(29148, 29444), (31238, 31430), (7463, 7611), (12902, 13020), (25123, 25364), (29778, 30000), (10596, 10732), (35360, 35999), (6270, 6387), (381, 539)}, (0, 3, 0, 5, 7): {(22000, 22223), (3624, 3933), (31541, 31760)}, (4, 2, 0): {(14778, 14869), (8828, 8961), (17432, 17612)}, (2, 0, 1): {(14965, 15090), (19761, 20068), (32488, 32748)}, (3, 5, 7, 6, 7, 1): {(29478, 29764), (6946, 7421), (9884, 10513), (14076, 14750), (30156, 30779), (13463, 13803), (4883, 5476), (32942, 33443), (15935, 16638), (11675, 11977), (4112, 4882)}, (6, 4, 2, 0, 1): {(17650, 17891), (606, 822), (34541, 34737), (1945, 2560), (33753, 33902), (5477, 5627), (33934, 34242), (28877, 29002)}, (6, 4, 2): {(34918, 35017), (17892, 17948), (27780, 28007), (8050, 8131), (16639, 16745), (868, 1274), (21768, 21853), (34779, 34917)}}, [((3, 0, 3, 5, 7), 198.39510609129724), ((7, 6, 1, 6, 4, 2, 0, 1), 194.96541125201068), ((3, 5, 7, 6, 7, 1), 179.94129461898518), ((6, 4, 2, 0, 1), 101.66799400946505), ((6, 4, 2, 0), 89.66273098881598), ((6, 4, 2), 37.687886495089366), ((0, 3, 0, 5, 7), 29.14667026791637), ((3, 5, 7), 26.545064756087598), ((5, 7, 6), 9.044865694933476), ((3, 0, 3), 8.860546414766082), ((4, 2, 0), 8.607370767074576), ((2, 0, 1), 6.851353083297655)]), (0.7, 3): ({(5, 7, 6): {(31804, 31973), (33528, 33662), (25904, 26000), (31597, 31779), (22100, 22260), (11800, 12035), (20406, 20626), (3690, 3969)}, (3, 0, 3, 5, 7): {(25402, 25903), (15296, 15451), (29001, 29147), (21169, 21767), (29444, 29615), (18002, 18126), (26891, 27196), (19251, 19502), (1317, 1944), (9825, 10159), (13087, 13243), (26491, 26890), (9603, 9791), (20152, 20405), (22510, 24277), (11580, 11674), (24647, 25014), (8996, 9332), (18805, 19250)}, (2, 0, 1): {(14795, 14964), (7525, 7693), (32488, 32748)}, (3, 5, 7, 6, 7, 1): {(6946, 7421), (14076, 14750), (30156, 30779), (13463, 13803), (4883, 5476), (2759, 3411), (32942, 33443), (15935, 16638), (4112, 4882)}, (6, 4, 2, 0, 1): {(17650, 17891), (381, 605), (5477, 5627), (18397, 18804), (32182, 32431), (12482, 12901), (606, 822), (34541, 34737), (1945, 2560), (7867, 8049), (33753, 33902), (31238, 31478), (17397, 17649), (26140, 26343), (33934, 34242), (17011, 17396), (28261, 28669)}, (6, 4, 2): {(34918, 35017), (17892, 17948), (27511, 27626), (16639, 16745), (27781, 28007), (8050, 8132), (868, 1274), (21768, 21853), (21902, 21999), (35360, 35999)}, (2, 0, 3, 0): {(26385, 26490), (3970, 4080), (8848, 8995), (24491, 24646), (3474, 3689)}, (3, 5, 7): {(13921, 13987), (11675, 11768), (15810, 15934), (8655, 8767), (15662, 15809), (15531, 15640)}, (6, 4, 2, 0): {(28877, 28963), (12902, 13020), (29148, 29443), (29778, 30000), (10596, 10732), (25123, 25363), (6270, 6387)}}, [((3, 0, 3, 5, 7), 258.2174556582623), ((6, 4, 2, 0, 1), 244.66131189425528), ((3, 5, 7, 6, 7, 1), 143.21532082396376), ((6, 4, 2, 0), 58.69834327256427), ((6, 4, 2), 52.10770136332576), ((2, 0, 3, 0), 35.73514257141529), ((5, 7, 6), 29.0868353301616), ((3, 5, 7), 20.236870397263232), ((2, 0, 1), 7.334740476609435)]), (0.8, 3): ({(1, 6, 4, 2, 0): {(28733, 28963), (27331, 27626), (18231, 18502), (16393, 16745), (21853, 22099)}, (5, 7, 6): {(31804, 31973), (33528, 33662), (25904, 26000), (31597, 31779), (22100, 22260), (24681, 25065), (11800, 12035), (13940, 14012), (20406, 20626), (4872, 5245), (3690, 3969), (1338, 2058), (30778, 31347)}, (6, 4, 2): {(34243, 34344), (34918, 35017), (27781, 28007)}, (3, 0, 3, 5, 7): {(15296, 15451), (21169, 21767), (14034, 14319), (25403, 25903), (26891, 27196), (8962, 9332), (6862, 7224), (9825, 10159), (19251, 19503), (26491, 26890), (29445, 29615), (29002, 29147), (22510, 24277), (11580, 11674), (20152, 20405)}, (2, 0, 1): {(2074, 2560), (14795, 14964), (31348, 31478), (7525, 7664)}, (3, 0, 3): {(4040, 4102), (11253, 11362), (18805, 19050), (24605, 24667)}, (6, 4, 2, 0, 1): {(381, 605), (10596, 10875), (606, 822), (8050, 8551), (34541, 34737), (33753, 33902), (32183, 32431), (6270, 6515), (17397, 17649), (12482, 12901), (33934, 34242), (17011, 17396), (28261, 28669)}, (1, 7, 6, 4, 2, 0): {(5355, 5552), (7716, 7974), (17731, 18030), (26001, 26260)}, (3, 5, 7): {(4113, 4840), (13088, 13243), (30157, 30580), (15642, 15809), (2759, 3325), (11675, 11768), (8656, 8767), (15810, 15934), (9629, 9791), (18031, 18126), (15531, 15641), (15935, 16287), (13463, 13593)}, (6, 4, 2, 0): {(8768, 8961), (29148, 29444), (17650, 17730), (12902, 13020), (868, 1337), (29778, 30002)}}, [((3, 0, 3, 5, 7), 201.96481349193718), ((6, 4, 2, 0, 1), 182.90796715015117), ((3, 5, 7), 64.8745306099196), ((1, 7, 6, 4, 2, 0), 61.263184278663324), ((1, 6, 4, 2, 0), 60.79410368439916), ((5, 7, 6), 56.437633424982636), ((6, 4, 2, 0), 49.87701952474127), ((2, 0, 1), 12.571292060118571), ((3, 0, 3), 11.213416529810699), ((6, 4, 2), 9.044935293804153)])}
    percTuples = [(0,1, None)]
    gammaVals = [0.99, 0.9, 0.7]
    for startPerc, endPerc, motifSubset in percTuples:
        for i, gamma in enumerate(gammaVals):
            plt.figure(i)
            plt.title("Gamma: %s" % gamma)
            motifValues, rankedList = driver6Motifs[(gamma, 3)]
            print(rankedList)
            motifColorKeys = motifValues.keys()
            motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
            create3Plots(startPerc, endPerc, 6, gamma, motifValues, motifColorMap, rankedList, motifSubset=motifSubset)
        plt.show()

    # anything past this is old

    # percTuples = [(0,1,None), (0.94, 1, None)]
    # driver6Motifs = {(0.99, 2): ({(0, 1): [(17609, 17649), (32529, 32748), (5548, 5627), (32236, 32431), (17717, 17891), (10702, 10872), (6350, 6515), (686, 822), (14835, 14964), (34120, 34242), (34658, 34736), (12986, 13088), (7967, 8049), (6516, 6861), (18479, 18804), (33829, 33902), (28941, 29005), (17101, 17396)], (6, 4): [(28877, 28906), (32183, 32235), (27665, 27682), (22293, 22425)], (3, 0, 3, 5, 7): [(15295, 15451), (21169, 21767), (25404, 25903)], (7, 6): [(33444, 33508), (33599, 33662), (13572, 13629), (9768, 9824), (31724, 31779), (11363, 11415), (7665, 7715), (18180, 18230), (11749, 11797), (13964, 14012), (23124, 24400), (30206, 30650), (11798, 12035), (1872, 2058), (12036, 12221), (29552, 29734), (20449, 20626), (25066, 25234), (9278, 9439), (27176, 27330), (7167, 7317), (3240, 3386), (10105, 10240), (16238, 16371), (24943, 25065), (19425, 19546), (22149, 22260), (14261, 14365), (18077, 18179), (31872, 31973), (3387, 3473), (7422, 7506), (33039, 33114), (25931, 26000), (4804, 4871), (5178, 5245)], (5, 7): [(5277, 5354), (4113, 4508), (30057, 30107), (26592, 26890), (20207, 20405), (19074, 19250)], (3, 5, 7): [(15531, 15641), (15810, 15934)], (2, 0): [(3474, 3689), (9440, 9629), (24491, 24680), (2074, 2244), (25235, 25403), (3970, 4112), (7525, 7611), (31348, 31430)], (1, 6): [(31431, 31495), (28389, 28693), (34344, 34387), (13704, 13848), (16746, 16885), (26261, 26353)], (4, 2): [(14778, 14834), (35018, 35301), (70, 251), (34782, 34917), (28261, 28388), (0, 69)], (0, 3): [(29374, 29511), (17949, 18076), (11268, 11362), (29006, 29093), (14066, 14144), (11583, 11647), (6935, 6992), (20140, 20194), (9873, 9913), (18874, 19050), (22448, 22621), (26408, 26564)], (6, 4, 2): [(868, 1274), (27781, 28007), (29148, 29373), (17397, 17608), (33934, 34119), (29778, 29948), (382, 539), (8768, 8901), (26140, 26260), (34541, 34657), (27511, 27626), (16639, 16745), (10596, 10701), (34243, 34343), (34918, 35017), (7867, 7966), (21903, 21999), (12482, 12576), (17011, 17100), (21768, 21852), (12902, 12985), (18397, 18478), (606, 685), (6270, 6349), (33753, 33828), (5477, 5547), (8050, 8120), (17650, 17716), (17892, 17948)], (3, 5): [(2759, 3239), (19251, 19424), (13089, 13242), (13314, 13462), (13463, 13571), (15935, 16237), (11675, 11748), (26891, 27175)]}, [((6, 4, 2), 50.32856449494014), ((3, 0, 3, 5, 7), 13.287885336665505), ((3, 5, 7), -15.016964965961634), ((6, 4), -28.064231185606566), ((4, 2), -33.05603143084397), ((2, 0), -50.81545374765616), ((0, 3), -58.63065307577004), ((1, 6), -62.83897719610762), ((5, 7), -64.53536550866656), ((3, 5), -66.840494952178), ((0, 1), -120.3742952458299), ((7, 6), -219.69391889475844)]), (0.6, 2): ({(2, 0, 1): [(14959, 15089), (7525, 7693)], (6, 4): [(5628, 5745), (7463, 7506), (27665, 27682), (34367, 34387), (382, 419)], (2, 0, 3, 0): [(26384, 26490), (25235, 25401)], (4, 2, 0, 1): [(14778, 14958), (5490, 5627), (32487, 32748), (70, 381)], (3, 0, 3, 5, 7): [(9603, 9791), (29001, 29147), (26491, 26890), (18805, 19250), (11579, 11674), (25402, 25903), (20152, 20405), (15296, 15451)], (2, 1, 6, 4, 2, 0, 1, 7, 6, 4, 2, 0): [(28278, 28964), (12222, 13020), (16886, 17612), (420, 1337)], (3, 0): [(22510, 22634), (29444, 29477), (22074, 22099), (19251, 19298), (20064, 20151)], (3, 0, 3, 5, 7, 6, 1): [(6860, 7421), (32942, 33443), (9825, 10241), (18002, 18220), (26891, 27510), (2758, 3387), (14026, 14749), (13087, 13321)], (5, 7, 6, 7, 1, 7, 6): [(19299, 19925), (22635, 24492), (29499, 29836), (11707, 12036), (20406, 21043), (15935, 16684), (24672, 25234), (4882, 5489), (13512, 13849)], (5, 7): [(5746, 6256), (15868, 15934)], (3, 5, 7): [(1338, 1944), (15531, 15640), (15662, 15809), (13921, 13987)], (6, 4, 2, 0): [(6270, 6387), (35360, 35999), (22293, 22509), (29148, 29443), (8767, 8961)], (0, 3, 0, 5, 7): [(31541, 31760), (3624, 3933)], (7, 6): [(33444, 33509), (22149, 22260), (25931, 26000), (28072, 28274), (11363, 11415), (31871, 31973)], (1, 6, 4, 2, 0): [(7740, 7974), (21855, 22073), (33712, 33886)], (4, 2): [(16685, 16745), (29837, 29949), (0, 69), (21816, 21854), (34294, 34345), (9422, 9558), (35018, 35302)], (1, 6, 4, 2, 0, 1): [(31974, 32431), (26001, 26343), (18231, 18804)], (6, 4, 2, 0, 1): [(10596, 11252), (33934, 34243), (17650, 17891), (8050, 8660), (34541, 34737), (1945, 2563), (31237, 31478)], (6, 4, 2): [(27780, 28007), (17892, 17948), (7694, 7739), (34918, 35017), (34779, 34917), (27511, 27626)], (3, 0, 3, 0, 5, 7, 5, 7, 6): [(30001, 30650), (8962, 9421), (4040, 4871), (21169, 21815)], (0, 3): [(24637, 24667), (24550, 24636), (11290, 11362)]}, [((2, 1, 6, 4, 2, 0, 1, 7, 6, 4, 2, 0), 156.84913472388843), ((3, 0, 3, 5, 7, 6, 1), 108.70817845581104), ((3, 0, 3, 0, 5, 7, 5, 7, 6), 89.7246838784228), ((5, 7, 6, 7, 1, 7, 6), 76.49475013218932), ((3, 0, 3, 5, 7), 53.93111517547479), ((6, 4, 2, 0, 1), 52.25253027191145), ((1, 6, 4, 2, 0, 1), 25.459636074203516), ((6, 4, 2, 0), 20.376044056648695), ((1, 6, 4, 2, 0), 17.31015438278169), ((0, 3, 0, 5, 7), 9.641692000285532), ((4, 2, 0, 1), 8.859106371672896), ((2, 0, 3, 0), 5.371350793385107), ((6, 4, 2), -10.332013885900778), ((2, 0, 1), -12.138987718397338), ((0, 3), -21.940982113558743), ((3, 5, 7), -24.466272468824926), ((5, 7), -25.62534278422474), ((3, 0), -31.460047284938), ((6, 4), -33.67528753049363), ((4, 2), -38.81037977302388), ((7, 6), -57.63185883016813)]), (0.7, 2): ({(6, 4): [(7463, 7506), (22293, 22425), (5628, 5745), (34367, 34387), (27665, 27682)], (0, 3, 5, 7, 6): [(9014, 9421), (29469, 29734), (13889, 14012), (18877, 19250), (2758, 3386)], (5, 7): [(11814, 11977), (30057, 30107), (15869, 15934), (5746, 6256)], (1, 7, 6, 4): [(255, 419), (34663, 34863)], (3, 0): [(15089, 15133), (25364, 25401), (30001, 30056), (26451, 26490), (20064, 20151)], (5, 7, 6, 4, 2, 0, 3): [(30779, 31437), (3690, 4071)], (6, 4, 2, 0): [(33753, 33886), (25123, 25363), (8767, 8961), (29148, 29443)], (1, 6): [(32290, 32455), (33510, 33662), (26261, 26353)], (4, 2): [(35018, 35302), (0, 69), (9422, 9557), (34294, 34345)], (6, 4, 2, 0, 1): [(10596, 11254), (18397, 18805), (6270, 6515), (17650, 17891), (8050, 8660), (33934, 34243), (5477, 5627), (1945, 2563)], (0, 3): [(11290, 11362), (24550, 24636), (8985, 9013)], (3, 0, 3, 5, 7): [(26491, 26890), (15296, 15451), (19251, 19502), (29001, 29147), (11579, 11675), (20152, 20405), (21169, 21767), (25402, 25903), (18002, 18126), (9603, 9791), (24647, 25014)], (6, 4, 2): [(27780, 28007), (34918, 35017), (35360, 35999), (21903, 21999), (27511, 27626), (21768, 21854), (16639, 16745), (17892, 17948), (7694, 7739)], (0, 3, 0, 5, 7, 6, 7): [(31541, 31802), (22000, 22292)], (2, 1, 6, 4, 2, 0, 1, 7, 6, 4, 2, 0): [(28278, 28963), (12222, 13020), (420, 1337)], (7, 6): [(28072, 28274), (18180, 18230), (25931, 26000), (12036, 12221), (33444, 33509), (3387, 3473), (31871, 31973), (11363, 11415)], (1, 7, 6, 4, 2, 0): [(31974, 32289), (26001, 26260), (16886, 17106), (34388, 34662), (17107, 17612), (29735, 30000), (7740, 7974)], (3, 5, 7): [(1338, 1944), (15662, 15809), (11676, 11768), (30157, 30580), (15531, 15640)], (5, 7, 6, 7, 1): [(20406, 21122), (33029, 33443), (15935, 16638), (4883, 5476), (13513, 13803)], (2, 0): [(3474, 3689), (26384, 26450)], (2, 0, 1): [(7525, 7693), (19761, 20063), (32488, 32748), (14795, 14958)], (3, 0, 3, 5, 7, 6, 7, 1): [(22510, 24492), (26891, 27510), (13086, 13321), (4081, 4882), (9825, 10513), (6861, 7421), (14026, 14749)], (3, 5): [(8661, 8747), (32758, 32941)]}, [((3, 0, 3, 5, 7, 6, 7, 1), 116.72180250507864), ((2, 1, 6, 4, 2, 0, 1, 7, 6, 4, 2, 0), 115.88837772587851), ((3, 0, 3, 5, 7), 77.72462692042664), ((1, 7, 6, 4, 2, 0), 77.70371008407545), ((6, 4, 2, 0, 1), 61.88390778932501), ((5, 7, 6, 4, 2, 0, 3), 28.371144913390342), ((0, 3, 0, 5, 7, 6, 7), 24.819760024284296), ((0, 3, 5, 7, 6), 22.374728029690495), ((6, 4, 2, 0), 14.45525077507617), ((5, 7, 6, 7, 1), -2.2607379575449182), ((1, 7, 6, 4), -4.409557193757577), ((6, 4, 2), -8.09764578140814), ((2, 0), -17.592942803715545), ((2, 0, 1), -18.69918248784475), ((0, 3), -22.472234753166887), ((3, 5), -22.540608546609572), ((4, 2), -26.667005592525335), ((3, 5, 7), -29.028092429811387), ((3, 0), -32.34546835095158), ((6, 4), -33.69841437750686), ((1, 6), -35.71436967966409), ((5, 7), -45.64429192820205), ((7, 6), -72.37126785904323)])}
    # gammaVals = [0.99, 0.7, 0.6]
    # for startPerc, endPerc, motifSubset in percTuples:
    #     for gamma in gammaVals:
    #         motifValues, rankedList = driver6Motifs[(gamma, 2)]
    #         print(rankedList)
    #         motifColorKeys = motifValues.keys()
    #         motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
    #         create3Plots(startPerc, endPerc, 6, gamma, motifValues, motifColorMap, rankedList, motifSubset=motifSubset)

    #0.99 turns or 0.6 -> 4,5,6,2
    # (0.99, 1), (0.74,0.772)
    # percTuples = [(0,1,None), (0.05, 0.135,None), (0.95,1, None), (0.857,0.865, None)]
    # percTuples = [(0.857,0.865, None)]
    # driver1Motifs = {(0.99, 2): ({(5, 6, 2): [(3662, 3954), (4509, 4786), (9702, 9965), (34759, 34999), (2717, 2941), (30870, 31069), (11751, 11917), (19812, 19972), (19277, 19425), (29643, 29786), (8274, 8415), (8682, 8815), (28356, 28487), (15999, 16120), (5007, 5118), (6734, 6835), (29301, 29400), (29110, 29204), (28717, 28811), (35714, 35805), (13656, 13746), (18381, 18468), (2421, 2500), (14086, 14159), (18797, 18866)], (5, 6): [(28041, 28136), (2135, 2216), (18012, 18084), (18914, 18979), (24533, 24583), (29072, 29109), (27524, 27705), (1924, 2051)], (6, 2): [(10452, 10646), (19123, 19222), (28155, 28290), (6910, 6997)], (2, 3): [(30121, 30691), (16311, 16588), (5209, 6463), (16589, 17012)], (4, 5): [(5163, 5208), (21887, 21919), (15170, 15193), (0, 1162), (10248, 10367), (30692, 30780), (18170, 18233)], (7, 2, 0): [(8816, 9444), (19973, 20319), (27706, 28040), (28812, 29045), (14866, 15043)], (2, 3, 4): [(20320, 21122), (17013, 17683), (12751, 13006), (10892, 11141), (14245, 14466), (13007, 13180)], (2, 0, 2, 3): [(21251, 21886), (12348, 12571), (7337, 8032)], (3, 4): [(33033, 34758), (31343, 33032), (35000, 35435), (21920, 21989), (15262, 15552), (4787, 5006), (31070, 31289), (10647, 10861), (35436, 35548), (3955, 4508), (3144, 3661), (31290, 31342)]}, [((5, 6, 2), 63.19690085799703), ((2, 0, 2, 3), -6.9573333810825915), ((7, 2, 0), -22.724035430579328), ((6, 2), -31.020524417130837), ((2, 3, 4), -32.02289714247833), ((5, 6), -39.37536273292384), ((4, 5), -40.96860044781156), ((2, 3), -43.07322583516956), ((3, 4), -109.46139566622152)]), (0.6, 2): ({(5, 6, 2): [(29110, 29204), (5008, 5118), (3667, 3954), (2135, 2216)], (7, 4, 5, 6, 2): [(18469, 18866), (15553, 16120), (19223, 19425), (1500, 2051), (18980, 19222), (19766, 19972), (11440, 11917), (18234, 18468), (0, 1399)], (4, 5, 6, 2, 3): [(30867, 31247), (10447, 10839)], (2, 0, 2, 3): [(16311, 16587), (10862, 11100), (12662, 12967), (21251, 21886), (14924, 15169), (34846, 35417)], (5, 6): [(29072, 29109), (24533, 24583), (18012, 18084), (34759, 34845), (28041, 28136)], (4, 5): [(18913, 18941), (21887, 21919), (18170, 18233), (5163, 5208), (15170, 15193)], (3, 4, 5, 6, 2, 3, 4): [(9692, 10317), (3955, 5007)], (7, 4, 5, 6, 2, 7, 0, 7, 4, 5, 6, 2, 0, 2, 3, 4): [(29205, 30760), (13336, 14466), (2321, 3666)], (7, 2, 0): [(6957, 7262), (7292, 7546), (19973, 20319)], (2, 3, 4): [(33033, 34758), (31343, 33032), (20320, 21122), (5209, 6465), (13008, 13180), (16588, 17683), (15237, 15552)], (0, 4, 5, 6, 2, 0, 4, 7, 5, 6, 2, 7, 2, 0): [(8032, 9444), (28291, 29045)], (3, 4): [(14467, 14833), (22376, 22709), (35436, 35548), (11142, 11374), (21920, 21989), (31290, 31342)], (5, 6, 2, 7): [(6734, 6909), (27524, 27757), (35714, 35892)]}, [((7, 4, 5, 6, 2, 7, 0, 7, 4, 5, 6, 2, 0, 2, 3, 4), 171.01476208938647), ((0, 4, 5, 6, 2, 0, 4, 7, 5, 6, 2, 7, 2, 0), 98.22130182080717), ((7, 4, 5, 6, 2), 83.73189122816916), ((3, 4, 5, 6, 2, 3, 4), 22.580121290470125), ((4, 5, 6, 2, 3), 8.846386674390233), ((5, 6, 2, 7), 7.513437795331147), ((5, 6, 2), -4.229518438890164), ((2, 0, 2, 3), -5.3430998224435875), ((7, 2, 0), -16.697948608525365), ((5, 6), -28.862298587288286), ((4, 5), -32.060364179843795), ((2, 3, 4), -35.43745532855345), ((3, 4), -63.192973130297005)]), (0.7, 2): ({(5, 6, 2): [(6734, 6835), (18797, 18866), (29110, 29204)], (7, 4, 5, 6, 2, 0): [(15553, 16124), (8517, 8843), (11440, 12347)], (2, 0, 2, 3): [(16311, 16587), (34846, 35417), (7337, 8032), (21251, 21886)], (5, 6): [(4509, 4702), (34759, 34845), (24533, 24583), (27524, 27706), (18932, 18979), (29072, 29109), (18012, 18084), (3667, 3778), (28041, 28136)], (4, 5): [(21887, 21919), (15170, 15193), (10248, 10367), (18170, 18233), (5163, 5208)], (7, 4, 5, 6, 2, 7, 0, 7, 4, 5, 6, 2, 0, 2, 3, 4): [(13337, 14466), (29205, 30760), (2321, 3666)], (7, 2, 0): [(6957, 7262), (28812, 29045), (14866, 15043), (19973, 20319)], (2, 3): [(4703, 4962), (12411, 12571)], (2, 3, 4): [(33033, 34758), (31343, 33032), (12751, 13007), (3779, 4508), (15237, 15552), (10892, 11141), (5209, 6465), (10567, 10861), (16588, 17683), (13008, 13180), (20320, 21122)], (5, 6, 2, 7): [(18381, 18796), (2135, 2307)], (3, 4): [(35436, 35548), (21920, 21989), (31290, 31342), (14467, 14833), (22376, 22709), (31070, 31289)], (4, 5, 6, 2): [(4963, 5118), (35699, 35806), (19813, 19972), (0, 1399), (8247, 8415), (28326, 28487), (9697, 9965), (19267, 19425), (19042, 19222), (1504, 2051), (28716, 28811), (30867, 31069)]}, [((7, 4, 5, 6, 2, 7, 0, 7, 4, 5, 6, 2, 0, 2, 3, 4), 171.3990492123976), ((4, 5, 6, 2), 69.15583341346237), ((7, 4, 5, 6, 2, 0), 33.164539717128555), ((5, 6, 2, 7), 3.5804867312318596), ((5, 6, 2), -4.824351630761974), ((2, 0, 2, 3), -7.04121302196384), ((7, 2, 0), -19.924825869422435), ((2, 3), -24.3321201396114), ((4, 5), -31.938041404391768), ((5, 6), -41.176770312342995), ((2, 3, 4), -45.7385420747872), ((3, 4), -63.20774973072246)])}
    # gammaVals = [0.99]
    # for startPerc, endPerc, motifSubset in percTuples:
    #     for gamma in gammaVals:
    #         motifValues, rankedList = driver1Motifs[(gamma, 2)]
    #         print(rankedList)
    #         motifColorKeys = motifValues.keys()
    #         motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
    #         create3Plots(startPerc, endPerc, 1, gamma, motifValues, motifColorMap, rankedList, motifSubset=motifSubset)


    # percTuples = [(0,1,None),(0, 0.03, None), (0.07, 0.1, None)]
    # driver2Motifs = {(0.99, 2): ({(1, 6, 3, 5): [(1383, 1964), (996, 1382)], (5, 4, 3, 4, 3): [(2961, 3222), (2775, 2902)], (1, 5, 1, 3, 5): [(460, 705), (2139, 2350)]}, [((1, 5, 1, 3, 5), 41.034167780174215), ((5, 4, 3, 4, 3), 37.61599983078587), ((1, 6, 3, 5), 25.452823433601697)]), (0.6, 2): ({(1, 6, 3, 5): [(1383, 1964), (996, 1382)], (5, 4, 3, 4, 3): [(157, 387), (2775, 2902), (2961, 3222)], (1, 5, 1, 3, 5): [(460, 705), (2139, 2350)]}, [((5, 4, 3, 4, 3), 58.782753078646444), ((1, 5, 1, 3, 5), 40.91124561581006), ((1, 6, 3, 5), 25.329901269237528)]), (0.7, 2): ({(1, 6, 3, 5): [(1383, 1964), (996, 1382)], (5, 4, 3, 4, 3): [(2775, 2902), (2961, 3222)], (1, 5, 1, 3, 5): [(460, 705), (2139, 2350)]}, [((1, 5, 1, 3, 5), 41.034167780174215), ((5, 4, 3, 4, 3), 37.61599983078587), ((1, 6, 3, 5), 25.452823433601697)])}
    # gammaVals = [0.99, 0.7, 0.6]
    # for startPerc, endPerc, motifSubset in percTuples:
    #     for gamma in gammaVals:
    #         motifValues, rankedList = driver2Motifs[(gamma, 2)]
    #         print(rankedList)
    #         motifColorKeys = motifValues.keys()
    #         motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
    #         create3Plots(startPerc, endPerc, 2, gamma, motifValues, motifColorMap, rankedList, motifSubset=motifSubset)

    # percTuples = [(0,1, None), (0.477, 0.586, None), (0.01, 0.042, None), (0.375, 0.399, None)]
    # percTuples = [(0.49, 0.54, [(3,7,6,3,7)]), (0.375, 0.395, [(3,7,6,3,7)])]
    # driver3Motifs = {(0.99, 2): ({(6, 4): [(2485, 2628), (5207, 5322), (21411, 21518), (3669, 3772), (742, 834), (20971, 21052), (20869, 20947), (22158, 22231), (2417, 2484), (17675, 17740), (22384, 22448), (14402, 14452), (14273, 14322), (21363, 21410), (197, 243), (15003, 15042)], (0, 3, 6): [(1024, 1113), (17741, 18052), (928, 1023)], (4, 5): [(16211, 16285), (2729, 2774), (1634, 2144), (16904, 17233), (16286, 16599), (21662, 21946), (15694, 15950), (3373, 3598), (4328, 4468), (3252, 3372), (15093, 15174)], (3, 7): [(18151, 20736), (7969, 9054), (6986, 7522), (15175, 15693), (6540, 6831), (13760, 14015), (5843, 6084), (17234, 17396), (13595, 13736), (20737, 20868), (17460, 17584), (2997, 3119), (14016, 14135), (14168, 14272), (5518, 5619), (18053, 18150), (6446, 6539), (6899, 6985), (16600, 16685), (1271, 1349), (6832, 6898), (677, 741)], (3, 2, 6): [(9220, 13366), (2145, 2330), (2900, 2984), (14818, 14953)], (4, 0): [(275, 398), (3131, 3251), (16128, 16210), (601, 676), (5152, 5206)]}, [((0, 3, 6), -6.457397153096583), ((3, 2, 6), -15.212207832528973), ((4, 0), -21.977448123608802), ((6, 4), -37.53155199761778), ((4, 5), -44.29949260707707), ((3, 7), -146.36105908582365)]), (0.6, 2): ({(6, 4): [(17675, 17740), (15003, 15042), (164, 243), (2417, 2485)], (6, 4, 0, 3): [(600, 721), (742, 884), (272, 438), (5207, 5399)], (7, 6, 0, 4, 5, 3): [(7916, 9179), (16666, 17320), (1370, 2236), (15266, 16115)], (4, 5, 3): [(4328, 4529), (15093, 15265), (2729, 2941)], (2, 7, 2): [(9180, 13340), (5444, 5677), (5787, 6445), (5693, 5749), (18220, 20726)], (4, 0): [(5152, 5206), (14319, 14372), (3131, 3251), (22417, 22620), (14439, 14763), (16128, 16210), (4705, 5015)], (3, 7, 6): [(20756, 20914), (1271, 1369)], (4, 5): [(3252, 3372), (16286, 16597), (16211, 16285)], (4, 3, 6, 4, 5): [(21407, 21584), (21585, 21945), (2511, 2728), (20915, 21292), (22084, 22283), (3597, 4225)], (3, 7, 1): [(22621, 35999), (6446, 7915)], (3, 7, 6, 3, 7): [(17772, 18161), (14024, 14277), (13605, 14023)]}, [((4, 3, 6, 4, 5), 78.5152903536966), ((7, 6, 0, 4, 5, 3), 57.64697303073657), ((6, 4, 0, 3), 24.56528794556924), ((3, 7, 6, 3, 7), 11.81762812209031), ((4, 5, 3), -3.8622300020648948), ((3, 7, 6), -10.786715603097706), ((3, 7, 1), -18.27158153298134), ((4, 5), -19.527317240939972), ((6, 4), -20.891849389808826), ((4, 0), -25.496810612218667), ((2, 7, 2), -33.16527287346045)]), (0.7, 2): ({(2, 4, 0, 4, 5, 4, 5, 3): [(3120, 3668), (9272, 13504), (16116, 16668)], (6, 4): [(14402, 14452), (196, 243), (21363, 21410), (15003, 15042), (742, 834), (20869, 20950), (2417, 2484), (17675, 17740), (2485, 2628)], (5, 4, 5, 3): [(4245, 4529), (2629, 2840)], (0, 3, 6): [(4705, 5057), (14453, 14864), (22055, 22183), (17741, 18052), (928, 1023), (4530, 4704), (1024, 1113)], (4, 5): [(15694, 15950), (15093, 15174)], (3, 7): [(14168, 14272), (5517, 5620), (18151, 20736), (7969, 9057), (6832, 6985), (14016, 14135), (6446, 6539), (13594, 13736), (86, 160), (1271, 1349), (20737, 20868), (6540, 6831), (677, 741), (6986, 7522), (2997, 3119), (15175, 15693), (17460, 17584), (13759, 14015), (5843, 6084)], (6, 0, 4, 5, 3): [(1442, 2236), (3669, 4244), (16686, 17322), (21650, 22018), (20971, 21351), (21411, 21649)], (6, 4, 0, 3): [(5207, 5443), (22384, 22736), (272, 464), (14273, 14401)], (4, 0): [(5152, 5206), (601, 676)]}, [((2, 4, 0, 4, 5, 4, 5, 3), 82.8894671744807), ((6, 0, 4, 5, 3), 68.9395465027618), ((6, 4, 0, 3), 23.02124578878652), ((5, 4, 5, 3), 4.353544067932783), ((0, 3, 6), -2.795074851092832), ((4, 0), -12.32595756827197), ((4, 5), -14.797463779352288), ((6, 4), -30.906766345980685), ((3, 7), -131.66445300654377)])}
    # gammaVals = [0.6]
    # # gammaVals = [0.99, 0.6, 0.7]
    # for startPerc, endPerc, motifSubset in percTuples:
    #     for gamma in gammaVals:
    #         motifValues, rankedList = driver3Motifs[(gamma, 2)]
    #         motifColorKeys = motifValues.keys()
    #         motifColorMap = {k:colormap[i] for i,k in enumerate(motifColorKeys)}
    #         create3Plots(startPerc, endPerc, 3, gamma, motifValues, motifColorMap, rankedList, motifSubset=motifSubset)


def create3Plots(startPerc, endPerc, driver_id, gamma, motifValues, motifColorMap, rankedList, motifSubset=None):
    if motifSubset != None:
        tempValues = {}
        for m in motifSubset:
            tempValues[m] = motifValues[m]
        motifValues = tempValues
    output_folder = "%s/%s/%s/%s_%s" % (outputDir, driver_id, gamma, startPerc, endPerc)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    clusterColorMap = colormap[:8]
    lats, longs = getNavData("%s/raw_data/%s/latlong.csv" % (resultsDir, driver_id))
    startVal = int(startPerc * len(lats))
    endVal = int(endPerc * len(lats))

    long_min = min(longs[startVal:endVal]) - 0.007
    long_max = max(longs[startVal:endVal]) + 0.002
    long_center = (long_min + long_max)/2
    lat_min = min(lats[startVal:endVal]) - 0.0001
    lat_max = max(lats[startVal:endVal]) + 0.002
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
    
    offset = 0.001
    #plot old assignments 
    '''
    old_assignment_name = "%s/old_assignments/%s_clust8_beta100.0.out" % (resultsDir,driver_id)
    assignMapLats, assignMapLongs = getClusterData(old_assignment_name, lats, longs, startVal, endVal, offset)
    plotClusterValues(m, assignMapLats, assignMapLongs, clusterColorMap)
    # if SHOW: plt.show()
    if SAVE: plt.savefig("%s/old_assignment.png"% output_folder)
    # plt.clf()

    lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,mode="expand", borderaxespad=0., ncol=8)

    # plot new assignments: 0.6
    new_assignment_name = "%s/new_assignments/%s/_clust8_beta100.0_gamma%s_req3.out" % (resultsDir, driver_id, gamma)
    assignMapLats, assignMapLongs = getClusterData(new_assignment_name, lats, longs, startVal, endVal, offset*2)
    plotClusterValues(m, assignMapLats, assignMapLongs, clusterColorMap)
    # if SHOW: plt.show()
    if SAVE: plt.savefig("%s/new_assignment.png"% output_folder)
    # plt.clf()

    '''
    # plot motifs    
    motifLats, motifLongs, endPointLats, endPointLongs=getMotifData(lats, longs, startVal, endVal, motifValues)    
    lgd = plotMotifValues(m, motifLats, motifLongs, endPointLats, endPointLongs, motifColorMap, rankedList)
    if SAVE:
        art = []
        art.append(lgd)
        plt.savefig("%s/new_motif.png"% output_folder, additional_artists=art, bbox_inches="tight")
    # if SHOW: plt.show()
    # plt.clf()

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

def getClusterData(input_name, lats, longs, startVal, endVal, longOffset):
    assignMapLongs = defaultdict(list)
    assignMapLats = defaultdict(list)
    with open(input_name, 'r') as instream:
        lines = instream.readlines()
        assigns = [int(val.strip()) for val in lines]
        for i,a in enumerate(assigns):
            if i >= startVal and i <= endVal:
                assignMapLongs[a].append(longs[i] - longOffset)
                assignMapLats[a].append(lats[i] + longOffset)
    return assignMapLats, assignMapLongs

def plotClusterValues(m, assignMapLats, assignMapLongs, clusterColorMap):
    keys = assignMapLongs.keys()
    plt.figure(1)
    for k in keys:
        x, y = m(assignMapLongs[k], assignMapLats[k])
        m.scatter(x,y,20,marker='.',color=clusterColorMap[k], label="%s" % k)
    lgd = None
    # lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,mode="expand", borderaxespad=0., ncol=8)
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

def plotMotifValues(m, motifLats, motifLongs, endPointLats, endPointLongs, motifColorMap, rankedList):
    keys = motifLats.keys()
    rankMap = {m:s for m,s in rankedList}
    for i,k in enumerate(keys):
        x, y = m(motifLongs[k], motifLats[k])
        endX, endY = m(endPointLongs[k], endPointLats[k])
        score = rankMap[k]
        motifString = "%s: %.2f" % (list(k), score)
        m.scatter(x,y,20,marker='.',color=motifColorMap[k], label=motifString)
        # m.scatter(endX, endY, 60, marker='|',color='black')
    numCols = min(len(keys), 2)
    lgd = None
    lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,mode="expand", borderaxespad=0., ncol=4)

    # lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.axes().set_aspect('equal', 'datalim')
    return lgd

run()