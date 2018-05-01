mapping1 = {3:7, 6:8, 2:9, 7:3, 8:6, 9:2}
mapping2 = {3:7, 1:8, 2:6, 7:3, 8:1, 6:2}
# mapping2 = {2:7, 9:2, 0:9, 5:6, 7:0, 6:5}
def transformMotif(m, mapping):
    result = []
    
    for i in m:
        if i in mapping:
            result.append(str(mapping[i]))
        else:
            result.append(str(int(i)))
    return result

def printResult(vals, mapping):
    for k in vals.keys():
        print("%s & %s\\\\" %(",".join(transformMotif(k, mapping)), len(vals[k])))