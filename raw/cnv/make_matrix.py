import sys

fil=open("genes.lst")
reacts = [i.strip() for i in fil.readlines()]
fil.close()

fil=open("cancer.lst")
cancer = [i.strip() for i in fil.readlines()]
fil.close()


matrix = {}
fil=open(sys.argv[1])
for lin in fil.readlines():
    flds = lin.strip().split(",")
    canc = flds[6]
    matrix[flds[6]+"\t"+flds[1]] = flds[-1]
fil.close()

print "Reaction\t",
for c in cancer:
    print c+"\t",
print
for r in reacts:
    found = []
    z=0
    found.append(r)
    for c in cancer:
        if(c+"\t"+r in matrix):
            found.append(matrix[c+"\t"+r])
            if(float(matrix[c+"\t"+r]) != 0.0):
                z=1
        else:
            found.append("1.000")
    if(z==1):
        print "\t".join(found)

