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

for c in cancer:
    wrf=open(c+".cnv","w")
    for r in reacts:
        if(c+"\t"+r in matrix):
            wrf.write(r+"\t"+matrix[c+"\t"+r]+"\n")
    wrf.close()
