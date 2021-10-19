import sys
inFile = sys.argv[1]
outFile = sys.argv[2]

with open(inFile,'r') as i:
    lines = i.readlines()

with open(outFile,'w') as o:
    for line in lines:
        if line.find("#") == -1 and line != "\n": 
            o.write(line)
        