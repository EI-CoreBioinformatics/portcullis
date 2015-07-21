__author__ = 'maplesod'

import sys

class Junction:

    seq = ""
    start = 0
    end = 0
    strand = "."
    cov = 0

    def __init__(self):
        self.data = []

    def display(self):
        print (self.seq + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + self.strand + "\t" + str(self.cov))



with open(sys.argv[1]) as f:
    content = f.readline()

    junctions_set = {}
    junctions = list();

    index = 0;
    for line in f:

        if not line.startswith("#"):
            words = line.split()

            if words[2] == "intron":

                start = int(words[3])
                end = int(words[4])

                dif = end - start
                key = words[0] + "_" + str(start) + "_" + str(end) + "_" + words[6]
                line = words[0] + "\t" + str(start) + "\t" + str(end) + "\tjunc_" + str(index) + "\t0\t" + words[6] + "\t" + str(start) + "\t" + str(end) + "\t255,0,0\t2\t0,0\t0," + str(dif)

                j = Junction()
                j.seq = words[0]
                j.start = start
                j.end = end
                j.strand = words[6]

                if key not in junctions_set:
                    junctions_set[key] = line
                    junctions.append(j)
                    index += 1

# Sort by
junctions.sort(key=lambda x: x.end)
junctions.sort(key=lambda x: x.start)
junctions.sort(key=lambda x: x.seq)

print ("track name=\"junctions\"")

index = 0;
for j in junctions:
    print (j.seq + "\t" + str(j.start - 2) + "\t" + str(j.end + 1) + "\tjunc_" + str(index) + "\t" + str(j.cov) + "\t" + j.strand + "\t" + str(j.start - 1) + "\t" + str(j.end) + "\t255,0,0\t2\t1,1\t0," + str(j.end - j.start))

    index += 1
    