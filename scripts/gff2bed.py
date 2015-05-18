__author__ = 'maplesod'

import sys

with open(sys.argv[1]) as f:
    content = f.readline()

    print "track name=\"junctions\""

    index = 0;
    for line in f:

        if not line.startswith("#"):
            words = line.split()

            if words[2] == "intron":

                start = int(words[3]) - 1
                end = int(words[4])

                dif = end - start

                print words[0] + "\t" + str(start) + "\t" + str(end) + "\tjunc_" + str(index) + "\t0\t" + words[6] + "\t" + str(start) + "\t" + str(end) + "\t255,0,0\t2\t0,0\t0," + str(dif)
                index += 1