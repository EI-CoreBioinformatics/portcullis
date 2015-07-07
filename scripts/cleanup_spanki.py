__author__ = 'maplesod'

import sys

with open(sys.argv[1]) as f:
    content = f.readline()

    print ("refid\tstart\tend\tstrand\tdinucleotide\tintron_size\tannostatus\tgmcode\tregcode\tgeneassign\tgeneassignL\tgeneassignR\tunfilt_cov\tcov\tnormcov\toffsets\tentropy\thamming3\thamming5\tMAXmmes\tMAXminanc\tlirt\trirt\tirt\tdatrans\tdncov\tancov")

    for line in f:
        
        if not line == "juncid":

            words = line.split("\t")
            
            parts1 = words[0].split(":")
            parts2 = parts1[1].split("_")

            print (parts1[0] + "\t" + parts2[0] + "\t" + parts2[1] + "\t" + parts1[2], end="")

            for i in range(1, len(words)):
                print ("\t" + words[i], end="")

            print()
