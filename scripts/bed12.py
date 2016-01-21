#!/usr/bin/env python3
__author__ = 'maplesod'

# Can't use bedtools as bedtools doesn't properly support bed12 files... specifically we need to base our intersections on the
# thickstart and thickend columns

class BedEntry:

    chrom = ""
    start = 0
    end = 0
    name = ""
    score = 0
    strand = "?"
    thick_start = 0
    thick_end = 0
    red = 0
    green = 0
    blue = 0
    block_count = 0
    block_sizes = list()
    block_starts = list()

    def __init__(self):
        self.data = []

    def __str__(self):
        return self.chrom + "\t" + str(self.start) + "\t" + str(self.end)+ "\t" + self.name \
                + "\t" + str(self.score) + "\t" + (self.strand if not self.strand == "" else "?") + "\t" + str(self.thick_start) + "\t" + str(self.thick_end) \
                + "\t" + str(self.red) + "," + str(self.green) + "," + str(self.blue) \
                + "\t" + str(self.block_count) + "\t" + ",".join(str(x) for x in self.block_sizes) + "\t" + ",".join(str(x) for x in self.block_starts)

    def __cmp__(self, other):
        if hasattr(other, 'chrom') and hasattr(other, 'start') and hasattr(other, 'end'):
            if self.__lt__(other):
                return 1
            elif self.__gt__(other):
                return -1
            else:
                return 0

    def __key__(self):
        return self.chrom + "_" + str(self.thick_start) + "_" + str(self.thick_end)
    def __hash__(self):
        return hash(self.__key__())

    def __lt__(self, other):
        if self.chrom.__lt__(other.chrom):
            return True
        else:
            if self.chrom.__gt__(other.chrom):
                return False
            else:
                # The same chrom
                if self.thick_start < other.thick_start:
                    return True
                else:
                    if self.thick_start > other.thick_start:
                        return False
                    else:
                        # Also the same start
                        if self.thick_end < other.thick_end:
                            return True
                        else:
                            if self.thick_end > other.thick_end:
                                return False
                            else:
                                # Same!
                                return False


    def __eq__(self, other):
        return not self<other and not other<self
    def __ne__(self, other):
        return self<other or other<self
    def __gt__(self, other):
        return other<self
    def __ge__(self, other):
        return not self<other
    def __le__(self, other):
        return not other<self


    @staticmethod
    def create_from_key(key):

        b = BedEntry()

        parts = key.split("_")

        b.chrom = parts[0]
        b.thick_start = int(parts[1])
        b.thick_end = int(parts[2])

        if len(parts) == 4:
            b.strand = parts[3]

        b.name = "junc"
        b.score = 0
        b.start = b.thick_start - 1
        b.end = b.thick_end + 1
        b.red = 255
        b.green = 0
        b.blue = 0
        b.block_count = 2
        b.block_sizes = list()
        b.block_sizes.append(1)
        b.block_sizes.append(1)
        b.block_starts = list()
        b.block_starts.append(0)
        b.block_starts.append(b.thick_end - b.thick_start)

        return b

    @staticmethod
    def create_from_line(key):

        b = BedEntry()

        parts = key.split("\t")

        b.chrom = parts[0]
        b.start = int(parts[1])
        b.end = int(parts[2])
        b.name = parts[3]
        b.strand = parts[4]
        b.score = int(parts[5])
        b.thick_start = int(parts[6])
        b.thick_end = int(parts[7])

        c_parts = parts[8].split(",")
        b.red = int(c_parts[0])
        b.green = int(c_parts[1])
        b.blue = int(c_parts[2])
        b.block_count = int(parts[9])

        bsize_parts = parts[10].split(",")
        for b in bsize_parts:
            b.block_sizes.append(int(b))

        bstart_parts = parts[11].split(",")
        for b in bstart_parts:
            b.block_starts.append(int(b))

        return b


def makekey(line, usestrand, tophat) :
    words = line.split()
    overhang = words[10]
    overhang_parts = overhang.split(",")
    lo = int(overhang_parts[0])
    ro = int(overhang_parts[1])
    chr = words[0]
    start = str(int(words[6]) + lo) if tophat else words[6]
    end = str(int(words[7]) - ro) if tophat else str(int(words[7]))
    strand = words[5]
    key = chr + "_" + start + "_" + end
    if usestrand:
        key += "_" + strand
    return key

def makekeyfromtab(line, usestrand) :
    words = line.split()
    chr = words[2]
    start = words[4]
    end = words[5]
    strand = words[11]
    key = chr + "_" + start + "_" + end
    if usestrand:
        key += "_" + strand
    return key

def loadbed(filepath, usestrand, tophat) :

    index = 0
    items = set()

    with open(filepath) as f:
        # Skip header
        f.readline()
        for line in f:
            key = makekey(line, usestrand, tophat)
            items.add(key)
            index += 1
    if len(items) != index :
        print ("duplicated items in bed file " + filepath)
    return items



def filterbed(filepath, refset, mode, usestrand, tophat, outfile) :

    o = open(outfile, 'w')

    o.write("track name=\"junctions\"\n")

    index = 0
    with open(filepath) as i:

        i.readline()    # Skip header
        for line in i:

            key = makekey(line, usestrand, tophat)

            if mode == 0:
                if key in refset:
                    o.write(line)
            elif mode == 1:
                if key not in refset:
                    o.write(line)

            index += 1
    o.close()

def filtertab(filepath, outfile, bed, mode, usestrand) :

    o = open(outfile, 'w')

    index = 0;
    with open(filepath) as f:

        o.write(f.readline())

        for line in f:

            line = line.strip()
            if line != "":
                key = makekeyfromtab(line, usestrand)
                if mode == 0:
                    if key in bed:
                        o.write(line + "\n")
                elif mode == 1:
                    if key not in bed:
                        o.write(line + "\n")

    o.close()

