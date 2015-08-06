.. _metrics:

Metrics
=======

During the junction making stage, Portcullis generates a list of potential splice
junctions by assuming any RNAseq alignment containing an 'N' cigar operation.  For
all potential junctions portcullis makes a number of observations about how the
reads align around the region.  The complete list of observations (or metrics) are
described in this section.  Before we do that, here is a small glossary of terms
that we will use throughout:

* Splice junction - Splice junctions are points on a DNA sequence at which `superfluous' DNA is removed during the process of protein creation in higher organisms.  For the purposes of this tool a splice junction is essentially the same as an intron.
* Splice site - A junction has both a 5' donor and 3' acceptor site, which mark the start and end of the junction.  Both donor and acceptor sites are 2bp long, and usually contain a canonical motif `GT*AG`, or its reverse complement on the negative strand.
* Intron - Any nucleotide sequence within a gene that is removed by RNA splicing while the final mature RNA product of a gene is being generated.  For the purposes of this tool an intron is essentially the same as a splice junction.
* Anchor - The part of a spliced read that aligns to the genome.  This will represent regions both upstream and downstream of the splice junction.
* Multiple Spliced Reads (MSRs) - Reads that cover more than a single spliced junction.  These types of reads will become more common as sequencers become capable of producing longer reads.


Metric 1 - Canonical splice site
--------------------------------

Most splice junctions contain distinctive dinucleotide pattern at the start and 
end of the splice junction.  These patterns are called `donors` and `acceptors`
respectively.  If the potential junction has the typical `GT_AG` motif (or its
reverse complement on the negative strand then it is considered a canonical splice
site.  Alternatively, if the junction contains either `AT_AC` or `GC_AG` or their
reverse complements then it is considered a semi-canonical splice site.  Other motifs
at the donor and acceptor sites are considered novel splice sites.

#TODO pretty pic


Metric 2 - Number of spliced reads
----------------------------------

This a count of the number of reads containing an 'N' cigar operation at exactly
the genomic position represented by this junction.  Because we use 'N' cigar operations
to derive all potential junctions each junction will have a value of 1 or more for
this metric.

#TODO pretty pic


Metric 3 - Number of Distinct Alignments
----------------------------------------

This is the number of distinct / non-redundant spliced reads supporting the junction.
Therefore duplicate reads are only counted as a single read for this metric.


Metric 4 - Number of Reliable Alignments
----------------------------------------

This is the number of spliced reads supporting this junction that probably do
not map elsewhere in the genome.  For us we try to gauge this by treating all spliced
alignments with a mapping score of 30 or greater as reliable alignments.


Metric 5 - Intron size
----------------------

This is the length of the region defined by the 'N' cigar operation.  Which should
also be the size of the intron if this is a genuine splice site.


Metric 6 / 7 - Left / Right Anchor Size
---------------------------------------

The size of each anchor measured in transcript coordinates.  These metrics will 
therefore have a maximum value of the readlength - 1.  Upstream and downstream
junctions contained within MSRs associated to this junction are collapsed for 
the purposes of this metric.


Metric 8 - Max-min-Anchor
-------------------------

This is the maximum of the minimum anchor sizes of all spliced reads associated with
this junction.  Again, all upstream and downstream junctions contained within MSRs
are collapsed for the purposes of this metric.  

#TODO pretty pic


Metric 9 - Difference between Anchors
-------------------------------------

This is the difference between the minimum L/R anchor and the maximum L/R anchor
with the smaller of the two values reported. Again, all upstream and downstream 
junctions contained within MSRs are collapsed for the purposes of this metric.  

#TODO pretty pic


Metric 10 - Distinct Anchors
----------------------------

For this metric the number of distinct left side anchors and the number of distinct 
right anchors are both recorded.  The value used for the metric is the smaller of 
the two values reported.


Metric 11 - Entropy
-------------------

This describes the entropy of the spliced reads associated with this junction. 
Higher entropy is generally more indicative of a genuine junction than a lower score.
The entropy score is a function of both the total number of reads that map to a 
given junction and the number of different offsets to which  those reads map and 
the number that map at each offset. Thus, junctions with multiple reads mapping 
at each of the possible windows across the junction will be assigned a higher 
entropy score, than junctions where many reads map to only one or two positions. 
     
Entropy was calculated using the following equations::

* p_i = nb_reads_at_offset_i / total_reads_in_junction_window 
* Entropy = - sum_i(p_i * log(pi) / log2) 

#TODO pretty pic


Metric 12 - Maximum of the Minimal Match of Either Side (MaxMMES)
-----------------------------------------------------------------

This metric takes into account mismatches in the anchors on either side of the junction.
For each spliced read associated with the junction, we look at both anchors.  The
score for each anchor is the anchor length minus any mismatches to the reference.
The minimal score from either the upstream or downstream anchor is taken.  Then from
these scores the maximum is taken from all spliced reads.  Caution must be taken
interpreting this result.  Some splice aware mappers such as tophat by default require
all spliced reads to have 0 mismatches.  In this case this metric is not very useful.

#TODO pretty pic

Metric 13 / 14 - 5' and 3' Hamming distance
---------------------------------------------

Aligners can often make incorrect alignments around repeated genomic locations.
In these instances it is good to know whether the region on the on the left side
of the donor site and the left side of the acceptor site, in addition to the region
on the right side of the donor site and the right side of the acceptor site are
similar.  In this is the case then it is likely that the false splice alignments
have been made.  We record both figures in terms of the hamming distances between
the regions.  Low scores indicate similarity, and therefore high change of alignment
to a repeat region, high scores indicate difference and therefore low chance of alignment
to a repeat region.

#TODO pretty pic

Metric 15 - Unspliced Coverage around junction
----------------------------------------------

When considering unspliced reads around a junction site, you would typically expect
to see a tailing off of reads towards the 5' junction boundary, and a ramping up
after the 3' junction boundary.  However, in practice this is complicated by MSRs,
alternative splicing and junctions near sequence ends.

#TODO pretty pic

Metric 16 - Unique Junction
---------------------------

This boolean metric determines whether or not there are any other junctions within
this junctions region.  In particular, whether any other junctions share it's donor
or acceptor sites.  This helps to determine if this junction might be involved
in alternative splicing.

Metric 17 - Primary Junction
----------------------------

If this is not a unique junction (see Metric 16), then this is a primary junction
if it has the most spliced reads when compared to the other junctions sharing its
donor or acceptor sites.  If this is a unique junction, then it is also a primary
junction.

Metric 18 - Multiple Mapping Score
----------------------------------

The multiple mapping score is the number of spliced reads associated with the junction
divided by the number of times those same reads are found mapped anywhere in the genome.
Therefore a score of 1 indicates that all spliced reads associated with the junction
are only found in this junction.  A low score would indicate that the those reads map
to multiple locations across the genome.

Originally described in TrueSight paper.


Metric 19 - Mean mismatches
--------------------------------

This is the mean number of mismatches found across all spliced reads supporting the
junction.  This includes any mismatches at any point along the spliced read, which
includes mismatches even if they are the otherside of another junction in the case 
of an MSR.

Originally described in TrueSight paper.


Metric 20 - Number of Multiple Spliced Reads
--------------------------------------------

This is a count of the number of spliced reads that support the junction that also
support another junction.


Metric 21 / 22 - Number of Upstream and Downstream Junctions
------------------------------------------------------------

The number of upstream and downstream junctions contained within any MSRs associated
with this junction.  Will be 0 for junctions without any MSRs.


Metric 23 / 24 - Number of Upstream and Downstream Alignments
---------------------------------------------------------------

This is a count of the number of unspliced reads aligning upstream of the splice 
junction, that overlap with the upstream anchor.  Caution must be taken interpreting
this metric closely packed introns could mean the presence of MSRs exclude the possibility
of getting any unspliced upstream alignments.  In addition, if the junction is close
to the sequence start, it maybe that no unspliced upstream alignments are possible
either.