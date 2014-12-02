Metrics
=======

During the junction making stage, Portcullis generates a list of potential splice
junctions by assuming any RNAseq alignment containing an 'N' cigar operation.  For
all potential junctions portcullis makes a number of observations about how the
reads align around the region.  The complete list of observations (or metrics) are
described in this section.  Before we do that, here is a small glossary of terms
that we will use throughout::

* Splice junction - Splice junctions are points on a DNA sequence at which `superfluous' DNA is removed during the process of protein creation in higher organisms.  For the purposes of this tool a splice junction is essentially the same as an intron.
* Splice site - A junction has both a 5' donor and 3' acceptor site, which mark the start and end of the junction.  Both donor and acceptor sites are 2bp long, and usually contain a canonical motif `GT*AG`, or its reverse complement on the negative strand.
* Intron - Any nucleotide sequence within a gene that is removed by RNA splicing while the final mature RNA product of a gene is being generated.  For the purposes of this tool an intron is essentially the same as a splice junction.
* Anchor - The part of a spliced read that aligns to the genome.  This will represent regions both upstream and downstream of the splice junction.




Metric 1 - Number of spliced reads
----------------------------------

This a count of the number of reads containing an 'N' cigar operation at exactly
this genomic position represented by this junction.  Because we use 'N' cigar operations
to derive all potential junctions each junction with have a value of 1 or more for
this metric.

#TODO pretty pic


Metric 2 - Canonical splice site
--------------------------------

Most splice junctions contain distinctive dinucleotide pattern at the start and 
end of the splice junction.  These patterns are called `donors` and `acceptors`
respectively.  If the potential junction has the typical `GT_AG` motif (or its
reverse complement on the negative strand then it is considered a canonical splice
site.  Alternatively, if the junction contains either `AT_AC` or `GC_AG` or their
reverse complements then it is considered a semi-canonical splice site.  Other motifs
at the donor and acceptor sites are considered novel splice sites.

#TODO pretty pic


Metric 3 - Intron size
----------------------

This is the length of the region defined by the 'N' cigar operation.  Which should
also be the size of the intron if this is a genuine splice site.


Metric 4 - Max-min-Anchor
-------------------------

This is the maximum of the minimum anchor size of all spliced reads associated with
this junction.  

#TODO pretty pic

Metric 5 - Difference between Anchors
-------------------------------------

This is the difference between the minimum L/R anchor and the maximum L/R anchor
with the smaller of the two values reported.

#TODO pretty pic


Metric 6 - Entropy
------------------

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

