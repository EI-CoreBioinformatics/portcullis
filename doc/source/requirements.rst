.. _requirements:


System requirements
===================

Portcullis supports Unix, linux or Mac systems.  Windows may work but hasn't been
tested.  Portcullis can take advantage of multiple cores via multi-threading where possible. 

A minimum of 4GB RAM, which will enable you to process small datasets, but
we suggest at least 16GB to handle most medium sized datasets.  
Large datasets will require more RAM, the largest dataset we tested required 20GB 
but the actual amount of memory required depends on how many splice junctions
are present in your BAM file and the number of accompanying reads that support them.

More specifically, the preparation, junction filtering and BAM filtering :ref:`subtools <using>` have relatively modest
memory requirements.  The amount of memory the portcullis junction analysis stage requires 
will depend on the quantity of RNAseq reads, the number of potential junctions, 
and the contiguity of your genome.  We have optimised the code so that we only keep reads associated with junctions in memory 
as long as necessary to accurately capture all the information about that junction.  
As soon as we move beyond scope of a junction by encountering a read that exists beyond
a junctions boundaries, we calculate the junctions metrics and discard all reads associated
to that junction.  Because, threading processes multiple target sequences in parallel, the higher number of threads
the higher the memory requirements. In addition, we also require a copy of the genome 
index in memory, although this copy is shared amongst threads.  This information
can help to estimate memory usage prior to running portcullis.
