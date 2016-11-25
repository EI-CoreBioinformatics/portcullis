.. _requirements:


System requirements
===================

Portcullis supports Unix, linux or Mac systems.  Windows may work but hasn't been
tested.  Portcullis can take advantage of multiple cores via multi-threading where possible. 

A minimum of 4GB RAM, which will enable you to process small datasets, but
we suggest at least 16GB to handle most medium sized datasets.  
Large datasets will require more RAM, the largest dataset we tested peaked at 50GB 
using 16 threads.  Generally, except for extremely deep datasets, if you have sufficient 
memory to align the reads, you should have sufficient memory to run portcullis.  

In terms of resource profiling the portcullis pipeline, the preparation, junction 
filtering and BAM filtering :ref:`subtools <using>` have relatively modest
memory requirements.  The amount of memory the portcullis junction analysis stage requires 
will depend on the quantity of RNAseq reads, the number of potential junctions, 
and the contiguity of your genome and the number of threads used.  We have optimised 
the code so that we only keep reads associated with junctions in memory 
as long as necessary to accurately capture all the information about that junction.  
As soon as we move beyond scope of a junction by encountering a read that exists beyond
a junctions boundaries, we calculate the junctions metrics and discard all reads associated
to that junction.  Therefore a single very highly expressed junction can increase peak
memory usage.  In addition, because, threading processes multiple target sequences in parallel, the higher number of threads
the higher the memory requirements. So try reducing the number of threads if you encounter
out of memory issues or segmentation faults.
