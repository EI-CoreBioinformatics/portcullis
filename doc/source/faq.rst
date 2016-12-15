.. _faq:


Frequently Asked Questions
==========================

Can I reduce portcullis' memory usage?
--------------------------------------

Portcullis is optimised to reduce memory usage where possible, so in most cases
memory usage should not be an issue if you were able to align your reads in the first
place.  However, if you are encountering out of memory issues, there are a few factors 
that can influence memory usage in portcullis, particularly during the junction
analysis stage.  The first is the number of threads used, more threads require more
memory.  You can therefore reduce memory usage at the expense of runtime.  Second,
we have an ability to process additional metrics in the junction analysis stage
called ``--extra``, this can require large amounts of memory, so unless you really
need this option on (there should be no direct impact on the default filtering
setup) the switch this off.  Finally, the main driver for memory usage is the depth
of the dataset, specifically highly covered junctions.  Should have have tens of 
thousands of reads supporting a single junction, all these reads must be kept in
memory.  You can therefore reduce memory usage by pre-processing the BAM files to
either cap the depth (you can use the `Kmer Analysis Toolkit <https://github.com/TGAC/KAT>`_ 
to identify high coverage kmers then remove reads associated with those kmers), 
or downsample using ``samtools view -s``. 

