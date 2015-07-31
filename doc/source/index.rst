.. portcullis documentation master file, created by
   sphinx-quickstart on Tue Dec  2 15:59:49 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to portcullis' documentation!
======================================

#.. image:: images/spectre.png
#    :scale: 50%

Portcullis is designed to predict junctions from aligned RNA-seq data.  We expect
the user to have already generated a BAM file using a splice aware aligner of their
choice.  For example, Tophat, Gsnap, STAR or HISAT will work fine.  Portcullis is designed to
be as portable as possible so where possible does not rely on esoteric SAM tags and other
artifacts that are not consistently present in all SAM/BAMs.  Portcullis
will then analyse the BAM file to look for alignments containing gaps (REFSKIP 'N' cigar ops)
and create a detailed analysis of all distinct gaps found in the BAM file, these
are considered as potential junctions.  Portcullis provides various means (both 
manual and/or automatic) of filtering these potential junctions in order to remove 
false positives.  Portcullis can also filter the original BAM file removing alignments 
associated with `bad` junctions.  Both the filtered junctions and BAM files are cleaner
and more usable resources which can more effectively be used to assist in downstream 
analyses such as gene prediction and genome annotation.  


Contents:

.. toctree::
    :numbered:
    :maxdepth: 2

    installation
    using
    metrics
    

.. _system:
System requirements
===================

Portcullis supports Unix, linux or Mac systems.  Windows may work but hasn't been
tested.  A minimum of 8GB RAM, which will enable you to process small - medium sized datasets.  
Large datasets will require more RAM (potentially a lot more), the actual amount of
memory required depends on how many spliced alignments are present in your BAM file.
Portcullis does not need to store the whole spliced alignment, just a hashcode from the name.


.. _issues:
Issues
======

Should you discover any issues with spectre, or wish to request a new feature please raise a ticket at https://github.com/maplesond/portcullis/issues.
Alternatively, contact Daniel Mapleson at: daniel.mapleson@tgac.ac.uk


.. _availability:
Availability and License
========================

Open source code available on github: https://github.com/maplesond/portcullis.git

Spectre is available under GNU GLP V3: http://www.gnu.org/licenses/gpl.txt


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

