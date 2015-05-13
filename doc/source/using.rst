.. _using:

Using Portcullis
================

Portcullis is a C++ program containing a number of subtools which can be used in
isolation or as part of a pipeline.  Typing ```portcullis --help``` will show a
list of the available subtools.  Each subtool has its own help system which you 
can access by typing ```portcullis <subtool> --help```.  

Quickstart
==========

For the impatient, just use ```portcullis full -t <threads> <genome> (<bam>)+``` 
to get up and running.  This runs the following subtools in order:

* Prepare
* Junction Analysis
* Filtering
* Bam Filtering


Subtools
========

Prepare
-------

This prepares all the input data into a format suitable for junction analysis.  Specifically,
this merges the input BAMs if more than one was required.  It then ensures the BAM
is both sorted and indexed.  It also indexes the genome.


Junction Analysis
-----------------

Collects all potential junctions identifiable within the input BAM.  The analyses
each one to collect a number of measurements (see _metrics for more details).


Junction Filtering
------------------

We provide several options for filtering out false positive junctions.  The first
is a rule-based approach which allows the user to decide what metric values constitute
reliable junctions.  The user must provide an JSON file describing their filter profile.
Additionally, we have some existing filter profiles, for the user, which can be used
for execution or just as a template for deriving their own filter profile.

We also provide a trained one-class SVM model which can make it's own decisions about
whether a given junction is genuine or not.  This normally provides good results
and gives the user a confidence score for each junction, however it does not provide
an easy to interpret justification for decisions made so the user may prefer the 
manual filtering method previously described.  

<JUSTIFY SVM HERE>  <DESCRIBE TRAINING SET>


Bam Filtering
-------------




Memory Usage
============

The amount of memory portcullis requires will depend on the quantity of RNAseq reads,
the number of potential junctions, and the contiguity of your genome.  Currently (V0.4.1),
each potential junction within a single reference sequence is kept in memory at the
same time.  
