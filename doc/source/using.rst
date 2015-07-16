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
is both sorted and indexed.  It also indexes the genome and writes any settings to used to disk.
The directory can then be used downstream for junction analysis, safe in the knowledge
that the input is in a state suitable for processing.

If requested we can either try to minimise work here where possible by avoiding 
resorting or indexing if the files are already in a suitable state, or we allow
the user to force resorting and reindexing.


Junction Analysis
-----------------

Collects all potential junctions identifiable within the input BAM.  The analyses
each one to collect a number of measurements (see _metrics for more details).

An important option with respect to runtime calculations is the "-f" option.  This
skips some of the less important metrics and saves a lot of runtime.  For the time
being (V0.8.2) we recommend the user has this option selected.  In the future
we might modify the logic here to make the user request additional metrics if they
want them.


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
the number of potential junctions, and the contiguity of your genome.  Currently (V0.8.2),
we have optimised the code so that we only keep reads associated with junctions in memory 
as long as necessary to accurately capture all the information about that junction.
As soon as we move beyond scope of a junction by encountering a read that exists beyond
a junctions boundaries, we calculate the junctions metrics and discard all reads associated
to that junction.
  
