.. _using:

Using Portcullis
================

Portcullis is a C++11 program containing a number of subtools which can be used in
isolation or as part of a pipeline.  Typing ``portcullis --help`` will show a
list of the available subtools.  Each subtool has its own help system which you 
can access by typing ``portcullis <subtool> --help``.  

The list of subtools present in portcullis are listed below in order.  A more detailed
description of each subtool follows in the subsequent sections:

* :ref:`prepare`
* :ref:`junc`
* :ref:`filt`
* :ref:`bamfilt` (optional)

However, it is possible to run all portcullis steps in one go by using the `full`
subtool.  The command line usage for this option is as follows::

    Usage: portcullis full [options] <genome-file> (<bam-file>)+

    System options:
      -t [ --threads ] arg (=1) The number of threads to use.  Note that increasing the number of threads will also increase memory requirements.  Default: 1
      -v [ --verbose ]          Print extra information
      --help                    Produce help message

    Output options:
      -o [ --output ] arg (="portcullis_out") Output directory. Default: portcullis_out
      -b [ --bam_filter ]                     Filter out alignments corresponding with false junctions.  Warning: this is time consuming; make sure you really want 
                                              to do this first!
      --exon_gff                              Output exon-based junctions in GFF format.
      --intron_gff                            Output intron-based junctions in GFF format.
      --source arg (=portcullis)              The value to enter into the "source" field in GFF files.

    Input options:
      --force               Whether or not to clean the output directory before processing, thereby forcing full preparation of the genome and bam files.  By 
                            default portcullis will only do what it thinks it needs to.
      --copy                Whether to copy files from input data to prepared data where possible, otherwise will use symlinks.  Will require more time and disk 
                            space to prepare input but is potentially more robust.
      --use_csi             Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it supports very long target sequences (probably not 
                            an issue unless you are working on huge genomes).  BAI has the advantage that it is more widely supported (useful for viewing in genome 
                            browsers).

    Analysis options:
      --orientation arg (=UNKNOWN)  The orientation of the reads that produced the BAM alignments: "F" (Single-end forward orientation); "R" (single-end reverse 
                                    orientation); "FR" (paired-end, with reads sequenced towards center of fragment -> <-.  This is usual setting for most Illumina 
                                    paired end sequencing); "RF" (paired-end, reads sequenced away from center of fragment <- ->); "FF" (paired-end, reads both 
                                    sequenced in forward orientation); "RR" (paired-end, reads both sequenced in reverse orientation); "UNKNOWN" (default, 
                                    portcullis will workaround any calculations requiring orientation information)
      --strandedness arg (=UNKNOWN) Whether BAM alignments were generated using a type of strand specific RNAseq library: "unstranded" (Standard Illumina); 
                                    "firststrand" (dUTP, NSR, NNSR); "secondstrand" (Ligation, Standard SOLiD, flux sim reads); "UNKNOWN" (default, portcullis will 
                                    workaround any calculations requiring strandedness information)
      --separate                    Separate spliced from unspliced reads.
      --extra                       Calculate additional metrics that take some time to generate.  Automatically activates BAM splitting mode (--separate).

    Filtering options:
      -r [ --reference ] arg Reference annotation of junctions in BED format.  Any junctions found by the junction analysis tool will be preserved if found in this 
                             reference file regardless of any other filtering criteria.  If you need to convert a reference annotation from GTF or GFF to BED format
                             portcullis contains scripts for this.
      --max_length arg (=0)  Filter junctions longer than this value.  Default (0) is to not filter based on length.
      --canonical arg (=OFF) Keep junctions based on their splice site status.  Valid options: OFF,C,S,N. Where C = Canonical junctions (GT-AG), S = Semi-canonical 
                             junctions (AT-AC, or  GC-AG), N = Non-canonical.  OFF means, keep all junctions (i.e. don't filter by canonical status).  User can 
                             separate options by a comma to keep two categories.
      --min_cov arg (=1)     Only keep junctions with a number of split reads greater than or equal to this number
      --save_bad             Saves bad junctions (i.e. junctions that fail the filter), as well as good junctions (those that pass)


This is the typical way to run portcullis but it's still helpful to know what each step
in the pipeline does in more detail.  Also the subtools offer some additional controls 
that can be useful in certain situations so please read on.


.. _prepare:

Prepare
-------

This prepares all the input data into a format suitable for junction analysis.  Specifically,
this merges the input BAMs if more than one was provided.  Using samtools, it then 
ensures the BAM is both sorted and both the sorted BAM and genome are indexed.
The prepare output directory contains all inputs in a state suitable for 
downstream processing by portcullis.

Normally we try to minimise the work done by avoiding re-sorting or indexing if 
the files are already in a suitable state.  However, options are provided should
the user wish to force re-sorting and re-indexing of the input.

Usage
~~~~~
::

    Usage: portcullis prep [options] <genome-file> (<bam-file>)+ 

    Options:
      -o [ --output ] arg (="portcullis_prep") Output directory for prepared files.
      --force                                  Whether or not to clean the output directory before processing, thereby forcing full 
                                               preparation of the genome and bam files.  By default portcullis will only do what it thinks 
                                               it needs to.
      --copy                                   Whether to copy files from input data to prepared data where possible, otherwise will use 
                                               symlinks.  Will require more time and disk space to prepare input but is potentially more 
                                               robust.
      -c [ --use_csi ]                         Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it 
                                               supports very long target sequences (probably not an issue unless you are working on huge 
                                               genomes).  BAI has the advantage that it is more widely supported (useful for viewing in 
                                               genome browsers).
      -t [ --threads ] arg (=1)                The number of threads to used to sort the BAM file (if required).  Default: 1
      -v [ --verbose ]                         Print extra information
      --help                                   Produce help message


.. _junc:

Junction Analysis
-----------------

Portcullis is designed to be as portable as possible so it does not rely on esoteric 
SAM tags and other artifacts that are not consistently present in all SAM/BAMs.  
In this stage, Portcullis analyses the BAM file to look for alignments containing 
gaps (REFSKIP 'N' cigar ops) and creates a detailed analysis of all distinct gaps 
detected, these are considered as potential junctions.  A number of observations 
are made for each junction such as number of supporting split reads, how those
reads are distributed around the junction, the actual nucleotides representing 
the splice sites and how repetitive the genomic region is around the splice sits 
(see _metrics for more details of all measurements taken).  Portcullis outputs
all junctions found in the BAM in a number of formats such as GFF3, BED, although
the most detailed information can be found in the tab file.

Usage
~~~~~
::

    Usage: portcullis junc [options] <prep_data_dir>

    System options:
      -t [ --threads ] arg (=1)     The number of threads to use.  Note that increasing the number of threads will also 
                                    increase memory requirements.
      -s [ --separate ]             Separate spliced from unspliced reads.
      --extra                       Calculate additional metrics that take some time to generate.  Automatically activates BAM
                                    splitting mode (--separate).
      --orientation arg (=UNKNOWN)  The orientation of the reads that produced the BAM alignments: "F" (Single-end forward 
                                    orientation); "R" (single-end reverse orientation); "FR" (paired-end, with reads sequenced
                                    towards center of fragment -> <-.  This is usual setting for most Illumina paired end 
                                    sequencing); "RF" (paired-end, reads sequenced away from center of fragment <- ->); "FF" 
                                    (paired-end, reads both sequenced in forward orientation); "RR" (paired-end, reads both 
                                    sequenced in reverse orientation); "UNKNOWN" (default, portcullis will workaround any 
                                    calculations requiring orientation information)
      --strandedness arg (=UNKNOWN) Whether BAM alignments were generated using a type of strand specific RNAseq library: 
                                    "unstranded" (Standard Illumina); "firststrand" (dUTP, NSR, NNSR); "secondstrand" 
                                    (Ligation, Standard SOLiD, flux sim reads); "UNKNOWN" (default, portcullis will workaround
                                    any calculations requiring strandedness information)
      -c [ --use_csi ]              Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage that it 
                                    supports very long target sequences (probably not an issue unless you are working on huge 
                                    genomes).  BAI has the advantage that it is more widely supported (useful for viewing in 
                                    genome browsers).
      -v [ --verbose ]              Print extra information
      --help                        Produce help message

    Output options:
      -o [ --output ] arg (=portcullis_junc/portcullis)
                                             Output prefix for files generated by this program.
      --exon_gff                             Output exon-based junctions in GFF format.
      --intron_gff                           Output intron-based junctions in GFF format.
      --source arg (=portcullis)             The value to enter into the "source" field in GFF files.


.. _filt:

Junction Filtering
------------------

Portcullis provides various means for filtering junctions detected in the input 
BAM file.  By default we use a machine learning approach, which trains on a high-
confidence subset of the data, and then applies the trained model to the full set
in order to score each junction.  By default scores of over 0.5 are classed as
genuine and those under as invalid.  The user can control the threshold value via
a command line option.  We use the a modified version of the ranger random forest 
code as the learner.  Portcullis outputs all junctions passing the filter in a number 
of formats such as GFF3, BED, and TSV format.  The user can also request all
junctions failing the filter are output into an additional set of GFF, BED and TSV files.


Reference annotations
~~~~~~~~~~~~~~~~~~~~~

Should the user have access to a reference annotation, they can supply that ( via the `--reference` command line option) so
that should portcullis filter out any junctions that are also found in the reference,
then those are put back into the set of genuine junctions.  This feature is useful
when working with model organisms where high-quality references are available.

The portcullis filter tool requires the reference junction annotation in BED format.
If this is not readily available Portcullis comes supplied with an addition toolkit 
called :ref:`junctools`, which can convert GTF annotation files to a set of junctions in BED format.


Validating results
~~~~~~~~~~~~~~~~~~

Should the user know whether each junction in the input set is genuine or not, that
can be provided to portcullis via the `--genuine` command line option.  This file
takes the format of a line separated list of either `1` indicating genuine and `0`
indicating invalid in the same order as the input junctions.  Portcullis, can then
measure the performance of its filtering strategy.


Rule-based filtering
~~~~~~~~~~~~~~~~~~~~

Alternatively, the user can filter junctions based on simple rules applied to the junction
metrics.  They do this via a JSON file describing their filter profile, which is
passed to the filter tool via the `--filter_file` command line option.  Examples
are provided in the `data` sub-directory, which can be used directly, or as a template 
for deriving a custom filter profile.  The rules can be combined 
using logic operations (and / or / not, etc) and applied to the full set of input 
junctions.

Here's an example set of rules that must all be satisfied to pass this filter::

    {
            "parameters": {
                    "M4-nb_rel_aln": {
                            "operator": "gte",
                            "value": 2
                    },
                    "M12-maxmmes": {
                            "operator": "gte",
                            "value": 10
                    },
                    "M11-entropy": {
                            "operator": "gte",
                            "value": 1.5
                    },
                    "M13-hamming5p": {
                            "operator": "gte",
                            "value": 2
                    },
                    "M14-hamming3p": {
                            "operator": "gte",
                            "value": 2
                    }
            },
            "expression": "M4-nb_rel_aln & M11-entropy & M12-maxmmes & M13-hamming5p & M14-hamming3p"  
    }


Filtering with a pre-made model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although it is generally not recommended, the user can re-use existing random forest
models to apply to new datasets.  This is done via the `--model_file` option.


Usage
~~~~~
::

    Usage: portcullis filter [options] <prep_data_dir> <junction_tab_file>

    System options:
      -t [ --threads ] arg (=1) The number of threads to use during testing (only applies if using forest model).
      -v [ --verbose ]          Print extra information
      --help                    Produce help message

    Output options:
      -o [ --output ] arg (="portcullis_filter/portcullis")
                                             Output prefix for files generated by this program.
      -b [ --save_bad ]                      Saves bad junctions (i.e. junctions that fail the filter), as well as good 
                                             junctions (those that pass)
      --exon_gff                             Output exon-based junctions in GFF format.
      --intron_gff                           Output intron-based junctions in GFF format.
      --source arg (=portcullis)             The value to enter into the "source" field in GFF files.

    Filtering options:
      -f [ --filter_file ] arg If you wish to custom rule-based filter the junctions file, use this option to provide a list 
                               of the rules you wish to use.  By default we don't filter using a rule-based method, we instead
                               filter via a self-trained random forest model.  See manual for more details.
      -r [ --reference ] arg   Reference annotation of junctions in BED format.  Any junctions found by the junction analysis 
                               tool will be preserved if found in this reference file regardless of any other filtering 
                               criteria.  If you need to convert a reference annotation from GTF or GFF to BED format 
                               portcullis contains scripts for this.
      -n [ --no_ml ]           Disables machine learning filtering
      --max_length arg (=0)    Filter junctions longer than this value.  Default (0) is to not filter based on length.
      --canonical arg (=OFF)   Keep junctions based on their splice site status.  Valid options: OFF,C,S,N. Where C = 
                               Canonical junctions (GT-AG), S = Semi-canonical junctions (AT-AC, or GC-AG), N = Non-canonical.
                                 OFF means, keep all junctions (i.e. don't filter by canonical status).  User can separate 
                               options by a comma to keep two categories.
      --min_cov arg (=1)       Only keep junctions with a number of split reads greater than or equal to this number
      --threshold arg (=0.5)   The threshold score at which we determine a junction to be genuine or not.  Increase value 
                               towards 1.0 to increase precision, decrease towards 0.0 to increase sensitivity.  We generally 
                               find that increasing sensitivity helps when using high coverage data, or when the aligner has 
                               already performed some form of junction filtering.

.. _bamfilt:

Bam Filtering
-------------

Portcullis can also filter the original BAM file removing alignments 
associated with `bad` junctions.  Both the filtered junctions and BAM files are cleaner
and more usable resources which can more effectively be used to assist in downstream 
analyses such as gene prediction and genome annotation. 

Usage
~~~~~
::

    Usage: portcullis bamfilt [options] <junction-file> <bam-file>

    Options:
      -o [ --output ] arg (="filtered.bam")   Output BAM file generated by this program.
      -s [ --strand_specific ] arg (=UNKNOWN) Whether BAM alignments were generated using a strand specific RNAseq library: 
                                              "unstranded" (Standard Illumina); "firststrand" (dUTP, NSR, NNSR); 
                                              "secondstrand" (Ligation, Standard SOLiD, flux sim reads)  Default: 
                                              "unstranded".  By default we assume the user does not know the strand specific 
                                              protocol used for this BAM file.  This has the affect that strand information is
                                              derived from splice site information alone, assuming junctions are either 
                                              canonical or semi-canonical in form.  Default: "unknown"
      -c [ --clip_mode ] arg (=HARD)          How to clip reads associated with bad junctions: "HARD" (Hard clip reads at 
                                              junction boundary - suitable for cufflinks); "SOFT" (Soft clip reads at junction
                                              boundaries); "COMPLETE" (Remove reads associated exclusively with bad junctions,
                                              MSRs covering both good and bad junctions are kept)  Default: "HARD"
      -m [ --save_msrs ]                      Whether or not to output modified MSRs to a separate file.  If true will output 
                                              to a file with name specified by output with ".msr.bam" extension
      -c [ --use_csi ]                        Whether to use CSI indexing rather than BAI indexing.  CSI has the advantage 
                                              that it supports very long target sequences (probably not an issue unless you 
                                              are working on huge genomes).  BAI has the advantage that it is more widely 
                                              supported (useful for viewing in genome browsers).
      -v [ --verbose ]                        Print extra information
      --help                                  Produce help message

