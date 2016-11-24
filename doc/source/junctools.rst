.. _junctools:

Junctools
=========

Portcullis contains a tools suite called junctools for interpretting, analysing 
and converting junction files.  A description of all tools within
junctools can be displayed by typing ``junctools --help`` at the command line:

::

    usage: This script contains a number of tools for manipulating junction files.
           [-h] {compare,convert,markup,set,split} ...

    optional arguments:
      -h, --help            show this help message and exit

    Junction tools:
      {compare,convert,markup,set,split}
        compare             Compares junction files.
        convert             Converts junction files between various formats.
        markup              Marks whether each junction in the input can be found in the reference or not.
        set                 Apply set operations to two or more junction files.
        split               Splits portcullis pass and fail juncs into 4 sets (TP, TN, FP, FN) based on whether or not the junctions are found in the reference or not.

If this help message does not appear then please follow the instructions in the
next section.

.. note:: Junctools is designed to support portcullis tab delimited files as well as
    many commonly used flavours of the bed format.  Junctools will generally auto-detect
    the format based on the file extension and interpret the file as necessary.



Installation
------------

If you have proceeded through the installation steps and have run ``make install`` 
then junctools will be installed along with portcullis as a python package to the
location specified by ``--prefix`` when running the ``configure`` script.  Should
prefix not be specified then junctools will be installed to ``/usr/local``.  The 
junctools library reside in the ``$(prefix)/lib/python3.x/site-packages`` (where x 
is the minor version of python installed on your system) directory and
the junctools executable will reside in ``$(prefix)/bin``.  If you wish to install
junctools into a specific python environment or if you wish to install junctools
without portcullis, first make sure that you chosen environment has precedence on 
your system then then ``cd <portcullis_dir>/scripts``, and then type ``python3 setup.py install``.

.. note:: When you type ``make uninstall`` for portcullis you will also uninstall
    junctools from the location specified by $(prefix).  Keep this in mind if you 
    installed to a system directory you will if you prefix ``sudo`` to the uninstall command.
    Alternatively, if you installed junctools manually to a specific python environment then you
    can uninstall junctools with: ``pip3 uninstall -y junctools``.
    

Compare
-------

This tool compares one or more junction files against a reference file.  

The usage information for this tool is as follows:

::

    usage: Compares junction files.
           [-h] [-s] [-l LABELS] [-m] reference input [input ...]

    positional arguments:
      reference             The junction file to treat as the reference
      input                 One or more junction files to compare against the
                            reference

    optional arguments:
      -h, --help            show this help message and exit
      -s, --use_strand      Whether to use strand information when building keys
      -l LABELS, --labels LABELS
                            Path to a file containing labels for the reference
                            indicating whether or not each reference junction is
                            genuine (as generated using the markup tool). If
                            provided this script produces a much richer
                            performance analysis. Not compatible with '--
                            multiclass'
      -m, --multiclass      Breakdown results into multiple classes: 1) Matching
                            intron 2) Two matching splice sites but no matching
                            intron (i.e. splice sites from different introns) 3)
                            One matching splice site 4) No matching splice sites

Essentially this tool can be run in one of three modes depending on the options
selected.  The first mode does a information retrieval style comparison of the input
files against a reference, as it's not possible to identify true negatives.  The
output contains counts of junctions, plus true positives, false positives and false
negatives, along with precision, recall and F1 scores.  The example below shows
output from comparing 6 bed files against a reference bed using the following command
``junctools compare ../ref.bed *.bed``:

::

    Reference:
    - # total junctions: 158156
    - # distinct junctions: 158156

    File     distinct    total     TP        FP      FN      REC    PRC    F1
    f1.bed   146800      146800    135570    11230   22586   85.72  92.35  88.91
    f2.bed   146778      146800    0         146778  158156  0.00   0.00   0.00
    f3.bed   130905      130905    129103    1802    29053   81.63  98.62  89.33
    f4.bed   130905      130905    129103    1802    29053   81.63  98.62  89.33
    f5.bed   118865      118865    117569    1296    40587   74.34  98.91  84.88
    f6.bed   117613      117613    113952    3661    44204   72.05  96.89  82.64

    Mean recall:  65.89
    Mean precision:  80.90
    Mean f1:  72.51

.. note:: Please keep in mind that if comparing predicted junctions from a real RNAseq
    sample against a reference that the sample may contain genuinely novel junctions
    that will appear as false positives.  Generally we use this method on simulated
    datasets where we know the definitive set of true junctions.



The second mode provides a more detailed analysis which is useful for comparing
portcullis results against marked up junctions from an alignment tool.  In this
case we markup junctions from an alignment tool using the :ref:`markup tool <markup>`
.  This is essentially a list specifying whether or not each junction found by the
aligner is present in a reference or not.  We can then compare results from portcullis
against the marked up alignment junctions.  This gives us a definitive set of false
negative junctions, i.e. junctions from the aligner that were genuine but incorrectly
marked as negative by portcullis.

Finally, the third mode is useful for comparing junctions from real RNAseq datasets
against a real reference.  This breaks down results into the 4 classes: 

 * 1 - Matching intron
 * 2 - Two matching splice sites but no matching intron (i.e. splice sites from different introns) 
 * 3 - One matching splice site 
 * 4 - No matching splice sites

This approach allows the user to better understand the set of false positives produced
from real datasets, and can give some indication of whether a junction is a novel
junction or a false positive.


.. note:: By default we do not use strand information when determining the location
    of a junction.  To clarify, it is possible that bed file contains multiple junctions with
    the same sequence, start and stop sites but with a different strand.  By default
    ``junctools compare`` will collapse these as a duplicate junction.  Although
    not immediately intuitive this allows us to circumvent problems from junctions
    that have unknown strand.  This is important as some tools do not output strand
    information.  However, should you wish to disable this feature you can do so
    with the ``--use_strand`` option.

Convert
-------

This can convert junction files between various commonly used formats.  The conversion
tool supports the following commonly used formats:

 * bed        = (Input only) BED format - we automatically determine if this is BED 6 or 12 format, as well as if it is intron, exon or tophat style).
 * ebed       = (Output only) Portcullis style exon-based BED12 format (Thick-start and end represent splice sites).
 * tbed       = (Output only) Tophat style exon-based BED12 format (splice sites derived from blocks).
 * ibed       = (Output only) Intron-based BED12 format.
 * bed6       = (Output only) BED6 format (BED6 files are intron-based).
 * gtf        = (Input only) Transcript assembly or gene model containing transcript and exon features.  NOTE: output will only contain junctions derived from this GTF.
 * gff        = (Input only) Transcript assembly or gene model containing introns to extract. NOTE: input must contain "intron" features, and output will only contain these introns represented as junctions.
 * egff       = (Output only) Exon-based junctions in GFF3 format, uses partial matches to indicate exon anchors.
 * igff       = (Output only) Intron-based junctions in GFF3 format

In addition we support the following application specific tab delimited formats:

 * portcullis = Portcullis style tab delimited output.
 * hisat      = HISAT style tab delimited format.
 * star       = STAR style tab delimited format.
 * finesplice = Finesplice style tab delimited format.
 * soapslice  = Soapsplice style tab delimited format.
 * spanki     = SPANKI style tab delimited format.
 * truesight  = Truesight style tab delimited format.

The usage information for the conversion tool looks like this::


    usage: Converts junction files between various formats.
           [-h] -if INPUT_FORMAT -of OUTPUT_FORMAT [-o OUTPUT] [-is] [-d] [-s]
           [-r] [--index_start INDEX_START] [--prefix PREFIX] [--source SOURCE]
           input

    positional arguments:
      input                 The input file to convert

    optional arguments:
      -h, --help            show this help message and exit
      -if INPUT_FORMAT, --input_format INPUT_FORMAT
                            The format of the input file to convert.
      -of OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                            The output format.
      -o OUTPUT, --output OUTPUT
                            Output to this file.  By default we print to stdout.
      -is, --ignore_strand  Whether or not to ignore strand when creating a key for the junction
      -d, --dedup           Whether or not to remove duplicate junctions
      -s, --sort            Whether or not to sort the junctions.  Note that sorting requires all junctions to be loaded into memory first.  This maybe an issue for very large input files.
      -r, --reindex         Whether or not to reindex the output.  The index is applied after prefix.
      --index_start INDEX_START
                            The starting index to apply if the user requested reindexing
      --prefix PREFIX       The prefix to apply to junction ids if the user requested reindexing
      --source SOURCE       Only relevant if output is GFF format, use this option to set the source column in the GFF


.. note:: The user can also use the conversion tool to deduplicate, sort and reindex junction files.

.. _markup:

Markup
------

This tool marks whether each junction in the input can be found in the reference 
or not.  Output from the tool is a line seperated list of 1's (indicating junction
is found in the reference) and 0's (indicating the junction is not found in the
reference).  Output is written to a file with the same name as the input except a 
'.res' extension is added.  Usage information follows::

    usage: Marks whether each junction in the input can be found in the reference or not.
           [-h] [-o OUTPUT_DIR] [-s] reference input [input ...]

    positional arguments:
      reference             The junction file to treat as the reference
      input                 One or more junction files to compare against the
                            reference

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                            If output dir is specified this will create output
                            files for each input file with a .res extension
                            indicating whether or not the junction was found in
                            the reference. By default we write out a .res file in
                            the same directory as the input file was found in.
      -s, --use_strand      Whether to use strand information when building keys


Set
---

Apply set operations to two or more junction files.  This tool supports various
different ways to apply set operations between junction files.  First you can merge
two or more junction files using the following modes:

 * intersection = Produces the intersection of junctions from multiple input files
 * union        = Produces the union of junctions from multiple input files
 * consensus    = If there are 3 or more input files, the consensus operation produces a merged set of junctions where those junctions are found across a user-defined number of input files

Output from these modes potentially involves mergeing multiple junctions from various
files into a single representative.  When this occurs junction anchors are extended
representing the most extreme extents found across all junctions at the given site.
In addition, the junction score is modified according to the setting selected by
the user, by default this involves summing the scores of all junctions, although
the user can alternatively choose to take the min, max or mean of the values.

The following modes only support two input files and produce an output file containing
junctions which are taken directly from the input without modification:

 * subtract     = Produces an output set of junctions containing all junctions present
                in the first input file that also are not found in the second file
 * filter       = Produces an output set of junctions containing all junctions present
                in the first input file that are also found in the second file.  This
                is similar to an intersection on two files except that duplicates and
                additional content assigned to junctions in the first file are retained
 * symmetric_difference =
                Produces an output set containing junctions from both input files
                that are not present in the intersection of both

In addition, these test modes also only support 2 input files and return either 
True or False depending on the test requested:

 * is_subset    = Returns True if all junctions in the first file are present in the second
 * is_superset  = Returns True if all junctions in the second file are present in the first
 * is_disjoint  = Returns True if there is a null intersection between both files

Usage::

    usage: Apply set operations to two or more junction files.
           [-h] [-m MIN_ENTRY] [--operator OPERATOR] [-o OUTPUT] [-p PREFIX] [-is]
           mode input [input ...]

    positional arguments:
      mode                  Set operation to apply.  Available options:
                             - intersection
                             - union
                             - consensus
                             - subtract
                             - symmetric_difference
                             - is_subset
                             - is_superset
                             - is_disjoint
      input                 List of junction files to merge (must all be the same file format)

    optional arguments:
      -h, --help            show this help message and exit
      -m MIN_ENTRY, --min_entry MIN_ENTRY
                            Minimum number of files the entry is require to be in.  0 means entry must be
                            present in all files, i.e. true intersection.  1 means a union of all input files
      --operator OPERATOR   Operator to use for calculating the score in the merged file.
                            This option is only applicable to 'intersection', 'union' and 'consensus' modes.
                            Available values:
                             - min
                             - max
                             - sum
                             - mean
      -o OUTPUT, --output OUTPUT
                            Output junction file.  Required for operations that produce an output file.
      -p PREFIX, --prefix PREFIX
                            Prefix to apply to name column in BED output file
      -is, --ignore_strand  Whether or not to ignore strand when creating a key for the junction


Split
-----

This tool splits portcullis pass and fail juncs into 4 sets (TP, TN, FP, FN) based 
on whether or not the junctions are found in the reference.  The pass and fail files
passed into this tool should be disjoint in order to get meaningful results.

Usage::

    usage: Splits portcullis pass and fail juncs into 4 sets (TP, TN, FP, FN) based on whether or not the junctions are found in the reference or not.
           [-h] [-o OUTPUT_PREFIX] reference passfile failfile

    positional arguments:
      reference             The reference junction file
      passfile              The junction file containing junctions that pass a
                            filter
      failfile              The junction file containing junctions failing a
                            filter

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Prefix for output files
