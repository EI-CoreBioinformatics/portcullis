.. _scripts

Scripts
=======

Portcullis contains a number of python3 scripts for interpretting, analysing and converting junction files.  We won't list them all but here are some of the more useful ones:

portcullis_convert.py
---------------------

This contains a number of sub commands for converting between various formats that contain or include junctions.  See the help message within the commands for a full list of options, however some of the more commonly used options we list in the sections below.

Portcullis outputs a number of different files representing the junctions in different formats.  The most commonly used format (for us) downstream is the BED format.  The output directly from portcullis uses an exon-based BED format in order to describe both the junction start and end site plus the upstream and downstream anchors. The way we use the BED format is slightly different from tophat.  In the tophat format they choose not to explicitly define the junctions start and stop sites and instead rely on the block_start and block_size columns combined with the start and end columns to determine this.  For us, decided that it would be useful to explictly define the junction start and stop sites and we do this with the thick_start and thick_end columns.  Despite the difference in text format, both tophat and portcullis BED files look identical in a genome browser.

pbed2ibed
~~~~~~~~~

This mode of the conversion tool removes the anchors from the portcullis BED file, resulting in a pure intron-based BED file.  This means the start and end columns will be the same as the thick_start and thick_end columns.

tbed2pbed
~~~~~~~~~

Converts a tophat BED file into a portcullis BED file, retaining anchors.

gtf2ibed
~~~~~~~~

Takes an existing annotation or transcript assembly and then extract introns from transcripts based on their exon composition, and then generates a portcullis intron-based BED format file from them.


Others
~~~~~~

Some other splice junction prediction tools don't directly generate bed files.  We provide some additional modes to convert into an intron-based portcullis BED file.  Hopefully the names are self explanatory: `star2ibed`, `finesplice2ibed`, `spanki2ibed`, `truesight2ibed`, `soapsplice2ibed`.


bed_v_ref.py
----------

Allows the user to compare one or more BED files against a reference BED file.  This is useful for measuring the performance of portcullis, or any other tool against a reference.


bed_merge.py
------------

Combines a bunch of BED files into a single BED file.  This can be used to create a, union, intersection, or consensus BED file from multiple BED files.  For example, after running multiple RNAseq mapping tools on the same data, we can then pass the BAMs through `portcullis junc` to get the unfiltered junctions in BED format.  Then by running through this tool we can create a consensus between aligners which is often a highly accurate set of junctions.


Snakefiles
----------

We created a couple of snakemake pipelines for testing portcullis.  These are `read_gen.snakefile`, which is responsible for generating simulated datasets and `predict.snakefile`, which is responsible for driving multiple RNAseq mapping tools, splice junction prediction tools and transcript reconstruction tools.  These tools are not intended to be run directly by the user but might provide an insight into how to build pipelines efficiently.  We also provide a script called `snakey.py` which can drive any snakefile in a more convienient manner.  For more information, see the snakemake documentation online.

