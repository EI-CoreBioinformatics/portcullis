.. _scripts

Scripts
=======

Portcullis contains a number of python3 scripts for interpretting, analysing and converting junction files.  We won't list them all but here are some of the more useful ones:

gtf2bed.py and gff2bed.py
-------------------------

These scripts take existing annotations commonly available for many model organisms and then extract introns (junctions) from transcripts based on their exon composition, and then generate a BED12 format file from them.


tophat2bed12.py
---------------

Tophat generates a BED file containing junctions and anchors, however the actual splice junction sites (intron edges), have to be interpretted by combining the start and end coords with the block sizes and block starts.  As far as we can tell there's no reason to represent junctions in this way.  We store the actual coords of the splice sites (the bits we are most interested in) directly in the thick_start and thick_end columns.  This makes the file easier to interpret and will display just fine in all genome browsers we are aware of.


Other tool to bed scripts
-------------------------

Some other splice junction prediction tools don't directly generate bed files.  We provide some scripts to help do exactly this: `star_tab2bed.py`, `fs2bed.py` (finesplice), `spanki2bed.py`, `truesight2bed.py`, `soapsplice2.bed.py`.


bed_v_ref.py
----------

Allows the user to compare one or more BED files against a reference BED file.  This is useful for measuring the performance of portcullis, or any other tool against a reference.


bed_merge.py
------------

Combines a bunch of BED files into a single BED file.  This can be used to create a, union, intersection, or consensus BED file from multiple BED files.  For example, after running multiple RNAseq mapping tools on the same data, we can then pass the BAMs through `portcullis junc` to get the unfiltered junctions in BED format.  Then by running through this tool we can create a consensus between aligners which is often a highly accurate set of junctions.


Snakefiles
----------

We created a couple of snakemake pipelines for testing portcullis.  These are `read_gen.snakefile`, which is responsible for generating simulated datasets and `predict.snakefile`, which is responsible for driving multiple RNAseq mapping tools, splice junction prediction tools and transcript reconstruction tools.  These tools are not intended to be run directly by the user but might provide an insight into how to build pipelines efficiently.  We also provide a script called `snakey.py` which can drive any snakefile in a more convienient manner.  For more information, see the snakemake documentation online.

