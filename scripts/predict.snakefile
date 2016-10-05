import sys

from os import listdir
from os.path import isfile, join, abspath

from snakemake.utils import min_version

min_version("3.4.2")

R1 = config["r1"]
R2 = config["r2"]
REF = config["ref_fa"]
REF_GTF = config["ref_gtf"]
NAME = config["name"]
OUT_DIR = config["out_dir"]
OUT_DIR_FULL = os.path.abspath(config["out_dir"])
MIN_INTRON = config["min_intron"]
MAX_INTRON = config["max_intron"]
STRANDEDNESS = config["strandedness"]
THREADS = int(config["threads"])
READ_LENGTH = config["min_read_length"]
READ_LENGTH_MINUS_1 = int(READ_LENGTH) - 1

RUN_TIME = config["run_time"]


INDEX_STAR_EXTRA = config["index_star_extra"]
ALIGN_TOPHAT_EXTRA = config["align_tophat_extra"]
ALIGN_GSNAP_EXTRA = config["align_gsnap_extra"]
ALIGN_STAR_EXTRA = config["align_star_extra"]
ASM_CUFFLINKS_EXTRA = config["asm_cufflinks_extra"]
ASM_STRINGTIE_EXTRA = config["asm_stringtie_extra"]
JUNC_PORTCULLIS_PREP_EXTRA = config["junc_portcullis_prep_extra"]
JUNC_PORTCULLIS_JUNC_EXTRA = config["junc_portcullis_junc_extra"]


INPUT_SETS = config["input_sets"]
SIM_INPUT_SET = [x for x in INPUT_SETS if x.startswith("sim")]
REAL_INPUT_SET = [x for x in INPUT_SETS if x.startswith("real")]
ALIGNMENT_METHODS = config["align_methods"]
ASSEMBLY_METHODS = config["asm_methods"]
ASSEMBLY_MODES = config["asm_modes"]
JUNC_METHODS = config["junc_methods"]

#Shortcuts
READS_DIR = "."
READS_DIR_FULL = os.path.abspath(READS_DIR)
TRUE_DIR = OUT_DIR + "/true"
ALIGN_DIR = OUT_DIR + "/alignments"
ALIGN_DIR_FULL = os.path.abspath(ALIGN_DIR)
ASM_DIR = OUT_DIR + "/assemblies"
ASM_DIR_FULL = os.path.abspath(ASM_DIR)
JUNC_DIR = OUT_DIR + "/junctions"
JUNC_DIR_FULL = os.path.abspath(JUNC_DIR)
PORTCULLIS_DIR = JUNC_DIR + "/portcullis"
PORTCULLIS_DIR_FULL = os.path.abspath(PORTCULLIS_DIR)
PORTCULLIS_DIR2 = JUNC_DIR + "/portcullis2"
PORTCULLIS_DIR2_FULL = os.path.abspath(PORTCULLIS_DIR2)

CWD = os.getcwd()


HISAT_STRAND = "--rna-strandness=RF" if STRANDEDNESS == "fr-firststrand" else "--rna-strandness=FR" if STRANDEDNESS == "fr-secondstrand" else ""
PORTCULLIS_STRAND = "firststrand" if STRANDEDNESS == "fr-firststrand" else "secondstrand" if STRANDEDNESS == "fr-secondstrand" else "unstranded"


ISOFORM_FRACTION = {"permissive": 0.01, "semipermissive": 0.03, "default": 0.1, "semistrict": 0.2, "strict": 0.5}

#########################
Rules

localrules: all, truesight2bed, soapsplice2bed

# Define 
rule all:
	input: 
		expand(ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.stats", aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/{aln_method}-{reads}-portcullis.bed", aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/tophat-{reads}-spanki.bed", reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/tophat-{reads}-finesplice.bed", reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/truesight-{reads}-truesight.bed", reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/soapsplice-{reads}-soapsplice.bed", reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/mapsplice-{reads}-mapsplice.bed", reads=INPUT_SETS),
		expand(ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.bed", asm_method=ASSEMBLY_METHODS, asm_mode=ASSEMBLY_MODES, aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS),
		#expand(PORTCULLIS_DIR2+"/{aln_method}-{reads}/junc/{aln_method}-{reads}.junctions.tab", aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS)
#		expand(ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.stats", asm_method=ASSEMBLY_METHODS, asm_mode=ASSEMBLY_MODES, aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS)
		

#rule clean:
#	shell: "rm -rf {out}"

	

rule align_bowtie_index:
	input: REF
	output: 
		idx=ALIGN_DIR + "/bowtie/index/"+NAME+".3.ebwt",
		fa_link=ALIGN_DIR+"/bowtie/index/"+NAME+".fa"
	params:
		out_prefix=ALIGN_DIR+"/bowtie/index/"+NAME,
		ref=os.path.abspath(REF),
		load=config["load"]["bowtie"]
	log: ALIGN_DIR + "/bowtie.index.log"
	threads: 1
	message: "Indexing genome with bowtie"
	shell: "{params.load} && bowtie-build {input} {params.out_prefix} > {log} 2>&1 && ln -sf {params.ref} {output.fa_link} && touch -h {output.fa_link}"


rule align_tophat_index:
	input: REF
	output: ALIGN_DIR +"/tophat/index/"+NAME+".4.bt2"
	params: load=config["load"]["tophat"]
	log: ALIGN_DIR + "/tophat.index.log"
	threads: 1
	message: "Indexing genome with tophat"
	shell: "{params.load} && bowtie2-build {REF} {ALIGN_DIR}/tophat/index/{NAME} > {log} 2>&1"


rule align_gsnap_index:
	input: REF
	output: ALIGN_DIR +"/gsnap/index/"+NAME+"/"+NAME+".sachildguide1024"
	params: load=config["load"]["gmap"]
	log: ALIGN_DIR +"/gsnap.index.log"
	threads: 1
	message: "Indexing genome with gsnap"
	shell: "{params.load} && gmap_build --dir={ALIGN_DIR}/gsnap/index --db={NAME} {input} > {log} 2>&1"


rule align_star_index:
	input: os.path.abspath(REF)
	output: ALIGN_DIR +"/star/index/SAindex"
	params: 
		indexdir=ALIGN_DIR_FULL+"/star/index",
		load=config["load"]["star"]
	log: ALIGN_DIR_FULL+"/star.index.log"
	threads: int(THREADS)
	message: "Indexing genome with star"
	shell: "{params.load} && cd {ALIGN_DIR_FULL}/star && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input} {INDEX_STAR_EXTRA} > {log} 2>&1 && cd {CWD}"



rule align_hisat_index:
	input: REF
	output: ALIGN_DIR+"/hisat/index/"+NAME+".4.ht2"
	params: load=config["load"]["hisat"]
	log: ALIGN_DIR+"/hisat.index.log"
	threads: 1
	message: "Indexing genome with hisat"
	shell: "{params.load} && hisat2-build {input} {ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1"





#####


rule align_tophat:
	input:
		r1=READS_DIR+"/{reads}.R1.fq",
		r2=READS_DIR+"/{reads}.R2.fq",
		index=rules.align_tophat_index.output
	output:
		bam=ALIGN_DIR_FULL+"/tophat/{reads}/accepted_hits.bam",
		link=ALIGN_DIR+"/output/tophat-{reads}.bam"
	params:
		load=config["load"]["tophat"],
		outdir=ALIGN_DIR+"/tophat/{reads}",
		indexdir=ALIGN_DIR+"/tophat/index/"+NAME,
		strand=lambda wildcards: STRANDEDNESS if wildcards.reads.startswith("real") else "fr-unstranded"
	log: ALIGN_DIR + "/tophat-{reads}.log"
	threads: THREADS
	message: "Aligning RNAseq data with tophat"
	shell: "{params.load} && {RUN_TIME} tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={params.strand} {ALIGN_TOPHAT_EXTRA} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {output.bam} {output.link} && touch -h {output.link}"



rule align_gsnap:
	input:
		r1=READS_DIR+"/{reads}.R1.fq",
		r2=READS_DIR+"/{reads}.R2.fq",
		index=rules.align_gsnap_index.output
	output:
		bam=ALIGN_DIR_FULL+"/gsnap/{reads}/gsnap.bam",
		link=ALIGN_DIR+"/output/gsnap-{reads}.bam"
	params:
		load_g=config["load"]["gmap"],
		load_s=config["load"]["samtools"]
	log: ALIGN_DIR+"/gsnap-{reads}.log"
	threads: THREADS
	message: "Aligning RNAseq with gsnap"
	shell: "{params.load_g} && {params.load_s} && {RUN_TIME} gsnap --dir={ALIGN_DIR}/gsnap/index --db={NAME} {ALIGN_GSNAP_EXTRA} --novelsplicing=1 --localsplicedist={MAX_INTRON} --nthreads={threads} --format=sam --npaths=20 {input.r1} {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {output.bam} {output.link} && touch -h {output.link}"




rule align_star:
	input:
		r1=READS_DIR_FULL+"/{reads}.R1.fq",
		r2=READS_DIR_FULL+"/{reads}.R2.fq",
		index=rules.align_star_index.output
	output:
		bam=ALIGN_DIR_FULL+"/star/{reads}/Aligned.out.bam",
		link=ALIGN_DIR+"/output/star-{reads}.bam"
	params:
		load=config["load"]["star"],
		outdir=ALIGN_DIR_FULL+"/star/{reads}",
		indexdir=ALIGN_DIR_FULL+"/star/index"
	log: ALIGN_DIR_FULL+"/star-{reads}.log"
	threads: int(THREADS)
	message: "Aligning input with star"
	shell: "{params.load} && cd {params.outdir} && {RUN_TIME} STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.indexdir} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON} --alignMatesGapMax 20000 --outFileNamePrefix {params.outdir}/ {ALIGN_STAR_EXTRA} > {log} 2>&1 && cd {CWD} && ln -sf {output.bam} {output.link} && touch -h {output.link}"



rule align_hisat:
	input:
		r1=READS_DIR+"/{reads}.R1.fq",
		r2=READS_DIR+"/{reads}.R2.fq",
		index=rules.align_hisat_index.output
	output:
		bam=ALIGN_DIR_FULL+"/hisat/{reads}/hisat.bam",
		link=ALIGN_DIR+"/output/hisat-{reads}.bam"
	params:
		load_h=config["load"]["hisat"],
		load_s=config["load"]["samtools"],
		indexdir=ALIGN_DIR+"/hisat/index/"+NAME,
		strand=lambda wildcards: HISAT_STRAND if wildcards.reads.startswith("real") else ""
	log: ALIGN_DIR+"/hisat-{reads}.log"
	threads: THREADS
	message: "Aligning input with hisat"
	shell: "{params.load_h} && {params.load_s} && {RUN_TIME} hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {strand} -x {params.indexdir} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {output.bam} {output.link} && touch -h {output.link}"

rule bam_sort:
	input: 
		bam=ALIGN_DIR+"/output/{aln_method}-{reads}.bam"
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam"
	params:
		load_s=config["load"]["samtools"]		
	threads: THREADS
	message: "Using samtools to sort {input.bam}"
	shell: "{params.load_s} && samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.aln_method}_{wildcards.reads} -@ {threads} {input.bam}"


rule bam_index:
	input: rules.bam_sort.output
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.bai"
	params:
		load_s=config["load"]["samtools"]		
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{params.load_s} && samtools index {input}"

rule bam_stats:
	input:
		bam=rules.bam_sort.output,
		idx=rules.bam_index.output
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.stats"
	params:
		load=config["load"]["samtools"],
		plot_out=ALIGN_DIR+"/output/plots/{aln_method}-{reads}/{aln_method}-{reads}"
	threads: 1
	message: "Using samtools to collected stats for: {input}"
	shell: "{params.load} && samtools stats {input.bam} > {output} && plot-bamstats -p {params.plot_out} {output}"

rule portcullis_prep:
	input:
		ref=REF,
		bam=rules.bam_sort.output,
		idx=rules.bam_index.output
	output: 
		bai=PORTCULLIS_DIR+"/{aln_method}-{reads}/prep/portcullis.sorted.alignments.bam.bai",
		fa=PORTCULLIS_DIR+"/{aln_method}-{reads}/prep/portcullis.genome.fa"
	params:
		outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/prep",
		load=config["load"]["portcullis"],
		strand=lambda wildcards: PORTCULLIS_STRAND if wildcards.reads.startswith("real") else "unstranded"
	log: PORTCULLIS_DIR+"/{aln_method}-{reads}-prep.log"
	threads: THREADS
	message: "Using portcullis to prepare: {input}"
	shell: "{params.load} && {RUN_TIME} portcullis prep -o {params.outdir} --strandedness={strand} -t {threads} {input.ref} {input.bam} > {log} 2>&1"


rule portcullis_junc:
	input:
		bai=rules.portcullis_prep.output
	output: PORTCULLIS_DIR+"/{aln_method}-{reads}/junc/{aln_method}-{reads}.junctions.tab"
	params:
		prepdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/prep",
		outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/junc",
		load=config["load"]["portcullis"],
		strand=lambda wildcards: PORTCULLIS_STRAND if wildcards.reads.startswith("real") else "unstranded"
	log: PORTCULLIS_DIR+"/{aln_method}-{reads}-junc.log"
	threads: THREADS
	message: "Using portcullis to analyse potential junctions: {input}"
	shell: "{params.load} && {RUN_TIME} portcullis junc -o {params.outdir}/{wildcards.aln_method}-{wildcards.reads} --strandedness={strand} -t {threads} {params.prepdir} > {log} 2>&1"


rule portcullis_filter:
	input: 	tab=rules.portcullis_junc.output,
		ref=rules.portcullis_prep.output.fa
	output:
		link=JUNC_DIR+"/output/{aln_method}-{reads}-portcullis.bed",
		unfilt_link=JUNC_DIR+"/output/{aln_method}-{reads}-all.bed",
		tab=PORTCULLIS_DIR+"/{aln_method}-{reads}/filt/{aln_method}-{reads}.pass.junctions.tab",
	params:
		outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/filt",
		load=config["load"]["portcullis"],
		bed=PORTCULLIS_DIR_FULL+"/{aln_method}-{reads}/filt/{aln_method}-{reads}.pass.junctions.bed",
		unfilt_bed=PORTCULLIS_DIR_FULL+"/{aln_method}-{reads}/junc/{aln_method}-{reads}.junctions.bed"
	log: PORTCULLIS_DIR+"/{aln_method}-{reads}-filter.log"
	threads: THREADS
	message: "Using portcullis to filter invalid junctions: {input}"
	shell: "{params.load} && portcullis filter -t {threads} -o {params.outdir}/{wildcards.aln_method}-{wildcards.reads} {input.ref} {input.tab} > {log} 2>&1 && ln -sf {params.bed} {output.link} && touch -h {output.link} && ln -sf {params.unfilt_bed} {output.unfilt_link} && touch -h {output.unfilt_link}"


rule portcullis_bamfilt:
	input: 
		bam=rules.bam_sort.output,
		tab=rules.portcullis_filter.output.tab
	output:
		bam=PORTCULLIS_DIR+"/{aln_method}-{reads}/bam/{aln_method}-{reads}-portcullis.bam"
	params:
		load=config["load"]["portcullis"]
	log: PORTCULLIS_DIR+"/{aln_method}-{reads}-bam.log"
	threads: THREADS
	message: "Using portcullis to filter alignments containing invalid junctions: {input.tab}"
	shell: "{params.load} && portcullis bamfilt --output={output.bam} {input.tab} {input.bam} > {log} 2>&1"




rule spanki:
    input:
        bam=rules.bam_sort.output,
        fa=REF,
        gtf=REF_GTF,
        idx=rules.bam_index.output
    output: link=JUNC_DIR+"/output/{aln_method}-{reads}-spanki.bed",
        bed=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki.bed"
    params:
        load_spanki=config["load"]["spanki"],
        load_portcullis=config["load"]["portcullis"],
        outdir=JUNC_DIR_FULL+"/spanki/{aln_method}-{reads}",
        bam=ALIGN_DIR_FULL+"/output/{aln_method}-{reads}.sorted.bam",
        fa=os.path.abspath(REF),
        gtf=os.path.abspath(REF_GTF),
        all_juncs=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/juncs.all",
        filt_juncs=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/juncs.filtered",
        bed=JUNC_DIR_FULL+"/spanki/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki.bed"
    log: JUNC_DIR_FULL+"/spanki/{aln_method}-{reads}-spanki.log"
    threads: 1
    message: "Using SPANKI to analyse junctions: {input.bam}"
    shell: "set +e && {params.load_spanki} && cd {params.outdir} && {RUN_TIME} spankijunc -i {params.bam} -g {params.gtf} -f {params.fa} > {log} 2>&1 && cd {CWD} && if [[ -s {params.all_juncs} ]] ; then {params.load_portcullis} && spanki_filter.py {params.all_juncs} > {params.filt_juncs} && portcullis_convert spanki2ibed {params.filt_juncs} > {output.bed} ; else touch {output.bed} ; fi && ln -sf {params.bed} {output.link} && touch -h {output.link}"

rule spanki_annot:
    input:
        bam=rules.bam_sort.output,
        fa=REF,
        gtf=REF_GTF,
        idx=rules.bam_index.output
    output: link=JUNC_DIR+"/output/{aln_method}-{reads}-spanki_annot.bed",
        bed=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki_annot.bed"
    params:
        load_spanki=config["load"]["spanki"],
        load_portcullis=config["load"]["portcullis"],
        outdir=JUNC_DIR_FULL+"/spanki_annot/{aln_method}-{reads}",
        bam=ALIGN_DIR_FULL+"/output/{aln_method}-{reads}.sorted.bam",
        fa=os.path.abspath(REF),
        gtf=os.path.abspath(REF_GTF),
        all_juncs=JUNC_DIR+"/spanki_annot/{aln_method}-{reads}/junctions_out/juncs.all",
        filt_juncs=JUNC_DIR+"/spanki_annot/{aln_method}-{reads}/junctions_out/juncs.filtered",
        bed=JUNC_DIR_FULL+"/spanki_annot/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki_annot.bed"
    log: JUNC_DIR_FULL+"/spanki_annot/{aln_method}-{reads}-spanki_annot.log"
    threads: 1
    message: "Using SPANKI to analyse junctions: {input.bam}"
    shell: "set +e && {params.load_spanki} && cd {params.outdir} && {RUN_TIME} spankijunc -i {params.bam} -g {params.gtf} -f {params.fa} -filter T > {log} 2>&1 && cd {CWD} && if [[ -s {params.all_juncs} ]] ; then {params.load_portcullis} && portcullis_convert spanki2ibed.py {params.all_juncs} > {output.bed} ; else touch {output.bed} ; fi && ln -sf {params.bed} {output.link} && touch -h {output.link}"
	
rule finesplice:
    input:
        bam=rules.bam_sort.output,
        idx=rules.bam_sort.output
    output: link=JUNC_DIR+"/output/{aln_method}-{reads}-finesplice.bed",
        bed=JUNC_DIR+"/finesplice/{aln_method}-{reads}/{aln_method}-{reads}-finesplice.bed"
    params:
        load_fs=config["load"]["finesplice"],
        load_portcullis=config["load"]["portcullis"],
        bam=ALIGN_DIR_FULL+"/output/{aln_method}-{reads}.sorted.bam",
        outdir=JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}",
        junc=JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}/{aln_method}-{reads}.sorted.accepted.junc",
        bed=JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}/{aln_method}-{reads}-finesplice.bed"
    log: JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}-finesplice.log"
    threads: 1
    message: "Using FineSplice to analyse junctions: {input.bam}"
    shell: "set +e && {params.load_fs} && cd {params.outdir} && if {RUN_TIME} FineSplice.py -i {params.bam} -l {READ_LENGTH} > {log} 2>&1 ; then cd {CWD} && {params.load_portcullis} && portcullis_convrt finesplice2ibed.py {params.junc} > {output.bed} ; else cd {CWD}; touch {output.bed} ; fi && ln -sf {params.bed} {output.link} && touch -h {output.link}"


rule truesight:
    input:
        idx=rules.align_bowtie_index.output,
        r1=READS_DIR+"/{reads}.R1.fq",
        r2=READS_DIR+"/{reads}.R2.fq"
    output: JUNC_DIR+"/truesight/{reads}/GapAli.junc"
    params:
        load_fs=config["load"]["truesight"],
	index=ALIGN_DIR_FULL + "/bowtie/index/"+NAME,
        outdir=JUNC_DIR+"/truesight/{reads}",
        junc=JUNC_DIR+"/truesight/{reads}/GapAli.junc",
	r1=READS_DIR_FULL+"/{reads}.R1.fq",
	r2=READS_DIR_FULL+"/{reads}.R2.fq"
    log: JUNC_DIR_FULL+"/truesight/{reads}-truesight.log"
    threads: int(THREADS)
    message: "Using Truesight to find junctions"
    shell: "{params.load_fs} && cd {params.outdir} && {RUN_TIME} truesight_pair.pl -i {MIN_INTRON} -I {MAX_INTRON} -v 1 -r {params.index} -p {threads} -o . -f {params.r1} {params.r2} > {log} 2>&1 && rm GapAli.sam && cd {CWD}"

rule truesight2bed:
	input: rules.truesight.output
	output:
		bed=JUNC_DIR+"/output/truesight-{reads}-truesight.bed"
	params:
		load_portcullis=config["load"]["portcullis"]
	threads: 1
	message: "Creating bed file from truesight output: {input}"
	shell: "{params.load_portcullis} && portcullis_convert.py truesight2ibed {input} > {output.bed}"

rule soapsplice_index:
	input: fa=REF
	output: JUNC_DIR + "/soapsplice/index/"+NAME+".index.bwt"
	log: JUNC_DIR + "/soapsplice/index.log"
	params:
        	load_ss=config["load"]["soapsplice"],
		index=JUNC_DIR + "/soapsplice/index/"+NAME
	threads: 1
	message: "Creating index for soapsplice"
	shell: "{params.load_ss} && 2bwt-builder {input.fa} {params.index} > {log} 2>&1"

rule soapsplice:
	input:
		r1=READS_DIR+"/{reads}.R1.fq",
		r2=READS_DIR+"/{reads}.R2.fq",
		idx=rules.soapsplice_index.output
	output: JUNC_DIR+"/soapsplice/{reads}/ss-{reads}.junc"
	params:
		load_ss=config["load"]["soapsplice"],
		index=JUNC_DIR + "/soapsplice/index/"+NAME+".index",
		outdir=JUNC_DIR+"/soapsplice/{reads}"
	log: JUNC_DIR+"/soapsplice/{reads}-soapsplice.log"
	threads: THREADS
	message: "Using soapsplice to find junctions"
	shell: "{params.load_ss} && {RUN_TIME} soapsplice -d {params.index} -1 {input.r1} -2 {input.r2} -I 450 -o {params.outdir}/ss-{wildcards.reads} -p {threads} -t {MAX_INTRON} -c 0 -f 2 -L {MAX_INTRON} -l {MIN_INTRON} > {log} 2>&1 && rm -f {params.outdir}/ss-{wildcards.reads}.sam"

rule soapsplice2bed:
	input: rules.soapsplice.output
	output:
		bed=JUNC_DIR+"/output/soapsplice-{reads}-soapsplice.bed"
	params:
		load_portcullis=config["load"]["portcullis"]
	threads: 1
	message: "Creating bed file from soapsplice output: {input}"
	shell: "{params.load_portcullis} && portcullis_convert.py soapsplice2ibed {input} > {output.bed}"

rule mapsplice_ref:
	input:
		ref=REF
	output:
		done=JUNC_DIR+"/mapsplice/ref/all.done"
	params:
		load=config["load"]["portcullis"],
		outdir=JUNC_DIR+"/mapsplice/ref"
	message: "Creating reference for mapsplice"
	threads: 1
	shell: "{params.load} && split.py -o {params.outdir} {input.ref} && touch {output.done}"

rule mapsplice:
	input:
		r1=READS_DIR+"/{reads}.R1.fq",
		r2=READS_DIR+"/{reads}.R2.fq",
		ref=rules.mapsplice_ref.output.done
	output:
		JUNC_DIR+"/mapsplice/{reads}/junctions.txt"
	params:
		refdir=JUNC_DIR+"/mapsplice/ref",
		idx=ALIGN_DIR+"/bowtie/index/"+NAME,
		load_ms=config["load"]["mapsplice"],
		outdir=JUNC_DIR+"/mapsplice/{reads}"
	log: JUNC_DIR+"/mapsplice/{reads}-mapsplice.log"
	threads: int(THREADS)
	message: "Using mapsplice to find junctions"
	shell: "{params.load_ms} && {RUN_TIME} mapsplice.py -c {params.refdir} -1 {input.r1} -2 {input.r2} -o {params.outdir}/mapsplice-{wildcards.reads} -p {threads} --bam -i {MIN_INTRON} -I {MAX_INTRON} > {log} 2>&1"

rule mapsplice2bed:
	input: rules.mapsplice.output
	output:
		bed=JUNC_DIR+"/output/mapsplice-{reads}-mapsplice.bed"
	params:
		load_p=config["load"]["portcullis"],
	threads: 1
	message: "Creating bed file from mapsplice output: {input}"
	shell: "{params.load_p} && portcullis_convert.py mapsplice2ibed {input} > {output.bed}"


###

rule asm_cufflinks:
	input:
		bam=rules.bam_sort.output,
		ref=REF
	output: 
		gtf=ASM_DIR+"/output/cufflinks_{asm_mode}-{aln_method}-{reads}.gtf"
	params:
		outdir=ASM_DIR+"/cufflinks_{asm_mode}-{aln_method}-{reads}",
		gtf=ASM_DIR+"/cufflinks_{asm_mode}-{aln_method}-{reads}/transcripts.gtf",
		link_src="../cufflinks_{asm_mode}-{aln_method}-{reads}/transcripts.gtf",
		load=config["load"]["cufflinks"],
		iso_frac=lambda wildcards: ISOFORM_FRACTION[wildcards.asm_mode]
	log: ASM_DIR+"/cufflinks_{asm_mode}-{aln_method}-{reads}.log"
	threads: int(THREADS)
	message: "Using cufflinks to assemble {input.bam}"
	shell: "{params.load} && cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} -F {params.iso_frac} --no-update-check {input.bam} > {log} 2>&1 && ln -sf {params.link_src} {output.gtf} && touch -h {output.gtf}"



rule asm_stringtie:
	input: 	bam=rules.bam_sort.output
	output: 
		link=ASM_DIR+"/output/stringtie_{asm_mode}-{aln_method}-{reads}.gtf",
		gtf=ASM_DIR+"/stringtie_{asm_mode}-{aln_method}-{reads}/stringtie_{asm_mode}-{aln_method}-{reads}.gtf"
	params:
		load=config["load"]["stringtie"],
		gtf=ASM_DIR+"/stringtie_{asm_mode}-{aln_method}-{reads}/stringtie_{asm_mode}-{aln_method}-{reads}.gtf",
		link_src="../stringtie_{asm_mode}-{aln_method}-{reads}/stringtie_{asm_mode}-{aln_method}-{reads}.gtf",
		name="Stringtie_{asm_mode}_{aln_method}_{reads}",
		iso_frac=lambda wildcards: ISOFORM_FRACTION[wildcards.asm_mode]
	log: ASM_DIR+"/stringtie_{asm_mode}-{aln_method}-{reads}.log"
	threads: int(THREADS)
	message: "Using stringtie to assemble: {input.bam}"
	shell: "{params.load} && {RUN_TIME} stringtie {input.bam} -l {params.name} -f {params.iso_frac} -m 200 -o {params.gtf} -p {threads} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule asm_class:
	input:
		bam=rules.bam_sort.output,
		ref=REF
	output:
		link=ASM_DIR+"/output/class_{asm_mode}-{aln_method}-{reads}.gtf",
		gtf=ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}/class_{asm_mode}-{aln_method}-{reads}.gtf"
	params:
		outdir=ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}",
		load=config["load"]["class"],
		gtf=ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}/class_{asm_mode}-{aln_method}-{reads}.gtf",
		link_src="../class_{asm_mode}-{aln_method}-{reads}/class_{asm_mode}-{aln_method}-{reads}.gtf",
		iso_frac=lambda wildcards: ISOFORM_FRACTION[wildcards.asm_mode]
	log: ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}.log"
	threads: int(THREADS)
	message: "Using class to assemble: {input.bam}"
	shell: "{params.load} && class_run.py -c '-F {params.isofrac}' -p {threads} {input.bam} > {output.gtf} 2> {log} && ln -sf {params.link_src} {output.link} && touch -h {output.link}"


rule gtf_2_bed:
	input:
		gtf=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.gtf"
	output:
		bed=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.bed"
	params:
		load_p=config["load"]["portcullis"]
	message: "Converting GTF to BED for: {input.gtf}"
	shell: "{params.load_p} && portcullis_convert.py gtf2ibed {input.gtf} > {output.bed}"

rule gtf_stats:
	input:
		gtf=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.gtf",
		ref=REF_GTF
	output: 
		comp=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.comp.stats",
		stats=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.stats"
	params: 
		load_mikado=config["load"]["mikado"],
		prefix=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.comp"
	message: "Calculating stats for: {input.gtf}"
	shell: "{params.load_mikado} && mikado.py util stats {input.gtf} > {output.stats} && mikado.py compare -r {input.ref} -p {input.gtf} -o {params.prefix}"


