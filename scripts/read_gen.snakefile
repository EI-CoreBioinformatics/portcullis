import sys

from os import listdir
from os.path import isfile, join, abspath

from snakemake.utils import min_version

min_version("3.4.2")

R1 = config["r1"]
R2 = config["r2"]
R1_NAME, ext1 = os.path.splitext(os.path.basename(R1))
R2_NAME, ext2 = os.path.splitext(os.path.basename(R2))
REF = config["ref_fa"]
REF_GTF = config["ref_gtf"]
NAME = config["name"]
OUT_DIR = config["out_dir"]
OUT_DIR_FULL = os.path.abspath(config["out_dir"])
MIN_INTRON = config["min_intron"]
MAX_INTRON = config["max_intron"]
STRANDEDNESS = config["strandedness"]
THREADS = config["threads"]
READ_LENGTH = config["min_read_length"]
READ_LENGTH_MINUS_1 = int(READ_LENGTH) - 1

LOAD_BOWTIE = config["load_bowtie"]
LOAD_FASTQC = config["load_fastqc"]
LOAD_SPANKI = config["load_spanki"]
LOAD_TRIMGALORE = config["load_trimgalore"]
LOAD_TOPHAT = config["load_tophat"]
LOAD_GMAP = config["load_gsnap"]
LOAD_STAR = config["load_star"]
LOAD_HISAT = config["load_hisat"]
LOAD_SAMTOOLS = config["load_samtools"]
LOAD_CUFFLINKS = config["load_cufflinks"]
LOAD_STRINGTIE = config["load_stringtie"]
LOAD_PORTCULLIS = config["load_portcullis"]
LOAD_PYTHON3 = config["load_python3"]

INDEX_STAR_EXTRA = config["index_star_extra"]
ALIGN_TOPHAT_EXTRA = config["align_tophat_extra"]
ALIGN_GSNAP_EXTRA = config["align_gsnap_extra"]
ALIGN_STAR_EXTRA = config["align_star_extra"]
ASM_CUFFLINKS_EXTRA = config["asm_cufflinks_extra"]
ASM_STRINGTIE_EXTRA = config["asm_stringtie_extra"]
JUNC_PORTCULLIS_PREP_EXTRA = config["junc_portcullis_prep_extra"]
JUNC_PORTCULLIS_JUNC_EXTRA = config["junc_portcullis_junc_extra"]


INPUT_SETS = ["real","sim_fixed","sim_var"]
ALIGNMENT_METHODS = ["tophat","star","gsnap","hisat"]
ALIGNMENT_WREF_METHODS = ["tophat","star"]
ASSEMBLY_METHODS = config["asm_methods"]

#Shortcuts
READS_DIR = OUT_DIR + "/reads"
READS_DIR_FULL = os.path.abspath(READS_DIR)
TRUE_DIR = OUT_DIR + "/true"
SIM_DIR = READS_DIR + "/sim"
SIM_DIR_FULL = os.path.abspath(SIM_DIR)
ALIGN_DIR = OUT_DIR + "/alignments"
ALIGN_DIR_FULL = os.path.abspath(ALIGN_DIR)
ASM_DIR = OUT_DIR + "/assemblies"
ASM_DIR_FULL = os.path.abspath(ASM_DIR)
PORT_DIR = OUT_DIR + "/portcullis"
PORT_DIR_FULL = os.path.abspath(PORT_DIR)

CWD = os.getcwd()


HISAT_STRAND = "--rna-strandness=RF" if STRANDEDNESS == "fr-firststrand" else "--rna-strandness=FR" if STRANDEDNESS == "fr-secondstrand" else "F"
PORTCULLIS_STRAND = "firststrand" if STRANDEDNESS == "fr-firststrand" else "secondstrand" if STRANDEDNESS == "fr-secondstrand" else "unstranded"


#########################
Rules

localrules: all, asm_cufflinks_cov

# Define 
rule all:
	input: 
		SIM_DIR_FULL+"/var/sim.bam",
		SIM_DIR_FULL+"/fixed/sim.bam",
		PORT_DIR+"/junc/portcullis_sim_var.junctions.bed",
		READS_DIR+"/trim/fastqc/real.R1_fastqc.html"
		

rule trim_reads:
	input: 
		r1=R1,
		r2=R2
	output:
		linkr1=READS_DIR+"/real.R1.fq",
		linkr2=READS_DIR+"/real.R2.fq"
	params:
		outdir=READS_DIR+"/trim",
		r1="trim/"+R1_NAME+"_val_1.fq",
		r2="trim/"+R2_NAME+"_val_2.fq",
	log: READS_DIR+"/trim/trim.log"		
	shell: "{LOAD_TRIMGALORE} && trim_galore --paired --length {READ_LENGTH} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"


rule fastqc:
	input:
		r1=rules.trim_reads.output.linkr1,
		r2=rules.trim_reads.output.linkr2
	output: 
		r1=READS_DIR+"/trim/fastqc/real.R1_fastqc.html",
		r2=READS_DIR+"/trim/fastqc/real.R2_fastqc.html"
	params:
		outdir=READS_DIR+"/trim/fastqc"
	log: READS_DIR+"/trim/fastqc.log"
	shell: "{LOAD_FASTQC} && fastqc -o {params.outdir} {input.r1} {input.r2}"


rule simgen_bowtie_index:
	input:
		REF
	output: 
		SIM_DIR+"/bowtie_index/"+REF+".1.ebwt"
	params: outdir=SIM_DIR+"/bowtie_index/"+REF
	log:
		SIM_DIR+"/bowtie_index/bowtie_index.log"
	threads: 1
	message: "Creating bowtie index of genome"
	shell: "{LOAD_BOWTIE} && bowtie-build {input} {params.outdir}"


rule simgen_bowtie_align:
	input:
		r1=rules.trim_reads.output.linkr1,
                r2=rules.trim_reads.output.linkr2,
                index=rules.simgen_bowtie_index.output
	output:
		map=SIM_DIR+"/bowtie/bowtie.map"
	params: indexdir=SIM_DIR+"/bowtie_index/"+REF
	log: SIM_DIR+"/bowtie.align.log"
		
	threads: int(THREADS)
	shell: "{LOAD_BOWTIE} && bowtie -p {threads} --chunkmbs 1000 -X 500 {params.indexdir} -1 {input.r1} -2 {input.r2} > {output.map} 2> {log}"
	


rule align_tophat_index:
        input: REF 
        output: ALIGN_DIR +"/tophat/index/"+NAME+".4.bt2"
        log: ALIGN_DIR + "/tophat.index.log"
        threads: 1
        message: "Indexing genome with tophat"
        shell: "{LOAD_TOPHAT} && bowtie2-build {REF} {ALIGN_DIR}/tophat/index/{NAME} > {log} 2>&1"




rule align_tophat_wref:
	input:
		r1=rules.trim_reads.output.linkr1,
		r2=rules.trim_reads.output.linkr2,
		index=rules.align_tophat_index.output,
		gtf=REF_GTF
	output:
		bam=ALIGN_DIR+"/output/tophat-real_wref.bam",
		bed=ALIGN_DIR+"/tophat_wref/junctions.bed"
	params:
		outdir=ALIGN_DIR+"/tophat_wref/",
		indexdir=ALIGN_DIR+"/tophat/index/"+NAME,
		bam="../tophat_wref/accepted_hits.bam",
	log: ALIGN_DIR + "/tophat-real_wref.log"
	threads: int(THREADS)
	message: "Aligning RNAseq data with tophat using GTF reference: {input.r1} {input.r2} {input.gtf}"
	shell: "{LOAD_TOPHAT} && tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={STRANDEDNESS} --GTF={input.gtf} {ALIGN_TOPHAT_EXTRA} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {params.bam} {output.bam} && touch -h {output.bam}"


rule bam_sort_wref:
	input: ALIGN_DIR+"/output/{align_ref_method}-real_wref.bam"
	output: ALIGN_DIR+"/output/{align_ref_method}-real_wref.sorted.bam"
	threads: int(THREADS)
	message: "Using samtools to sort {input}"
	shell: "{LOAD_SAMTOOLS} && samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.align_ref_method}_wref -@ {threads} {input}"

rule bam_index_wref:
	input: rules.bam_sort_wref.output
	output: ALIGN_DIR+"/output/{align_ref_method}-real_wref.sorted.bam.bai"
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{LOAD_SAMTOOLS} && samtools index {input}"

rule asm_cufflinks_wref:
	input:
		bam=ALIGN_DIR+"/output/tophat-real_wref.sorted.bam",
		ref=REF_GTF
	output:	ASM_DIR+"/cufflinks-tophat_wref-real/isoforms.fpkm_tracking"
	params:	outdir=ASM_DIR+"/cufflinks-tophat_wref-real",
	log: ASM_DIR+"/cufflinks-tophat_wref-real.log"
	threads: int(THREADS)
	message: "Using cufflinks to assemble: {input.bam}"
	shell: "{LOAD_CUFFLINKS} && cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} -G {input.ref} -g {input.ref} --no-update-check {ASM_CUFFLINKS_EXTRA} {input.bam} > {log} 2>&1"

rule asm_cufflinks_cov:
	input: rules.asm_cufflinks_wref.output
	output: ASM_DIR+"/cufflinks-tophat_wref-real/transcripts.cov"
	params: fpkminref=ASM_DIR+"/cufflinks-tophat_wref-real/isoforms.fpkm_tracking_inref"
	threads: 1
	message: "Creating transcript coverage output file"
	shell: "{LOAD_PORTCULLIS} && {LOAD_PYTHON3} && grep -v ^CUFF {input} > {params.fpkminref} && cufflinks2spankicov.py {params.fpkminref} > {output}"

rule simgen_model:
	input:
		map=rules.simgen_bowtie_align.output.map
	output: SIM_DIR+"/model/logfile.txt"
	params: outdir=SIM_DIR+"/model"
	log: SIM_DIR+"/spanki_model.log"
	threads: 1
	shell: "{LOAD_SPANKI} && spankisim_models -i {input.map} -e 2 -l {READ_LENGTH} -o {params.outdir} > {log} 2>&1"


rule sim_fixed_reads:
	input:
		gtf=REF_GTF,
                fa=REF,
		model=rules.simgen_model.output
	output: 
		bam=SIM_DIR_FULL+"/fixed/sim.bam",
		linkr1=READS_DIR+"/sim_fixed.R1.fq",
		linkr2=READS_DIR+"/sim_fixed.R2.fq"
	params: 
		outdir=SIM_DIR+"/fixed",
		mdir=SIM_DIR+"/model",
		r1=SIM_DIR_FULL+"/fixed/sim_1.fastq",
		r2=SIM_DIR_FULL+"/fixed/sim_2.fastq",
	log: SIM_DIR+"/spanki_readgen_fixed.log"
	threads: 1	
	shell: "{LOAD_SPANKI} && {LOAD_CUFFLINKS} && spankisim_transcripts -g {input.gtf} -f {input.fa} -o {params.outdir} -m custom -mdir {params.mdir} -cov 30 -bp {READ_LENGTH} -ends 2 -frag 250 > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"
	


rule sim_var_reads:
	input:
		gtf=REF_GTF,
                fa=REF,
                transcript_cov=rules.asm_cufflinks_cov.output,
		model=rules.simgen_model.output
	output:
		bam=SIM_DIR_FULL+"/var/sim.bam",
		linkr1=READS_DIR+"/sim_var.R1.fq",
		linkr2=READS_DIR+"/sim_var.R2.fq"
	params: 
		outdir=SIM_DIR+"/var",
		mdir=SIM_DIR+"/model",
		r1=SIM_DIR_FULL+"/var/sim_1.fastq",
		r2=SIM_DIR_FULL+"/var/sim_2.fastq",
	log: SIM_DIR+"/spanki_readgen_var.log"
	threads: 1	
	shell: "{LOAD_SPANKI} && {LOAD_CUFFLINKS} && spankisim_transcripts -g {input.gtf} -f {input.fa} -o {params.outdir} -m custom -mdir {params.mdir} -t {input.transcript_cov} -bp {READ_LENGTH} -ends 2 -frag 250 > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"

rule sim_var_portcullis_prep:
	input:
		bam=rules.sim_var_reads.output.bam,
		ref=REF,
	output: PORT_DIR+"/prep/portcullis.sorted.alignments.bam.bai"
	params:
		outdir=PORT_DIR+"/prep",
		load=LOAD_PORTCULLIS
	log: PORT_DIR+"/prep.log"
	threads: int(THREADS)
	message: "Using portcullis to prepare: {input}"
	shell: "{params.load} && portcullis prep -o {params.outdir} -l --strandedness={PORTCULLIS_STRAND} -t {threads} {input.ref} {input.bam} > {log} 2>&1"


rule sim_var_portcullis_junc:
	input:
		bai=rules.sim_var_portcullis_prep.output
	output: PORT_DIR+"/junc/portcullis_sim_var.junctions.bed"
	params:
		prepdir=PORT_DIR+"/prep",
		outdir=PORT_DIR+"/junc",
		load=LOAD_PORTCULLIS
	log: PORT_DIR+"/junc.log"
	threads: int(THREADS)
	message: "Using portcullis to analyse potential junctions: {input}"
	shell: "{params.load} && portcullis junc -o {params.outdir} -p portcullis_sim_var -t {threads} {params.prepdir} > {log} 2>&1"


