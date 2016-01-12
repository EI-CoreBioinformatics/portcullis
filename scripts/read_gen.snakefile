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
THREADS = config["threads"]
READ_LENGTH = config["min_read_length"]
READ_LENGTH_MINUS_1 = int(READ_LENGTH) - 1

LOAD_BOWTIE = config["load_bowtie"]
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

# Define 
rule all:
	input: 
		SIM_DIR_FULL+"/var/sim.sam",
		SIM_DIR_FULL+"/fixed/sim.sam"
		#expand(PORT_DIR + "/output/portcullis-{align_method}-{reads}.unfiltered.bed", align_method=ALIGNMENT_METHODS, reads=INPUT_SETS)
		

#rule clean:
#	shell: "rm -rf {out}"


rule trim_reads:
	input: 
		r1=R1,
		r2=R2
	output:
		linkr1=READS_DIR+"/real.R1.fq",
		linkr2=READS_DIR+"/real.R2.fq"
	params:
		outdir=READS_DIR+"/trim",
		r1=READS_DIR_FULL+"/trim/r1_val_1.fq",
		r2=READS_DIR_FULL+"/trim/r2_val_2.fq",
	log: READS_DIR+"/trim/trim.log"		
	shell: "{LOAD_TRIMGALORE}; trim_galore --paired --length {READ_LENGTH} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"



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
	shell: "{LOAD_BOWTIE}; bowtie-build {input} {params.outdir}"


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
	shell: "{LOAD_BOWTIE}; bowtie -p {threads} --chunkmbs 1000 -X 500 {params.indexdir} -1 {input.r1} -2 {input.r2} > {output.map} 2> {log}"
	


rule align_tophat_index:
        input: REF 
        output: ALIGN_DIR +"/tophat/index/"+NAME+".4.bt2"
        log: ALIGN_DIR + "/tophat.index.log"
        threads: 1
        message: "Indexing genome with tophat"
        shell: "{LOAD_TOPHAT}; bowtie2-build {REF} {ALIGN_DIR}/tophat/index/{NAME} > {log} 2>&1"




rule align_gsnap_index:
        input: REF
        output: ALIGN_DIR +"/gsnap/index/"+NAME+"/"+NAME+".sachildguide1024"
        log: ALIGN_DIR +"/gsnap.index.log"
        threads: 1
        message: "Indexing genome with gsnap"
        shell: "{LOAD_GMAP}; gmap_build --dir={ALIGN_DIR}/gsnap/index --db={NAME} {input} > {log} 2>&1"



rule align_star_index:
        input: os.path.abspath(REF)
        output: ALIGN_DIR +"/star/index/SAindex"
        params: indexdir=ALIGN_DIR_FULL+"/star/index"
        log: ALIGN_DIR_FULL+"/star.index.log"
        threads: int(THREADS)
        message: "Indexing genome with star"
        shell: "{LOAD_STAR}; cd {ALIGN_DIR_FULL}/star; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input} {INDEX_STAR_EXTRA} > {log} 2>&1; cd {CWD}"

rule align_star_index_wref:
        input:
                fa=os.path.abspath(REF),
                gtf=os.path.abspath(REF_GTF)	
        output: ALIGN_DIR +"/star/index_wref/SAindex"
        params: indexdir=ALIGN_DIR_FULL+"/star/index_wref"
        log: ALIGN_DIR_FULL+"/star.index_wref.log"
        threads: int(THREADS)
        message: "Indexing genome with star using GTF reference"
        shell: "{LOAD_STAR}; cd {ALIGN_DIR_FULL}/star; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang {READ_LENGTH_MINUS_1} {INDEX_STAR_EXTRA} > {log} 2>&1; cd {CWD}"



rule align_hisat_index:
	input: REF
	output: ALIGN_DIR+"/hisat/index/"+NAME+".4.ht2"
	log: ALIGN_DIR+"/hisat.index.log"
	threads: 1
	message: "Indexing genome with hisat"
	shell: "{LOAD_HISAT}; hisat2-build {input} {ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1"



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
		bam=ALIGN_DIR_FULL+"/tophat_wref/accepted_hits.bam",
	log: ALIGN_DIR + "/tophat-real_wref.log"
	threads: int(THREADS)
	message: "Aligning RNAseq data with tophat using GTF reference: {input.r1} {input.r2} {input.gtf}"
	shell: "{LOAD_TOPHAT}; tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={STRANDEDNESS} --GTF={input.gtf} {ALIGN_TOPHAT_EXTRA} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {params.bam} {output.bam} && touch -h {output.bam}"


rule align_star_wref:
        input:
                r1=os.path.abspath(rules.trim_reads.output.linkr1),
		r2=os.path.abspath(rules.trim_reads.output.linkr2),
                index=rules.align_star_index_wref.output
        output:
                bam=ALIGN_DIR+"/output/star-real_wref.bam",
		bed=ALIGN_DIR+"/star_wref/star_wref.bed"
        params:
                outdir=ALIGN_DIR_FULL+"/star_wref",
                indexdir=ALIGN_DIR_FULL+"/star/index_wref",
		tab=ALIGN_DIR+"/star_wref/SJ.out.tab",
		bam=ALIGN_DIR_FULL+"/star_wref/Aligned.out.bam"
        log: ALIGN_DIR_FULL+"/star-real_wref.log"
        threads: int(THREADS)
        message: "Aligning input with star using GTF indexed genome"
        shell: "{LOAD_STAR}; cd {params.outdir}; STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.indexdir} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON} --alignMatesGapMax 20000 --outFileNamePrefix {params.outdir}/ {ALIGN_STAR_EXTRA} > {log} 2>&1 && cd {CWD} && ln -sf {params.bam} {output.bam} && touch -h {output.bam} && {LOAD_PORTCULLIS} && {LOAD_PYTHON3} && star_tab2bed {params.tab} > {output.bed}"


rule bam_sort_wref:
	input: ALIGN_DIR+"/output/{align_ref_method}-real_wref.bam"
	output: ALIGN_DIR+"/output/{align_ref_method}-real_wref.sorted.bam"
	threads: int(THREADS)
	message: "Using samtools to sort {input}"
	shell: "{LOAD_SAMTOOLS}; samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.align_ref_method}_wref -@ {threads} {input}"

rule bam_index_wref:
	input: rules.bam_sort_wref.output
	output: ALIGN_DIR+"/output/{align_ref_method}-real_wref.sorted.bam.bai"
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{LOAD_SAMTOOLS}; samtools index {input}"

rule asm_cufflinks_wref:
	input:
		bam=ALIGN_DIR+"/output/tophat-real_wref.sorted.bam",
		ref=REF_GTF
	output:
		isoforms=ASM_DIR+"/cufflinks-tophat_wref-real/isoforms.fpkm_tracking"
	params: 
		outdir=ASM_DIR+"/cufflinks-tophat_wref-real",
	log: ASM_DIR+"/cufflinks-tophat_wref-real.log"
	threads: int(THREADS)
	message: "Using cufflinks to assemble: {input.bam}"
	shell: "{LOAD_CUFFLINKS}; cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --GTF={input.ref}  --no-update-check {ASM_CUFFLINKS_EXTRA} {input.bam} > {log} 2>&1"

rule convert_isoform_fpkm_2_spanki:
	input: isoforms=ASM_DIR+"/cufflinks-tophat_wref-real/isoforms.fpkm_tracking"
	output: cov=ASM_DIR+"/cufflinks-tophat_wref-real/transcripts.cov"
	threads: 1
	message: "Converting cufflinks isoforms file to spanki format"
	shell: "{LOAD_PORTCULLIS} && {LOAD_PYTHON3} && cufflinks2spankicov.py {input.isoforms} > {output.cov}"


rule simgen_model:
	input:
		map=rules.simgen_bowtie_align.output.map
	output: SIM_DIR+"/model/logfile.txt"
	params: outdir=SIM_DIR+"/model"
	log: SIM_DIR+"/spanki_model.log"
	threads: 1
	shell: "{LOAD_SPANKI}; spankisim_models -i {input.map} -e 2 -l {READ_LENGTH} -o {params.outdir} > {log} 2>&1"


rule sim_fixed_reads:
	input:
		gtf=REF_GTF,
                fa=REF,
		model=rules.simgen_model.output
	output: 
		sam=SIM_DIR_FULL+"/fixed/sim.sam",
		linkr1=READS_DIR+"/sim_fixed.R1.fq",
		linkr2=READS_DIR+"/sim_fixed.R2.fq"
	params: 
		outdir=SIM_DIR+"/fixed",
		mdir=SIM_DIR+"/model",
		r1=SIM_DIR_FULL+"/fixed/sim_1.fastq",
		r2=SIM_DIR_FULL+"/fixed/sim_2.fastq",
	log: SIM_DIR+"/spanki_readgen_fixed.log"
	threads: 1	
	shell: "{LOAD_SPANKI}; {LOAD_CUFFLINKS}; spankisim_transcripts -g {input.gtf} -f {input.fa} -o {params.outdir} -m custom -mdir {params.mdir} -cov 30 -bp {READ_LENGTH} -ends 2 -frag 250 > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"
	


rule sim_var_reads:
	input:
		gtf=REF_GTF,
                fa=REF,
                transcript_cov=rules.convert_isoform_fpkm_2_spanki.output.cov,
		model=rules.simgen_model.output
	output:
		sam=SIM_DIR_FULL+"/var/sim.sam",
		linkr1=READS_DIR+"/sim_var.R1.fq",
		linkr2=READS_DIR+"/sim_var.R2.fq"
	params: 
		outdir=SIM_DIR+"/var",
		mdir=SIM_DIR+"/model",
		r1=SIM_DIR_FULL+"/var/sim_1.fastq",
		r2=SIM_DIR_FULL+"/var/sim_2.fastq",
	log: SIM_DIR+"/spanki_readgen_var.log"
	threads: 1	
	shell: "{LOAD_SPANKI}; {LOAD_CUFFLINKS}; spankisim_transcripts -g {input.gtf} -f {input.fa} -o {params.outdir} -m custom -mdir {params.mdir} -t {input.transcript_cov} -bp {READ_LENGTH} -ends 2 -frag 250 > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"


'''

rule align_tophat:
        input:
                r1=READS_DIR+"/{reads}.R1.fq",
                r2=READS_DIR+"/{reads}.R2.fq",
                index=rules.align_tophat_index.output
        output:
                bam=ALIGN_DIR_FULL+"/tophat/{reads}/accepted_hits.bam",
                link=ALIGN_DIR+"/output/tophat-{reads}.bam"
        params: 
                outdir=ALIGN_DIR+"/tophat/{reads}",
                indexdir=ALIGN_DIR+"/tophat/index/"+NAME
        log: ALIGN_DIR + "/tophat-{reads}.log"
        threads: int(THREADS)
        message: "Aligning RNAseq data with tophat: {input.r1}; {input.r2}"
        shell: "{LOAD_TOPHAT}; tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={STRANDEDNESS} {ALIGN_TOPHAT_EXTRA} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {output.bam} {output.link} && touch -h {output.link}"


rule align_gsnap:
        input:
                r1=READS_DIR+"/{reads}.R1.fq",
                r2=READS_DIR+"/{reads}.R2.fq",
                index=rules.align_gsnap_index.output
        output:
                bam=ALIGN_DIR_FULL+"/gsnap/{reads}/gsnap.bam",
                link=ALIGN_DIR+"/output/gsnap-{reads}.bam"
        log: ALIGN_DIR+"/gsnap-{reads}.log"
        threads: int(THREADS)
        message: "Aligning RNAseq with gsnap: {input.r1}; {input.r2}"
        shell: "{LOAD_GMAP}; {LOAD_SAMTOOLS}; gsnap --dir={ALIGN_DIR}/gsnap/index --db={NAME} {ALIGN_GSNAP_EXTRA} --novelsplicing=1 --localsplicedist={MAX_INTRON} --nthreads={threads} --format=sam --npaths=20 {input.r1} {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {output.bam} {output.link} && touch -h {output.link}"




rule align_star:
        input:
                r1=READS_DIR_FULL+"/{reads}.R1.fq",
                r2=READS_DIR_FULL+"/{reads}.R2.fq",
                index=rules.align_star_index.output
        output:
                bam=ALIGN_DIR_FULL+"/star/{reads}/Aligned.out.bam",
                link=ALIGN_DIR+"/output/star-{reads}.bam"
        params:
                outdir=ALIGN_DIR_FULL+"/star/{reads}",
                indexdir=ALIGN_DIR_FULL+"/star/index"
        log: ALIGN_DIR_FULL+"/star-{reads}.log"
        threads: int(THREADS)
        message: "Aligning input with star"
        shell: "{LOAD_STAR}; cd {params.outdir}; STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.indexdir} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON} --alignMatesGapMax 20000 --outFileNamePrefix {params.outdir}/ {ALIGN_STAR_EXTRA} > {log} 2>&1; cd {CWD}; ln -sf {output.bam} {output.link} && touch -h {output.link}"



rule align_hisat:
	input: 
		r1=READS_DIR+"/{reads}.R1.fq",
		r2=READS_DIR+"/{reads}.R2.fq",
		index=rules.align_hisat_index.output
	output:
		bam=ALIGN_DIR_FULL+"/hisat/{reads}/hisat.bam",
		link=ALIGN_DIR+"/output/hisat-{reads}.bam"
	params:
		indexdir=ALIGN_DIR+"/hisat/index/"+NAME
	log: ALIGN_DIR+"/hisat-{reads}.log"
        threads: int(THREADS)
        message: "Aligning input with hisat"
	shell: "{LOAD_HISAT}; {LOAD_SAMTOOLS}; hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {HISAT_STRAND} -x {params.indexdir} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam}; ln -sf {output.bam} {output.link} && touch -h {output.link}"



rule bam_sort:
	input: ALIGN_DIR+"/output/{align_method}-{reads}.bam"
	output: ALIGN_DIR+"/output/{align_method}-{reads}.sorted.bam"
	threads: int(THREADS)
	message: "Using samtools to sort {input}"
	shell: "{LOAD_SAMTOOLS}; samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.align_method}_{wildcards.reads} -@ {threads} {input}"


rule bam_index:
	input: rules.bam_sort.output
	output: ALIGN_DIR+"/output/{align_method}-{reads}.sorted.bam.bai"
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{LOAD_SAMTOOLS}; samtools index {input}"



rule asm_cufflinks:
        input:
                bam=rules.bam_sort.output,
                ref=REF
        output:
                gtf=ASM_DIR_FULL+"/cufflinks-{align_method}-real_wref/transcripts.gtf",
                link=ASM_DIR+"/output/cufflinks-{align_method}-real_wref.gtf"
        params: outdir=ASM_DIR+"/cufflinks-{align_method}-real_wref"
        log: ASM_DIR+"/cufflinks-{align_method}-real_wref.log"
        threads: int(THREADS)
        message: "Using cufflinks to assemble: {input.bam}"
        shell: "{LOAD_CUFFLINKS}; cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --no-update-check {ASM_CUFFLINKS_EXTRA} {input.bam} > {log} 2>&1 && ln -sf {output.gtf} {output.link} && touch -h {output.link}"




rule asm_stringtie:
        input: 
                bam=rules.bam_sort.output,
                ref={REF}
        output: 
                gtf=ASM_DIR_FULL+"/stringtie-{align_method}-{reads}/stringtie-{align_method}-{reads}.gtf",
                link=ASM_DIR+"/output/stringtie-{align_method}-{reads}.gtf"
        log: ASM_DIR+"/stringtie-{align_method}-{reads}.log"
        threads: int(THREADS)
        message: "Using stringtie to assemble: {input.bam}"
        shell: "{LOAD_STRINGTIE}; stringtie {input.bam} -l Stringtie_{wildcards.align_method}_{wildcards.reads} -f 0.05 -m 200 {ASM_STRINGTIE_EXTRA} -o {output.gtf} -p {threads} > {log} 2>&1 && ln -sf {output.gtf} {output.link} && touch -h {output.link}"





rule portcullis_prep:
        input:
                ref={REF},
                bam=rules.bam_sort.output
        output: PORT_DIR+"/{align_method}-{reads}-prep/portcullis.sorted.alignments.bam.bai"
        params: outdir=PORT_DIR+"/{align_method}-{reads}-prep"
        log: PORT_DIR+"/{align_method}-{reads}-prep.log"
        threads: int(THREADS)
        message: "Using portcullis to prepare: {input}"
        shell: "{LOAD_PORTCULLIS}; portcullis prep -o {params.outdir} -l -s {PORTCULLIS_STRAND} -t {threads} {input.ref} {input.bam} > {log} 2>&1"


rule portcullis_junc:
        input:
                bai=rules.portcullis_prep.output
	output:
		bed=PORT_DIR_FULL+"/{align_method}-{reads}-junc/{align_method}-{reads}.junctions.bed",
		link=PORT_DIR+"/output/portcullis-{align_method}-{reads}.unfiltered.bed"
	params:
		prepdir=PORT_DIR+"/{align_method}-{reads}-prep",
		outdir=PORT_DIR+"/{align_method}-{reads}-junc"
	log: PORT_DIR+"/{align_method}-{reads}-junc.log"
	threads: int(THREADS)
	message: "Using portcullis to analyse potential junctions: {input}"
	shell: "{LOAD_PORTCULLIS}; portcullis junc -o {params.outdir} -p {wildcards.align_method}-{wildcards.reads} -t {threads} {params.prepdir} > {log} 2>&1 && ln -sf {output.bed} {output.link} && touch -h {output.link}"

'''
