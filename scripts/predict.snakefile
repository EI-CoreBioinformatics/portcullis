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
LOAD_SPANKI = config["load_spanki"]
LOAD_FINESPLICE = config["load_finesplice"]
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
SIM_INPUT_SET = ["sim_fixed","sim_var"]
REAL_INPUT_SET = ["real"]
ALIGNMENT_METHODS = config["align_methods"]
ASSEMBLY_METHODS = config["asm_methods"]
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

CWD = os.getcwd()


HISAT_STRAND = "--rna-strandness=RF" if STRANDEDNESS == "fr-firststrand" else "--rna-strandness=FR" if STRANDEDNESS == "fr-secondstrand" else ""
PORTCULLIS_STRAND = "firststrand" if STRANDEDNESS == "fr-firststrand" else "secondstrand" if STRANDEDNESS == "fr-secondstrand" else "unstranded"




#########################
Rules

# Define 
rule all:
	input: 
		expand(JUNC_DIR+"/output/{aln_method}-{reads}-{junc_method}.bed", aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS, junc_method=JUNC_METHODS)
		

#rule clean:
#	shell: "rm -rf {out}"

	


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



rule align_hisat_index:
	input: REF
	output: ALIGN_DIR+"/hisat/index/"+NAME+".4.ht2"
	log: ALIGN_DIR+"/hisat.index.log"
	threads: 1
	message: "Indexing genome with hisat"
	shell: "{LOAD_HISAT}; hisat2-build {input} {ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1"





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
        outdir=ALIGN_DIR+"/tophat/{reads}",
        indexdir=ALIGN_DIR+"/tophat/index/"+NAME
    log: ALIGN_DIR + "/tophat-{reads}.log"
    threads: int(THREADS)
    message: "Aligning RNAseq data with tophat: {input.r1}; {input.r2}"
    run:
        strand = STRANDEDNESS if wildcards.reads.startswith("real") else "fr-unstranded"
        shell("{LOAD_TOPHAT}; tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={strand} {ALIGN_TOPHAT_EXTRA} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {output.bam} {output.link} && touch -h {output.link}")



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
    message: "Aligning input with star: {input.r1}; {input.r2}"
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
    message: "Aligning input with hisat: {input.r1}; {input.r2}"
    run:
        strand = HISAT_STRAND if wildcards.reads.startswith("real") else ""
        shell("{LOAD_HISAT}; {LOAD_SAMTOOLS}; hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {strand} -x {params.indexdir} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam}; ln -sf {output.bam} {output.link} && touch -h {output.link}")

rule bam_sort:
	input: 
		bam=ALIGN_DIR+"/output/{aln_method}-{reads}.bam"
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam"
	threads: int(THREADS)
	message: "Using samtools to sort {input.bam}"
	shell: "{LOAD_SAMTOOLS}; samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.aln_method}_{wildcards.reads} -@ {threads} {input.bam}"


rule bam_index:
	input: rules.bam_sort.output
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.bai"
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{LOAD_SAMTOOLS}; samtools index {input}"



rule asm_cufflinks:
    input:
        bam=rules.bam_sort.output,
        ref=REF
    output:
        gtf=ASM_DIR_FULL+"/cufflinks-{aln_method}-real_wref/transcripts.gtf",
        link=ASM_DIR+"/output/cufflinks-{aln_method}-real_wref.gtf"
    params: outdir=ASM_DIR+"/cufflinks-{aln_method}-real_wref"
    log: ASM_DIR+"/cufflinks-{aln_method}-real_wref.log"
    threads: int(THREADS)
    message: "Using cufflinks to assemble: {input.bam}"
    shell: "{LOAD_CUFFLINKS}; cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --no-update-check {ASM_CUFFLINKS_EXTRA} {input.bam} > {log} 2>&1 && ln -sf {output.gtf} {output.link} && touch -h {output.link}"


rule asm_stringtie:
    input:
        bam=rules.bam_sort.output,
        ref={REF}
    output:
        gtf=ASM_DIR_FULL+"/stringtie-{aln_method}-{reads}/stringtie-{aln_method}-{reads}.gtf",
        link=ASM_DIR+"/output/stringtie-{aln_method}-{reads}.gtf"
    log: ASM_DIR+"/stringtie-{aln_method}-{reads}.log"
    threads: int(THREADS)
    message: "Using stringtie to assemble: {input.bam}"
    shell: "{LOAD_STRINGTIE}; stringtie {input.bam} -l Stringtie_{wildcards.aln_method}_{wildcards.reads} -f 0.05 -m 200 {ASM_STRINGTIE_EXTRA} -o {output.gtf} -p {threads} > {log} 2>&1 && ln -sf {output.gtf} {output.link} && touch -h {output.link}"

rule portcullis_prep:
    input:
        ref=REF,
        bam=rules.bam_sort.output,
        idx=rules.bam_index.output
    output: PORTCULLIS_DIR+"/{aln_method}-{reads}/prep/portcullis.sorted.alignments.bam.bai"
    params:
        outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/prep",
        load=LOAD_PORTCULLIS
    log: PORTCULLIS_DIR+"/{aln_method}-{reads}-prep.log"
    threads: int(THREADS)
    message: "Using portcullis to prepare: {input}"
    shell: "{params.load}; portcullis prep -o {params.outdir} -l -s {PORTCULLIS_STRAND} -t {threads} {input.ref} {input.bam} > {log} 2>&1"


rule portcullis_junc:
    input:
        bai=rules.portcullis_prep.output
    output: PORTCULLIS_DIR+"/{aln_method}-{reads}/junc/{aln_method}-{reads}.junctions.tab"
    params:
        prepdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/prep",
        outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/junc",
        load=LOAD_PORTCULLIS
    log: PORTCULLIS_DIR+"/{aln_method}-{reads}-junc.log"
    threads: int(THREADS)
    message: "Using portcullis to analyse potential junctions: {input}"
    shell: "{params.load}; portcullis junc -o {params.outdir} -p {wildcards.aln_method}-{wildcards.reads} -t {threads} {params.prepdir} > {log} 2>&1"


rule portcullis_filter:
    input: rules.portcullis_junc.output
    output:
        link=JUNC_DIR+"/output/{aln_method}-{reads}-portcullis.bed",
        unfilt_link=JUNC_DIR+"/output/{aln_method}-{reads}-all.bed"
    params:
        outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/filt",
        load=LOAD_PORTCULLIS,
        bed=PORTCULLIS_DIR_FULL+"/{aln_method}-{reads}/filt/{aln_method}-{reads}.pass.junctions.bed",
        unfilt_bed=PORTCULLIS_DIR_FULL+"/{aln_method}-{reads}/junc/{aln_method}-{reads}.junctions.bed"
    log: PORTCULLIS_DIR+"/{aln_method}-{reads}-filter.log"
    threads: 1
    message: "Using portcullis to filter invalid junctions: {input}"
    shell: "{params.load}; portcullis filter -o {params.outdir} -p {wildcards.aln_method}-{wildcards.reads} {input} > {log} 2>&1 && ln -sf {params.bed} {output.link} && touch -h {output.link} && ln -sf {params.unfilt_bed} {output.unfilt_link} && touch -h {output.unfilt_link}"


rule spanki:
    input:
        bam=rules.bam_sort.output,
        fa=REF,
        gtf=REF_GTF,
        idx=rules.bam_index.output
    output: link=JUNC_DIR+"/output/{aln_method}-{reads}-spanki.bed",
        bed=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki.bed"
    params:
        load_spanki=LOAD_SPANKI,
        load_portcullis=LOAD_PORTCULLIS,
        outdir=JUNC_DIR_FULL+"/spanki/{aln_method}-{reads}",
        bam=ALIGN_DIR_FULL+"/output/{aln_method}-{reads}.sorted.bam",
        fa=os.path.abspath(REF),
        gtf=os.path.abspath(REF_GTF),
        all_juncs=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/juncs.all",
        filt_juncs=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/juncs.filtered",
        bed=JUNC_DIR_FULL+"/spanki/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki.bed"
    log: JUNC_DIR_FULL+"/spanki/{aln_method}-{reads}-spanki.log"
    threads: 1
    message: "Using SPANKI to anlayse junctions: {input.bam}"
    shell: "{params.load_spanki} ; cd {params.outdir} ; spankijunc -i {params.bam} -g {params.gtf} -f {params.fa} > {log} 2>&1 ; cd {CWD} ; if [[ -s {params.all_juncs} ]] ; then {params.load_portcullis} && spanki_filter.py {params.all_juncs} > {params.filt_juncs} && spanki2bed.py {params.filt_juncs} > {output.bed} ; else touch {output.bed} ; fi ; ln -sf {params.bed} {output.link} && touch -h {output.link}"

	
rule finesplice:
    input:
        bam=rules.bam_sort.output,
        idx=rules.bam_sort.output
    output: link=JUNC_DIR+"/output/{aln_method}-{reads}-finesplice.bed",
        bed=JUNC_DIR+"/finesplice/{aln_method}-{reads}/{aln_method}-{reads}-finesplice.bed"
    params:
        load_fs=LOAD_FINESPLICE,
        load_portcullis=LOAD_PORTCULLIS,
        bam=ALIGN_DIR_FULL+"/output/{aln_method}-{reads}.sorted.bam",
        outdir=JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}",
        junc=JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}/{aln_method}-{reads}.sorted.accepted.junc",
        bed=JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}/{aln_method}-{reads}-finesplice.bed"
    log: JUNC_DIR_FULL+"/finesplice/{aln_method}-{reads}-finesplice.log"
    threads: 1
    message: "Using FineSplice to anlayse junctions: {input.bam}"
    shell: "{params.load_fs} && cd {params.outdir} ; if FineSplice.py -i {params.bam} -l {READ_LENGTH} > {log} 2>&1 ; then cd {CWD}; {params.load_portcullis} && fs2bed.py {params.junc} > {output.bed} ; else cd {CWD}; touch {output.bed} ; fi ; ln -sf {params.bed} {output.link} && touch -h {output.link}"


