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

RUN_TIME = config["run_time"]

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
LOAD_CLASS = config["load_class"]
LOAD_PORTCULLIS = config["load_portcullis"]
LOAD_SPANKI = config["load_spanki"]
LOAD_FINESPLICE = config["load_finesplice"]
LOAD_TRUESIGHT = config["load_truesight"]
LOAD_PYTHON3 = config["load_python3"]
LOAD_GT = config["load_gt"]
LOAD_MIKADO = config["load_mikado"]

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
ASSEMBLY_MODES = ["permissive","default","strict"]
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

# Min isoform fraction
PERMISSIVE = 0.01
DEFAULT = 0.1
STRICT = 0.5


#########################
Rules

localrules: all, truesight2bed

# Define 
rule all:
	input: 
		expand(ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.stats", aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/{aln_method}-{reads}-portcullis.bed", aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/tophat-{reads}-spanki.bed", reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/tophat-{reads}-finesplice.bed", reads=INPUT_SETS),
		expand(JUNC_DIR+"/output/truesight-{reads}-truesight.bed", reads=INPUT_SETS),
		expand(ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.bed", asm_method=ASSEMBLY_METHODS, asm_mode=ASSEMBLY_MODES, aln_method=ALIGNMENT_METHODS, reads=INPUT_SETS)
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
		ref=os.path.abspath(REF)
	log: ALIGN_DIR + "/bowtie.index.log"
	threads: 1
	message: "Indexing genome with bowtie"
	shell: "{LOAD_BOWTIE} && bowtie-build {input} {params.out_prefix} > {log} 2>&1 && ln -sf {params.ref} {output.fa_link} && touch -h {output.fa_link}"


rule align_tophat_index:
    input: REF
    output: ALIGN_DIR +"/tophat/index/"+NAME+".4.bt2"
    log: ALIGN_DIR + "/tophat.index.log"
    threads: 1
    message: "Indexing genome with tophat"
    shell: "{LOAD_TOPHAT} && bowtie2-build {REF} {ALIGN_DIR}/tophat/index/{NAME} > {log} 2>&1"




rule align_gsnap_index:
    input: REF
    output: ALIGN_DIR +"/gsnap/index/"+NAME+"/"+NAME+".sachildguide1024"
    log: ALIGN_DIR +"/gsnap.index.log"
    threads: 1
    message: "Indexing genome with gsnap"
    shell: "{LOAD_GMAP} && gmap_build --dir={ALIGN_DIR}/gsnap/index --db={NAME} {input} > {log} 2>&1"



rule align_star_index:
    input: os.path.abspath(REF)
    output: ALIGN_DIR +"/star/index/SAindex"
    params: indexdir=ALIGN_DIR_FULL+"/star/index"
    log: ALIGN_DIR_FULL+"/star.index.log"
    threads: int(THREADS)
    message: "Indexing genome with star"
    shell: "{LOAD_STAR} && cd {ALIGN_DIR_FULL}/star && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input} {INDEX_STAR_EXTRA} > {log} 2>&1 && cd {CWD}"



rule align_hisat_index:
	input: REF
	output: ALIGN_DIR+"/hisat/index/"+NAME+".4.ht2"
	log: ALIGN_DIR+"/hisat.index.log"
	threads: 1
	message: "Indexing genome with hisat"
	shell: "{LOAD_HISAT} && hisat2-build {input} {ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1"





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
    message: "Aligning RNAseq data with tophat"
    run:
        strand = STRANDEDNESS if wildcards.reads.startswith("real") else "fr-unstranded"
        shell("{LOAD_TOPHAT} && {RUN_TIME} tophat2 --output-dir={params.outdir} --num-threads={threads} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} --microexon-search --library-type={strand} {ALIGN_TOPHAT_EXTRA} {params.indexdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {output.bam} {output.link} && touch -h {output.link}")



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
    message: "Aligning RNAseq with gsnap"
    shell: "{LOAD_GMAP} && {LOAD_SAMTOOLS} && {RUN_TIME} gsnap --dir={ALIGN_DIR}/gsnap/index --db={NAME} {ALIGN_GSNAP_EXTRA} --novelsplicing=1 --localsplicedist={MAX_INTRON} --nthreads={threads} --format=sam --npaths=20 {input.r1} {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {output.bam} {output.link} && touch -h {output.link}"




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
    shell: "{LOAD_STAR} && cd {params.outdir} && {RUN_TIME} STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.indexdir} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMin {MIN_INTRON} --alignIntronMax {MAX_INTRON} --alignMatesGapMax 20000 --outFileNamePrefix {params.outdir}/ {ALIGN_STAR_EXTRA} > {log} 2>&1 && cd {CWD} && ln -sf {output.bam} {output.link} && touch -h {output.link}"



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
    run:
        strand = HISAT_STRAND if wildcards.reads.startswith("real") else ""
        shell("{LOAD_HISAT} && {LOAD_SAMTOOLS} && {RUN_TIME} hisat2 -p {threads} --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {strand} -x {params.indexdir} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -b -@ {threads} - > {output.bam} && ln -sf {output.bam} {output.link} && touch -h {output.link}")

rule bam_sort:
	input: 
		bam=ALIGN_DIR+"/output/{aln_method}-{reads}.bam"
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam"
	threads: int(THREADS)
	message: "Using samtools to sort {input.bam}"
	shell: "{LOAD_SAMTOOLS} && samtools sort -o {output} -O bam -m 1G -T sort_{wildcards.aln_method}_{wildcards.reads} -@ {threads} {input.bam}"


rule bam_index:
	input: rules.bam_sort.output
	output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.bai"
	threads: 1
	message: "Using samtools to index: {input}"
	shell: "{LOAD_SAMTOOLS} && samtools index {input}"

rule bam_stats:
        input:
                bam=rules.bam_sort.output,
                idx=rules.bam_index.output
        output: ALIGN_DIR+"/output/{aln_method}-{reads}.sorted.bam.stats"
        params: 
                load=LOAD_SAMTOOLS,
                plot_out=ALIGN_DIR+"/output/plots/{aln_method}-{reads}/{aln_method}-{reads}"
        threads: 1
        message: "Using samtools to collected stats for: {input}"
        shell: "{params.load} && samtools stats {input.bam} > {output} && plot-bamstats -p {params.plot_out} {output}"

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
	run:
		strand = PORTCULLIS_STRAND if wildcards.reads.startswith("real") else "unstranded"
		shell("{params.load} && portcullis prep -o {params.outdir} -l --strandedness={strand} -t {threads} {input.ref} {input.bam} > {log} 2>&1")


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
	run:
		strand = PORTCULLIS_STRAND if wildcards.reads.startswith("real") else "unstranded"		
		shell("{params.load} && {RUN_TIME} portcullis junc -o {params.outdir} -p {wildcards.aln_method}-{wildcards.reads} --strandedness={strand} -t {threads} {params.prepdir} > {log} 2>&1")


rule portcullis_filter:
	input: rules.portcullis_junc.output
	output:
		link=JUNC_DIR+"/output/{aln_method}-{reads}-portcullis.bed",
		unfilt_link=JUNC_DIR+"/output/{aln_method}-{reads}-all.bed",
		tab=PORTCULLIS_DIR+"/{aln_method}-{reads}/filt/{aln_method}-{reads}.pass.junctions.tab",
	params:
		outdir=PORTCULLIS_DIR+"/{aln_method}-{reads}/filt",
		load=LOAD_PORTCULLIS,
		bed=PORTCULLIS_DIR_FULL+"/{aln_method}-{reads}/filt/{aln_method}-{reads}.pass.junctions.bed",
		unfilt_bed=PORTCULLIS_DIR_FULL+"/{aln_method}-{reads}/junc/{aln_method}-{reads}.junctions.bed"
	log: PORTCULLIS_DIR+"/{aln_method}-{reads}-filter.log"
	threads: int(THREADS)
	message: "Using portcullis to filter invalid junctions: {input}"
	shell: "{params.load} && portcullis filter -t {threads} -o {params.outdir} -p {wildcards.aln_method}-{wildcards.reads} {input} > {log} 2>&1 && ln -sf {params.bed} {output.link} && touch -h {output.link} && ln -sf {params.unfilt_bed} {output.unfilt_link} && touch -h {output.unfilt_link}"


rule portcullis_bamfilt:
	input: 
		bam=rules.bam_sort.output,
		tab=rules.portcullis_filter.output.tab
	output:
		bam=PORTCULLIS_DIR+"/{aln_method}-{reads}/bam/{aln_method}-{reads}-portcullis.bam"
	params:
		load=LOAD_PORTCULLIS,
	log: PORTCULLIS_DIR+"/{aln_method}-{reads}-bam.log"
	threads: int(THREADS)
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
    message: "Using SPANKI to analyse junctions: {input.bam}"
    shell: "set +e && {params.load_spanki} && cd {params.outdir} && {RUN_TIME} spankijunc -i {params.bam} -g {params.gtf} -f {params.fa} > {log} 2>&1 && cd {CWD} && if [[ -s {params.all_juncs} ]] ; then {params.load_portcullis} && spanki_filter.py {params.all_juncs} > {params.filt_juncs} && spanki2bed.py {params.filt_juncs} > {output.bed} ; else touch {output.bed} ; fi && ln -sf {params.bed} {output.link} && touch -h {output.link}"

rule spanki_annot:
    input:
        bam=rules.bam_sort.output,
        fa=REF,
        gtf=REF_GTF,
        idx=rules.bam_index.output
    output: link=JUNC_DIR+"/output/{aln_method}-{reads}-spanki_annot.bed",
        bed=JUNC_DIR+"/spanki/{aln_method}-{reads}/junctions_out/{aln_method}-{reads}-spanki_annot.bed"
    params:
        load_spanki=LOAD_SPANKI,
        load_portcullis=LOAD_PORTCULLIS,
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
    shell: "set +e && {params.load_spanki} && cd {params.outdir} && {RUN_TIME} spankijunc -i {params.bam} -g {params.gtf} -f {params.fa} -filter T > {log} 2>&1 && cd {CWD} && if [[ -s {params.all_juncs} ]] ; then {params.load_portcullis} && spanki2bed.py {params.all_juncs} > {output.bed} ; else touch {output.bed} ; fi && ln -sf {params.bed} {output.link} && touch -h {output.link}"
	
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
    message: "Using FineSplice to analyse junctions: {input.bam}"
    shell: "set +e && {params.load_fs} && cd {params.outdir} && if {RUN_TIME} FineSplice.py -i {params.bam} -l {READ_LENGTH} > {log} 2>&1 ; then cd {CWD} && {params.load_portcullis} && fs2bed.py {params.junc} > {output.bed} ; else cd {CWD}; touch {output.bed} ; fi && ln -sf {params.bed} {output.link} && touch -h {output.link}"


rule truesight:
    input:
        idx=rules.align_bowtie_index.output,
        r1=READS_DIR+"/{reads}.R1.fq",
        r2=READS_DIR+"/{reads}.R2.fq"
    output: JUNC_DIR+"/truesight/{reads}/GapAli.junc"
    params:
        load_fs=LOAD_TRUESIGHT,
	index=ALIGN_DIR_FULL + "/bowtie/index/"+NAME,
        outdir=JUNC_DIR+"/truesight/{reads}",
        junc=JUNC_DIR+"/truesight/{reads}/GapAli.junc",
	r1=READS_DIR_FULL+"/{reads}.R1.fq",
	r2=READS_DIR_FULL+"/{reads}.R2.fq"
    log: JUNC_DIR_FULL+"/truesight/{reads}-truesight.log"
    threads: int(THREADS)
    message: "Using Truesight to find junctions"
    shell: "{params.load_fs} && cd {params.outdir} && {RUN_TIME} truesight_pair.pl -i {MIN_INTRON} -I {MAX_INTRON} -v 1 -r {params.index} -p {threads} -o . -f {params.r1} {params.r2} > {log} 2>&1 && cd {CWD}"

rule truesight2bed:
	input: rules.truesight.output
	output:
		link=JUNC_DIR+"/output/truesight-{reads}-truesight.bed",
	        bed=JUNC_DIR+"/truesight/{reads}/truesight-{reads}-truesight.bed"
	params:
		load_portcullis=LOAD_PORTCULLIS,
		bed="../truesight/{reads}/truesight-{reads}-truesight.bed"
	threads: 1
	message: "Creating bed file from truesight output: {input}"
	shell: "{params.load_portcullis} && truesight2bed.py {input} > {output.bed} && ln -sf {params.bed} {output.link} && touch -h {output.link}"


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
		load=LOAD_CUFFLINKS
	log: ASM_DIR+"/cufflinks_{asm_mode}-{aln_method}-{reads}.log"
	threads: int(THREADS)
	message: "Using cufflinks to assemble {input.bam}"
	run:
		mode = PERMISSIVE if wildcards.asm_mode == "permissive" else STRICT if wildcards.asm_mode == "strict" else DEFAULT
		shell("{params.load} && cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} -F {mode} --no-update-check {input.bam} > {log} 2>&1 && ln -sf {params.link_src} {output.gtf} && touch -h {output.gtf}")



rule asm_stringtie:
	input: 	bam=rules.bam_sort.output
	output: 
		link=ASM_DIR+"/output/stringtie_{asm_mode}-{aln_method}-{reads}.gtf",
		gtf=ASM_DIR+"/stringtie_{asm_mode}-{aln_method}-{reads}/stringtie_{asm_mode}-{aln_method}-{reads}.gtf"
	params:
		load=LOAD_STRINGTIE,
		gtf=ASM_DIR+"/stringtie_{asm_mode}-{aln_method}-{reads}/stringtie_{asm_mode}-{aln_method}-{reads}.gtf",
		link_src="../stringtie_{asm_mode}-{aln_method}-{reads}/stringtie_{asm_mode}-{aln_method}-{reads}.gtf",
		name="Stringtie_{asm_mode}_{aln_method}_{reads}"
	log: ASM_DIR+"/stringtie_{asm_mode}-{aln_method}-{reads}.log"
	threads: int(THREADS)
	message: "Using stringtie to assemble: {input.bam}"
	run:
		mode = PERMISSIVE if wildcards.asm_mode == "permissive" else STRICT if wildcards.asm_mode == "strict" else DEFAULT
		shell("{params.load} && {RUN_TIME} stringtie {input.bam} -l {params.name} -f {mode} -m 200 -o {params.gtf} -p {threads} > {log} 2>&1 && ln -sf {params.link_src} {output.link} && touch -h {output.link}")


rule asm_class:
	input:
		bam=rules.bam_sort.output,
		ref=REF
	output:
		link=ASM_DIR+"/output/class_{asm_mode}-{aln_method}-{reads}.gtf",
		gtf=ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}/class_{asm_mode}-{aln_method}-{reads}.gtf"
	params:
		outdir=ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}",
		load=LOAD_CLASS,
		gtf=ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}/class_{asm_mode}-{aln_method}-{reads}.gtf",
		link_src="../class_{asm_mode}-{aln_method}-{reads}/class_{asm_mode}-{aln_method}-{reads}.gtf"
	log: ASM_DIR+"/class_{asm_mode}-{aln_method}-{reads}.log"
	threads: int(THREADS)
	message: "Using class to assemble: {input.bam}"
	run:
		mode = PERMISSIVE if wildcards.asm_mode == "permissive" else STRICT if wildcards.asm_mode == "strict" else DEFAULT
		shell("{params.load} && class_run.py -c '-F {mode}' -p {threads} {input.bam} > {output.gtf} 2> {log} && ln -sf {params.link_src} {output.link} && touch -h {output.link}")


rule gtf_2_bed:
	input:
		gtf=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.gtf"
	output:
		bed=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.bed"
	params:
		load_gt=LOAD_GT,
		load_p=LOAD_PORTCULLIS,
		gff=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.gff3",
		gffi=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.introns.gff3"
	message: "Converting GTF to BED for: {input.gtf}"
	shell: "{params.load_gt} && gt gtf_to_gff3 -tidy -force -o {params.gff} {input.gtf} 2> /dev/null && gt gff3 -sort -tidy -addintrons -force -o {params.gffi} {params.gff} && {params.load_p} && gff2bed.py {params.gffi} > {output.bed} && rm {params.gff} {params.gffi}"

rule gtf_stats:
	input:
		gtf=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.gtf",
		ref=REF_GTF
	output: 
		comp=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.comp.stats",
		stats=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.stats"
	params: 
		load_mikado=LOAD_MIKADO,
		prefix=ASM_DIR+"/output/{asm_method}_{asm_mode}-{aln_method}-{reads}.comp"
	message: "Calculating stats for: {input.gtf}"
	shell: "{params.load_mikado} && mikado.py util stats {input.gtf} > {output.stats} && mikado.py compare -r {input.ref} -p {input.gtf} -o {params.prefix}"


