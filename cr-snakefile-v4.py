import os
import glob

shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; echo END at $(date)")

configfile: "config.json"

localrules: all

samples=config["samples"]

wildcard_constraints:
    #cutoff="^[+-]?((\d+(\.\d+)?)|(\.\d+))$",
    #qv="^a-zA-z+",
    merge="\d+",
    t="\w+"
    

#rule 1: fastqc
fastq_file_dir=config["fastq_file_dir"]
fastq_prefix=config["fastq_prefix"]
fastq_suffix=config["fastq_suffix"]
fastq_mid=config["fastq_mid"]

#rule 1.5: cutadapt
cutadapt_fwd=config["cutadapt_a"]
cutadapt_rev=config["cutadapt_A"]

#rule 2: align
assembly=config["assembly"]
bowtie2_params=config["bowtie2_params"]

#rule 3a,3b: qualtrim and nodup
mapq=config["qualtrim_mapq"]
picard_params=config["picard_params"]
samtools_nodup_params=config["samtools_nodup_params"]

#rule 4: process bams
genomeSizeFile=config["genomeSizeFile"]
bamCoverage_params=config["bamCoverage_params"]

#rule 5a, 5b: seacr, macs2
seacr_path=config["seacr_params"]["seacr_path"]
seacr_params_type=config["seacr_params"]["type"]
seacr_params_cutoff=config["seacr_params"]["cutoff"]
seacr_params_merge=config["seacr_params"]["merge"]
macs2_q=config["macs2_params"]["q"]
macs2_genome=config["macs2_params"]["genome"]
macs2_type=config["macs2_params"]["type"]
bOn=""
if (macs2_type=="broad"):
    bOn="--broad" 

#rule 6: count
macs2_merge=config["macs2_params"]["merge"]

#rule 7abc: deeptools solo plots
list_of_samples='-'.join(map(str,samples))
bamPE_params=config["deeptools"]["bamPE_params"]
plotCov_params=config["deeptools"]["plotCov_params"]
blacklist_file=config["deeptools"]["blacklist_file"]
mBS_params=config["deeptools"]["multiBamSum_params"]
plotCorr_params=config["deeptools"]["plotCorr_params"]

#rule 8ab: deeptools compute Matrix and plotHeatmap/profile macs2 and seacr
cM_params=config["deeptools"]["computeMatrix_params"]
plotHeatmap_params=config["deeptools"]["plotHeatmap_params"]





###for testing###
##lols=["lol","lolol"]
##fs=["f","ff"]
##rule temp:
##    input:
##        bdg="SM/bdg/{sample}.NoDup.sorted.bedpe.bdg"
##    output:
##        macs="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_control_lambda.txt"
##    params:
##        p=seacr_path,
##        q="{qv}"
##    shell:
##        "echo {params.q} > {output.macs}"
###end testing###

rule all:
    input:
        expand("SM/analysis_tars/{sample}-analysis.tar.gz",sample=samples),
    run:
        print(list_of_samples)



rule make_archive:
    input:
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+"_fastqc.html",
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+"_fastqc.zip",
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R2"+fastq_suffix+"_fastqc.html",
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R2"+fastq_suffix+"_fastqc.zip",
	"SM/cutadapt/{sample}.cutadapt.R1.fastq.gz",
	"SM/cutadapt/{sample}.cutadapt.R2.fastq.gz",
        "SM/sambam/{sample}.metrics",
        "SM/sambam/{sample}.bam",
        "SM/qualtrim_nodup/{sample}.QualTrimmed.bam",
        "SM/qualtrim_nodup/{sample}.QualTrimmed.Picard.bam",
        "SM/qualtrim_nodup/{sample}.metrix",
        "SM/qualtrim_nodup/{sample}.NoDup.bam",
        "SM/qualtrim_nodup/{sample}.NoDup.sorted.bam",
        "SM/qualtrim_nodup/{sample}.NoDup.sorted.bam.bai",
        "SM/bed/{sample}.NoDup.bed",
        "SM/bed/{sample}.NoDup.sorted.bed",
        "SM/bed/{sample}.NoDup.sorted.clipped.bed",
        "SM/bdg/{sample}.NoDup.sorted.clipped.bdg",
        "SM/bw/{sample}.NoDup.sorted.clipped.bw",
        "SM/bedpe/{sample}.NoDup.sorted.bedpe",
        "SM/bw/{sample}.NoDup.SeqDepthNorm.bw",
        "SM/bdg/{sample}.NoDup.sorted.bedpe.bdg",
        expand("SM/seacr/{sample}_SEACR_{t}-{cutoff}.{t}.bed",t=seacr_params_type,cutoff=seacr_params_cutoff,allow_missing=True),
        expand("SM/seacr/merged/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.bed",t=seacr_params_type,cutoff=seacr_params_cutoff,merge=seacr_params_merge,allow_missing=True),
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_control_lambda.bdg",qv=macs2_q,allow_missing=True),
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_control_lambda.bw",qv=macs2_q,allow_missing=True),
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_peaks."+macs2_type+"Peak",qv=macs2_q,allow_missing=True),
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_peaks.xls",qv=macs2_q,allow_missing=True),
#        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_summits.bed",qv=macs2_q,allow_missing=True),# comment out for broad Peaks!
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_treat_pileup.bdg",qv=macs2_q,allow_missing=True),
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_treat_pileup.bw",qv=macs2_q,allow_missing=True),
        expand("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/merged_"+macs2_merge+"bp_.sorted.macs2.bed",qv=macs2_q,allow_missing=True),
        expand("SM/peakcounts/{sample}.q{qv}."+macs2_type+"/merged_"+macs2_merge+"bp_.macs2.counts",qv=macs2_q,allow_missing=True),
        expand("SM/peakcounts/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.counts",t=seacr_params_type,cutoff=seacr_params_cutoff,merge=seacr_params_merge,allow_missing=True),
        "SM/deeptools/coverage.{sample}.pdf",
        "SM/deeptools/fragmentSize.{sample}.pdf",
        "SM/deeptools/fragmentSize.{sample}.log.pdf",
        "SM/deeptools/coverage."+list_of_samples+".pdf",
        "SM/deeptools/fragmentSize."+list_of_samples+".pdf",
        "SM/deeptools/fragmentSize."+list_of_samples+".log.pdf",
        "SM/deeptools/"+list_of_samples+"_mBS.npz",
        "SM/deeptools/"+list_of_samples+"_mBS.heatmap.pearson.pdf",
        "SM/deeptools/"+list_of_samples+"_mBS.heatmap.spearman.pdf",
        "SM/deeptools/"+list_of_samples+"_mBS.scatterplot.pearson.pdf",
        "SM/deeptools/"+list_of_samples+"_mBS.scatterplot.spearman.pdf",
        expand("SM/deeptools/{sample}.q{qv}."+macs2_type+".merged_"+macs2_merge+"bp_.computedMatrix.gz",qv=macs2_q,allow_missing=True),
        expand("SM/deeptools/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.computedMatrix.gz",t=seacr_params_type,cutoff=seacr_params_cutoff,merge=seacr_params_merge,allow_missing=True),
    output:
        tar="SM/analysis_tars/{sample}-analysis.tar.gz",
    shell:
        "tar -czvf {output.tar} {input}"
rule fastqc:
    input:
        fq1=fastq_file_dir+"/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+".fastq.gz",
        fq2=fastq_file_dir+"/"+fastq_prefix+"{sample}"+fastq_mid+"R2"+fastq_suffix+".fastq.gz"
    output:
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+"_fastqc.html",
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+"_fastqc.zip",
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R2"+fastq_suffix+"_fastqc.html",
        "SM/fastqc/"+fastq_prefix+"{sample}"+fastq_mid+"R2"+fastq_suffix+"_fastqc.zip"
    log:
        "SM/fastqc-{sample}-%j.log"
    message:
        "Running fastqc on {wildcards.sample} \\n"
    shell:
        """
        module load fastqc
        [ ! -d SM/fastqc ] && mkdir -p SM/fastqc 
        fastqc {input.fq1} --outdir=SM/fastqc/
        fastqc {input.fq2} --outdir=SM/fastqc/
        """

rule cutadapt:
    input:
        fq1c=fastq_file_dir+"/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+".fastq.gz",
        fq2c=fastq_file_dir+"/"+fastq_prefix+"{sample}"+fastq_mid+"R2"+fastq_suffix+".fastq.gz"
    output:
        fq1_cut="SM/cutadapt/{sample}.cutadapt.R1.fastq.gz",
	fq2_cut="SM/cutadapt/{sample}.cutadapt.R2.fastq.gz"
    log:
        "SM/cutadapt-{sample}-%j.log"
    message:
        "running cutadapt on {wildcards.sample} \\n"
    params:
        "-a " + cutadapt_fwd + " -A " + cutadapt_rev + " "
    shell:
        '''
        module load python/3.6.1
        cutadapt {params} -o {output.fq1_cut} -p {output.fq2_cut} --cores=0 {input.fq1c} {input.fq2c}
        '''

rule align:
    input:
        read1="SM/cutadapt/{sample}.cutadapt.R1.fastq.gz",
        read2="SM/cutadapt/{sample}.cutadapt.R2.fastq.gz"
    output:
        sam=temp("SM/sambam/{sample}.sam"),
        met="SM/sambam/{sample}.metrics",
        bam=protected("SM/sambam/{sample}.bam")
    log:
        "SM/sambam-{sample}-%j.log"
    message:
        "running bowtie2 on {wildcards.sample} \\n"
    params:
        " " + bowtie2_params + " -x $" + assembly + " "
    threads: 20
    shell:
        '''
        module load bowtie2
        module load samtools
        [ ! -d SM/sambam ] && mkdir -p SM/sambam 
        bowtie2 {params} -p {threads}  -1 {input.read1} -2 {input.read2} -S {output.sam} --met-file {output.met}
        samtools view -b -h {output.sam} > {output.bam}
        '''

rule qualtrim:
    input:
        bam="SM/sambam/{sample}.bam"
    output:
        qt="SM/qualtrim_nodup/{sample}.QualTrimmed.bam"
    log:
        "SM/qualtrim_nodup-{sample}_qualtrim-%j.log"
    message:
        "starting qualtrim for {wildcards.sample}\\n"
    params:
        m=mapq
    shell:
        '''
        module load samtools
        module load java
        module load picard
        module load bedtools
        samtools view -h {input.bam} | samtools view -b -h -f 2 | samtools view -h | 
            awk -v MAPQ="{params.m}" 'BEGIN {{NR1=0;}} {{
            if(/^@/){{ print $0; }} else {{ 
            NR1++; 
            if(/XS:i:/){{ iflag=1; }} else {{ iflag=0;}}
            if(NR1%2==1){{
             read1=$0; name_r1=$1; mapq_r1=$5; iflag1=iflag; 
            }} else {{
             read2=$0; name_r2=$1; mapq_r2=$5; iflag2=iflag;
             if(name_r1!=name_r2){{
               echo "Record names \\t" name_r1 "\\t" name_r2 "\\t are not paired. Aborting execution";  
               exit 1;
             }}
             if(mapq_r1>=MAPQ && mapq_r2>=MAPQ){{
              if(iflag1==1 && iflag2==1){{}} else {{print read1 "\\n" read2}}
             }} 
            }}
            }}
            }}' | samtools sort -O BAM -T {wildcards.sample} > {output.qt}
        '''

rule nodup:
    input:
        qt="SM/qualtrim_nodup/{sample}.QualTrimmed.bam"
    output:
        picard="SM/qualtrim_nodup/{sample}.QualTrimmed.Picard.bam",
        metrix="SM/qualtrim_nodup/{sample}.metrix",
        nodup=protected("SM/qualtrim_nodup/{sample}.NoDup.bam")
    log:
        "SM/qualtrim_nodup-{sample}_nodup-%j.log"
    message:
        "starting nodup for {wildcards.sample}\\n"
    params:
        picard=picard_params,
        st_nodup=samtools_nodup_params
    shell:
        '''
        module load java
        module load picard
	module load samtools
        java {params.picard} MarkDuplicates I={input.qt} O={output.picard} M={output.metrix}
        samtools view {params.st_nodup} {output.picard} | samtools sort -n -O BAM -T {wildcards.sample} > {output.nodup}
        ''' 

#sort and index bam by read coordinate
rule sort_index_bam:
    input:
        bam="SM/qualtrim_nodup/{sample}.NoDup.bam"
    output:
        sorted_bam=protected("SM/qualtrim_nodup/{sample}.NoDup.sorted.bam"),
        sorted_bai=protected("SM/qualtrim_nodup/{sample}.NoDup.sorted.bam.bai")
    message:
        "sorting and indexing {wildcards.sample}.NoDup.bam\n"
    log:
        "SM/qualtrim_nodup-sort_index-{sample}-%j.log"
    shell:
        '''
        module load samtools
        samtools sort -O BAM -T {wildcards.sample} {input.bam} -@ 8 > {output.sorted_bam}
        samtools index {output.sorted_bam}
        '''

#process bam creates bed, bdg, bw, normalized bw, and bedpe files for macs2 and seacr

rule process_bam:
    input:
        bam="SM/qualtrim_nodup/{sample}.NoDup.sorted.bam",
        bam_name_sorted="SM/qualtrim_nodup/{sample}.NoDup.bam",
        bai="SM/qualtrim_nodup/{sample}.NoDup.sorted.bam.bai"
    output:
        bed="SM/bed/{sample}.NoDup.bed",
        sorted_bed="SM/bed/{sample}.NoDup.sorted.bed",
        sorted_clipped_bed="SM/bed/{sample}.NoDup.sorted.clipped.bed",
        sorted_clipped_bdg="SM/bdg/{sample}.NoDup.sorted.clipped.bdg",
        sorted_clipped_bw="SM/bw/{sample}.NoDup.sorted.clipped.bw",
        bedpe="SM/bedpe/{sample}.NoDup.sorted.bedpe", #for both macs2 and seacr
        norm_bw=protected("SM/bw/{sample}.NoDup.SeqDepthNorm.bw"),
        bedpe_bdg="SM/bdg/{sample}.NoDup.sorted.bedpe.bdg" #for seacr
    params:
        gF=genomeSizeFile,
        bc=bamCoverage_params,
        bl=blacklist_file
    message:
        "starting processing bam for {wildcards.sample}\n"
    log:
        "SM/process_bams-{sample}-%j.log"
    shell:
        '''
        module load bedtools
        module load ucsc_toolkits
        module load deeptools
        bedtools bamtobed -i {input.bam} > {output.bed} &&
        LC_COLLATE=C sort -k1,1 -k2,2n {output.bed} > {output.sorted_bed} &&
        bedClip {output.sorted_bed} {params.gF} {output.sorted_clipped_bed} && 
        bedtools genomecov -i {output.sorted_clipped_bed} -g {params.gF} -bg > {output.sorted_clipped_bdg} &&
        bedGraphToBigWig {output.sorted_clipped_bdg} {params.gF} {output.sorted_clipped_bw} &&
        bedtools bamtobed -i {input.bam_name_sorted} -bedpe | awk '{{print $1 "\\t" $2 "\\t" $6}}' | sort -k1,1 -k2,2n -k3,3n > {output.bedpe} &&
        bedtools genomecov -bg -i {output.bedpe} -g {params.gF} > {output.bedpe_bdg} &&
        bamCoverage --bam {input.bam} {params.bc} -bl {params.bl} -o {output.norm_bw} 
        '''


rule deeptools_solo:
    input:
        bam="SM/qualtrim_nodup/{sample}.NoDup.sorted.bam"
    output:
        coverage="SM/deeptools/coverage.{sample}.pdf",
        size="SM/deeptools/fragmentSize.{sample}.pdf",
        logsize="SM/deeptools/fragmentSize.{sample}.log.pdf"
    params:
        bpe=bamPE_params,
        c=plotCov_params,
        bl=blacklist_file
    message:
        "starting solo plots for plotCoverage and bamPEFragmentSize for {wildcards.sample}\n"
    log:
        "SM/deeptools-{sample}.solo-%j.log"
    shell:
        '''
        module load deeptools
        [ ! -d SM/deeptools ] && mkdir -p SM/deeptools
        plotCoverage --bamfiles {input.bam} --plotFile {output.coverage} -bl {params.bl} -p max {params.c} &&
        bamPEFragmentSize --bamfiles {input.bam} -o {output.size} -bl {params.bl} -p 8 {params.bpe} &&
        bamPEFragmentSize --bamfiles {input.bam} -o {output.logsize} -bl {params.bl} -p 8 {params.bpe}
        '''
rule deeptools_multi:
    input:
        bam=expand("SM/qualtrim_nodup/{sample}.NoDup.sorted.bam",sample=samples)
    output:
        coverage="SM/deeptools/coverage."+list_of_samples+".pdf",
        size="SM/deeptools/fragmentSize."+list_of_samples+".pdf",
        logsize="SM/deeptools/fragmentSize."+list_of_samples+".log.pdf"
    params:
        bpe=bamPE_params,
        c=plotCov_params,
        bl=blacklist_file       
    message:
        "starting multi plots for plotCoverage, bamPEFragmentSize for "+list_of_samples+"\n"
    log:
        "SM/deeptools-"+list_of_samples+"multi-coverage-bamPE-%j.log"
    shell:
        '''
        module load deeptools
        [ ! -d SM/deeptools ] && mkdir -p SM/deeptools
        plotCoverage --bamfiles {input.bam} --plotFile {output.coverage} -bl {params.bl} -p max {params.c} &&
        bamPEFragmentSize --bamfiles {input.bam} -o {output.size} -bl {params.bl} -p 8 {params.bpe} &&
        bamPEFragmentSize --bamfiles {input.bam} -o {output.logsize} -bl {params.bl} -p 8 {params.bpe}
        '''
        
rule deeptools_multiBamSummary:
    input:
        bam=expand("SM/qualtrim_nodup/{sample}.NoDup.sorted.bam",sample=samples)
    output:
        npz="SM/deeptools/"+list_of_samples+"_mBS.npz"
    params:
        mp=mBS_params,
        bl=blacklist_file       
    message:
        "starting multi plots for multiBamSummary for "+list_of_samples+"\n"
    log:
        "SM/deeptools-"+list_of_samples+"multiBamSummary-%j.log"
    shell:
        '''
        module load deeptools
        [ ! -d SM/deeptools ] && mkdir -p SM/deeptools
        multiBamSummary bins -b {input.bam} -o {output.npz} {params.mp} -bl {params.bl} -p max
        '''
rule deeptools_multiCorr:
    input:
        npz="SM/deeptools/"+list_of_samples+"_mBS.npz"
    output:
        hp="SM/deeptools/"+list_of_samples+"_mBS.heatmap.pearson.pdf",
        hs="SM/deeptools/"+list_of_samples+"_mBS.heatmap.spearman.pdf",
        sp="SM/deeptools/"+list_of_samples+"_mBS.scatterplot.pearson.pdf",
        ss="SM/deeptools/"+list_of_samples+"_mBS.scatterplot.spearman.pdf"
    params:
        cp=plotCorr_params      
    message:
        "starting multi corr plots for plotCorr for "+list_of_samples+"\n"
    log:
        "SM/deeptools-"+list_of_samples+"multiCorr-%j.log"
    shell:
        '''
        module load deeptools
        [ ! -d SM/deeptools ] && mkdir -p SM/deeptools
        plotCorrelation --corData {input.npz} {params.cp} --whatToPlot heatmap --corMethod pearson --plotFile {output.hp} --plotNumbers &&
        plotCorrelation --corData {input.npz} {params.cp} --whatToPlot heatmap --corMethod spearman --plotFile {output.hs} --plotNumbers &&
        plotCorrelation --corData {input.npz} {params.cp} --whatToPlot scatterplot --corMethod pearson --plotFile {output.sp} --log1p &&
        plotCorrelation --corData {input.npz} {params.cp} --whatToPlot scatterplot --corMethod spearman --plotFile {output.ss} --log1p
        ''' 
#########
rule seacr_1of3:
    input:
        bdg="SM/bdg/{sample}.NoDup.sorted.bedpe.bdg"
    output:
        bed="SM/seacr/{sample}_SEACR_{t}-{cutoff}.{t}.bed"
    params:
        p=seacr_path,
        t="{t}",
        c="{cutoff}"
    log:
        "SM/seacr1of3-{sample}_SEACR_{t}-{cutoff}-%j.log"
    message:
        "starting SEACR with type={wildcards.t} and cutoff={wildcards.cutoff} for {wildcards.sample}\\n"
    shell:
        '''
        module load R
        [ ! -d SM/seacr ] && mkdir -p SM/seacr
        bash {params.p} {input.bdg} {params.c} norm {params.t} SM/seacr/{wildcards.sample}_SEACR_{params.t}-{params.c}
        '''
rule seacr_merge_2of3:
    input:
       bed="SM/seacr/{sample}_SEACR_{t}-{cutoff}.{t}.bed"
    output:
       merged_bed="SM/seacr/merged/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.bed"
    params:
        m="{merge}"
    log:
        "SM/seacr-2of3-{sample}_SEACR_{t}-{cutoff}.merge-{merge}-%j.log"
    message:
        "starting merge seacr for {wildcards.sample}, type={wildcards.t}, cutoff={wildcards.cutoff}\\n"
    shell:
        '''
        module load bedtools
        bedtools merge -i {input.bed} -d {params.m} > {output.merged_bed}
        '''
rule seacr_count_3of3:
    input:
        bedpe="SM/bedpe/{sample}.NoDup.sorted.bedpe",
        seacr_merged_bed="SM/seacr/merged/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.bed"
    output:
        seacr_counts="SM/peakcounts/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.counts"
    message:
        "starting peak counts for seacr (type={wildcards.t}, cutoff= {wildcards.cutoff}, merge={wildcards.merge} for {wildcards.sample}\\n"
    log:
        "SM/peakcounts-{sample}.seacr_{t}-{cutoff}.merge-{merge}-%j.log"
    shell:
        '''
        module load bedtools
        bedtools intersect -a {input.seacr_merged_bed} -b {input.bedpe} -sorted -c > {output.seacr_counts} 
        sed -i '/chrY/d' {output.seacr_counts}
        sed -i '/chrM/d;/chrUn/d;/random/d;/alt/d' {output.seacr_counts}
        '''

rule macs2_1of2:
    input:
        bedpe="SM/bedpe/{sample}.NoDup.sorted.bedpe"
    output:
        lam_bdg="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_control_lambda.bdg",
        lam_temp_bed=temp("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_control_lambda.temp.bed"),
        lam_bw="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_control_lambda.bw",
        peak="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_peaks."+macs2_type+"Peak",
        xls="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_peaks.xls",
        #summits="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_summits.bed", #comment out for broad peaks!!
        pileup_bdg="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_treat_pileup.bdg",
        pileup_temp_bed=temp("SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_treat_pileup.temp.bed"),
        pileup_bw="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_treat_pileup.bw"
    params:
        gF=genomeSizeFile,
        q="{qv}",
        g=macs2_genome,
        t=macs2_type,
        b=bOn
    log:
        "SM/macs2-{sample}-q{qv}-%j.log"
    message:
        "starting macs2 call peaks for {wildcards.sample} for q={wildcards.qv}\\n"
    shell:
        '''
        module load bedtools
        module load ucsc_toolkits
        module load macs2
        macs2 callpeak -t {input.bedpe} -q {params.q} -g {params.g} --nolambda -n {wildcards.sample} --format BEDPE -B {params.b} --outdir SM/macs2/{wildcards.sample}.paired.q{params.q}.{params.t} &&
#macs2 callpeak -g 2700000000.0 --name CYC103_2_BAF155 --treatment output/bams/deduped/CYC103_2_BAF155.noMT.filtered.deduped.bam --outdir output/peaks --format BAMPE --nomodel --call-summits --nolambda --keep-dup all -p 0.01 -B --SPMR

        LC_COLLATE=C sort -k1,1 -k2,2n {output.pileup_bdg} | bedtools slop -i stdin -g {params.gF} -b 0 > {output.pileup_temp_bed} &&
        bedClip {output.pileup_temp_bed} {params.gF} {output.pileup_bdg} &&
        bedGraphToBigWig {output.pileup_bdg} {params.gF} {output.pileup_bw} 

        LC_COLLATE=C sort -k1,1 -k2,2n {output.lam_bdg} | bedtools slop -i stdin -g {params.gF} -b 0 > {output.lam_temp_bed} &&
        bedClip {output.lam_temp_bed} {params.gF} {output.lam_bdg} &&
        bedGraphToBigWig {output.lam_bdg} {params.gF} {output.lam_bw}
        '''

rule macs2_count_2of2:
    input:
        macs2_peakfile="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_peaks."+macs2_type+"Peak",
        bedpe="SM/bedpe/{sample}.NoDup.sorted.bedpe"
    output:
        macs2_merged_bed="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/merged_"+macs2_merge+"bp_.sorted.macs2.bed",
        macs2_counts="SM/peakcounts/{sample}.q{qv}."+macs2_type+"/merged_"+macs2_merge+"bp_.macs2.counts"
    params:
        macs2_m=macs2_merge
    message:
        "starting peak counts for macs2 (q{wildcards.qv}) for {wildcards.sample}\\n"
    log:
        "SM/peakcounts-{sample}.macs2_q{qv}-%j.log"
    shell:
        '''
        module load bedtools
        cat {input.macs2_peakfile} >> {output.macs2_merged_bed}.temp.bed &&
        sort -k1,1 -k2,2n {output.macs2_merged_bed}.temp.bed | bedtools merge -i stdin -d {params.macs2_m} > {output.macs2_merged_bed} &&
        bedtools intersect -a {output.macs2_merged_bed} -b {input.bedpe} -sorted -c > {output.macs2_counts} &&
        sed -i '/chrY/d' {output.macs2_counts} &&
        sed -i '/chrM/d;/chrUn/d;/random/d;/alt/d' {output.macs2_counts}
        '''

rule computeMatrix_macs2:
    input:
        macs2_peakfile="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/merged_"+macs2_merge+"bp_.sorted.macs2.bed",
        bw="SM/bw/{sample}.NoDup.SeqDepthNorm.bw"
    output:
        cM_macs2="SM/deeptools/{sample}.q{qv}."+macs2_type+".merged_"+macs2_merge+"bp_.computedMatrix.gz",
    params:
        cm_p=cM_params,
    message:
        "starting compute Matrix  for macs2 (q{wildcards.qv}) for {wildcards.sample}\\n"
    log:
        "SM/deeptools-computeMatrix-{sample}.macs2_q{qv}-%j.log"
    shell:
        '''
        module load deeptools
        computeMatrix reference-point -S {input.bw} -R {input.macs2_peakfile} {params.cm_p} -o {output.cM_macs2}
        '''

rule computeMatrix_seacr:
    input:
        seacr_peakfile="SM/seacr/merged/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.bed",
        bw="SM/bw/{sample}.NoDup.SeqDepthNorm.bw"
    output:
        cM_seacr="SM/deeptools/{sample}_SEACR_{t}-{cutoff}.merge-{merge}bp.computedMatrix.gz",
    params:
        cm_p=cM_params,
    message:
        "starting compute Matrix for seacr  (type={wildcards.t}, cutoff= {wildcards.cutoff}, merge={wildcards.merge} for {wildcards.sample}\\n"
    log:
        "SM/deeptools-computeMatrix-{sample}.seacr_{t}-{cutoff}.merge-{merge}-%j.log"
    shell:
        '''
        module load deeptools
        computeMatrix reference-point -S {input.bw} -R {input.seacr_peakfile} {params.cm_p} -o {output.cM_seacr}
        '''
   









##rule find_memes:
##    input:
##        seacr_peakfile="SM/seacr/{sample}_SEACR_{t}-{cutoff}.{t}.bed"
##        macs2_peakfile="SM/macs2/{sample}.paired.q{qv}."+macs2_type+"/{sample}_peaks."+macs2_type+"Peak",
##    output:
##        meme_chip=directory("SM/meme/{sample}_MEME/"),
##        seacr_filt_peak=temp("SM/meme/{sample}_MEME/{sample}.seacr.filt.peak"),
##        macs2_filt_peak=temp("SM/meme/{sample}_MEME/{sample}.macs2.filt.peak"),
##    params:
##        bl=blacklist_file,
##        bedops_bin=config["meme"]["bedops_bin"],
##        meme_chip_bin=config["meme"]["meme_bin"],
##        meme_chip_params=config["meme"]["params"],
##        genome_sequence=config["meme"]["genomesequence"]
##    shell:
##        '''
##        module load bedtools
##        cat {seacr_peakfile} | grep -v -e "chrM" | {params.bedops_bin}/sort-bed - | {params.bedops_bin}/bedops -n 1 - {params.bl} > {output.seacr_filt_peak} \
##            sort -t "   " -g -k8 -r | head -n 5000 | shuf | head -n 1000 | {params.bedops_bin}/sort-bed - > {output.seacr_filt_peak} &&
##        cat {macs2_peakfile} | grep -v -e "chrM" | {params.bedops_bin}/sort-bed - | {params.bedops_bin}/bedops -n 1 - {params.bl}  \
##            sort -t "   " -g -k8 -r | head -n 5000 | shuf | head -n 1000 | {params.bedops_bin}/sort-bed - > {output.macs2_filt_peak} &&
##
##        bedtools getfasta -fi {params.genome_sequence} -bed {input.{params.meme_chip_bin} -oc {output.meme_chip} {params.meme_chip_params} 

