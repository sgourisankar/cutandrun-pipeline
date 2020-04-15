# cutandrun-pipeline
This is a pipeline I wrote using Python's Snakemake to automatically QC, align, filter, call peaks, and generate plots for data generated from paired end CUTANDRUN experiments.

It requires a cluster with slurm (although you can edit the run-SM-workflow file to use a different jobs manager, such as SGE/Torque), and is currently written to use the "module load" function to run various prerequisites such as bedtools, macs2, samtools, picard, etc. 

The pipeline will run both seacr and macs2 to call peaks with as many parameters (stringent vs relaxed, different cutoffs and base pairs of merging nearby peaks, q values) as you'd like. It will also generate diagnostic plots using deeptools: the size distribution of fragments, the coverage, and the read depth correlation between files. 

To use it, first generate fastq.gz files. 

Move them to your SCRATCH workspace along with the pipeline (all files in the pipeline should be in the same directory, but need not be in the same directory as the fastq.gz files). 

Then, fill out the config.json file and the cluster.json files with the appropriate parameters. Instructions for both are below. 

Then, run "run-SM-workflow.sh" by using $ bash run-SM-workflow.sh 
You'll first generate a dry run of what commmands will be executed, and a graph of the dependencies in an .svg format.. Check the graph in your browser to make sure all the files are being processed appropriately. If there are errors at this point, it is most likely because the parameters in the config.json file (sample names, parameters, etc.) were not appropriately filled in. 

Then, type "Y" to run the pipeline. Jobs will be queued for submission automatically, and output will be appended to the file "nohup.out". You can check progression periodically by, for example, $ tail nohup.out -n 100 
(which will display the last 100 lines in nohup.out). 

You will know the run is concluded when a) each of your samples has a tar.gz archive created in the folder "analysis_tars" with all of its processed files; and b) the rule "all" has been run (which you can see in nohup.out). 

######HOW TO EDIT config.json#########
The Number One error here is forgetting commas at the end of each line but the last one in a {} group. If you get an error, check this first. 

Make extra sure to change anything below labeled "->"!!

-> 1. "fastq_file_dir": directory with fastq files in it. 
-> 2. "samples", "fastq_prefix", "fastq_mid", and "fastq_suffix": these use the following format: 
    fastq_file_dir+"/"+fastq_prefix+"{sample}"+fastq_mid+"R1"+fastq_suffix+".fastq.gz" (or R2)
    In "samples", you can have a comma-separated list [] of all your samples you'd like to 
    process. All following files will be named with just the name of the sample. 
-> 3. "assembly": mm10, mm9, hg38, hg19, etc. Used as a parameter in bowtie2. 
4. "bowtie2_params": the rest of the bowtie2 parameters
5. "qualtrim_mapq": the minimum mapq quality score for filtering reads
6. "picard_params": Picard parameters INCLUDING the path of the picard.jar file.
7. "samtools_nodup_params": parameters to process picard-marked duplicates and sort them
-> 8. "genomeSizeFile": path of genome size file for appropriate organism
-> 9. "bamCoverage_params": used to create depth-normalized bigwigs. Make sure to change the 
      effectiveGenomeSize parameter for the appropriate organism.
-> 10. "seacr_params": make sure to change the path to the .sh file for SEACR where you have it. 
       Here you can list all the different trials of seacr parameters you'd like, where 
       "type" can be "relaxed" and/or "stringent"; "cutoff" can be any list of cutoffs you'd like 
       to try; "merge" is any list of base-pairs you'd like to try merging peaks. Watch out - the
       more you list here, the longer your analysis will take (e.g. it will submit (# of types x # 
       of cutoffs x # of merge) total jobs) 
-> 11. "macs2_params": "q" is any list of q-values you'd like to try. Make sure to change "genome" 
       to the correct organism. "type" can EITHER be "narrow" or "broad", you cannot list both. 
       "merge" can only be one value for the bp to merge peaks, it cannot be a list. 
-> 12. "deeptools": these are parameters for the deeptools functions run. "bamPE_params" are 
       "parameters" for bamPEFragmentSize; "plotCov_params" are for plotCoverage; 
       "multiBamSum_params" are for multiBamSummary; and "plotCorr_params" are for 
       plotCorrelation. Look at deeptools docs for an explanation. 
       You can change colors for example of the plots here. 
       Make sure to change "blacklist_file" to the path where you have the appropriate organism's 
       blacklist bed file of regions to skip. 
       

######HOW TO EDIT cluster.json#########
Under "__default__",
1. "mem" is --mem-per-cpu 
2. "time" is -t (time)
3. "email" is --mail-user (if you'd like to be emailed with job updates)
4. "mail_type" is --mail-type and can be ALL, FAIL, BEGIN, or END
5. Don't change "error" or "output"
You can change job resources for each rule in a similar format, to suit how large your files are.

######OUTPUT of Pipeline#########
The main directory created will be __/SM__
__SM/fastqc__ has fastqc files for all samples processed
__SM/sambam__ has .sam, .metrics, and .bam files after alignment
__SM/qualtrim_nodup__ has: 
    QualTrimmed.bam = after quality trimming
    QualTrimmed.Picard.bam = after quality trimming and Picard marking,
    NoDup.bam is a read-name sorted file after filtering duplicates and low quality reads
    NoDup.sorted.bam(.bai) is a coordinate-sorted version (and index) of No.Dup.bam
    (all further analyses use NoDup.sorted.bam)
__SM/bed__ has sorted and clipped bedfiles
__SM/bdg__ has sorted and clipped bedgraphs
__SM/bedpe__ has sorted bedpe files for macs2 and seacr 
__SM/bw__ has depth normalized bigwigs (after sorging and clipping)
__SM/seacr__ has all the seacr beds and merged beds
__SM/macs2__ has all the macs2 files
__SM/peakcounts__ has all the counts files from macs2 and seacr peak calling
__SM/deeptools__ has all the depptools plots. Specifically, 
    coverage.{sample}.pdf is depth from plotCoverage for one sample
    coverage."{sample1-sample2-sample3...}.pdf is the same but with all the samples together
    fragmentSize.{sample}.pdf/.log.pdf is log-scale or normal scale fragment size distribution for 
      that sample
    fragmentSize.{sample1-sample2-sample3...}.pdf/.log.pdf is the same but w/ all together
    .npz is for deeptools, an internal matrix of samples
    .heatmap/scatterplot.pearson/spearman are from multiBamSummary for all the samples comparison
__SM/analysis_tars__ has, for each sample, a tar.gz archive of all the relevant plots and files. 

######What you can change#########
Snakemake will automatically detect if any files are missing or changed, and rerun the appropriate jobs. If the dry run and graph doesn't show those rules executing, try $ touch FILE.TO.TOUCH 
You can always directly run rules as targets by changing the $ snakemake  command in the run-SM-workflow.sh file. 
######Updates to make#########
Probably none to this besides fixing bugs. The python file is already too long/unwieldy and for downstream analyses like motif finding etc I'll just write new workflows. 

