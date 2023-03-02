import os
import os.path
from os import listdir
from os.path import isfile, join

################################################################################################################################################################################################################################################################################################################################################
################################################################### USEFULL FUNCTIONS

    
    
######## GET LOCATION OF SAMPLE ################################################

def LocationOfPairedSample1(Sample):
    FASTQ=RUNS_TO_FASTQS[Sample][0]
    
    if FASTQ in DIR_LOC_FASTQ.keys():   ## Is sample available locally?
        LOC=DIR_LOC_FASTQ[FASTQ]
        
    if FASTQ not in DIR_LOC_FASTQ.keys(): ## If sample not in any of the folders given, will try to download it
        LOC='NEWLY_DOWNLOADED_FASTQS/PE1/'
        FASTQ=F'{Sample}_1.fastq.gz'
    
    return "{}{}".format(LOC,FASTQ)
    
    
    
def LocationOfPairedSample2(Sample):
    FASTQ=RUNS_TO_FASTQS[Sample][1]
    
    
    if FASTQ in DIR_LOC_FASTQ.keys():   ## Is sample available locally?
        LOC=DIR_LOC_FASTQ[FASTQ]
        
    if FASTQ not in DIR_LOC_FASTQ.keys(): ## If sample not in any of the folders given, will try to download it
        LOC='NEWLY_DOWNLOADED_FASTQS/PE2/'
        FASTQ=F'{Sample}_2.fastq.gz'
    
    return "{}{}".format(LOC,FASTQ)



def LocationOfNonPairedSample(Sample):
    FASTQ=RUNS_TO_FASTQS[Sample][0]
    
    if FASTQ in DIR_LOC_FASTQ.keys():   ## Is sample available locally?
        LOC=DIR_LOC_FASTQ[FASTQ]
    
    if FASTQ not in DIR_LOC_FASTQ.keys(): ## If sample not in any of the folders given, will try to download it
        LOC='NEWLY_DOWNLOADED_FASTQS/NPE/'
        FASTQ=F'{Sample}.fastq.gz'
    
    
    return "{}{}".format(LOC,Sample)







########### DECIDE WHICH PATH OF THE PIPELINE TO TAKE: PAIR END OR NOT  ############################################

def SAM_PE_OR_NOT(Sample):

    if PAIRED_END[Sample]=='YES':
        FILE="{}_PE.sam".format(Sample)
    if PAIRED_END[Sample]=='NO':
        FILE="{}_NONPE.sam".format(Sample)
        
    return FILE



########### Function to return list of sorted files corresponding to each sample

def GetSortedForSample(Sample):
    
    SORTED_BAMS=SAMPLE_TO_RUNS[Sample]
    SORTED_BAMS=[ F'{X}.sorted' for X in SORTED_BAMS ]
    return SORTED_BAMS














################################################################################################################ SET UP ################################################################################################################################################################################################################################







SAMPLES=[]
FASTQ_SAMPLES=[]
SORTED_SAMPLES=[]
SORTED_SAMPLES_WITH_NO_MERGING=[]


DIR_LOC_FASTQ={}
DIR_LOC_SORTED={}
PAIRED_END={}
SAMPLE_TO_RUNS={}
RUNS_TO_FASTQS={}
RUNS_TO_DOWNLOAD_LINKS={}
SAMPLE_TO_SPECIES_ID={}


FASTQ_FOLDERS=['../PONGO_STUDY_FASTQ/'] ##### Set by user, should end with a '/', If you already have the Fastq files somewhere in your machine, add the link here so they don't need to be downloaded




#### REFERENCE TO MAP TO
REF_LOC='../../REFERENCE_GENOMES/'      ##### Set by user
REF=REF_LOC+'Pongo_abelii.Susie_PABv2.dna.toplevel.fa' ##### Set by user



############################################################################################################################################################################################################################
############################################################################################################################################################################################################################
##### The most important file!
if os.path.exists('METADATA'):          #### IMPORTANT NOTE, READ BELOW
    META_DATA_FILE=open('METADATA','r') #### Set by user, should be a table file from ENA. Required, tab seperated columns: 'sample_accession', 'run_accession', 'scientific_name', 'fastq_ftp'
if os.path.exists('METADATA')!=True:
    print('Error - No Metadata txt file identified!')
    

if (os.path.exists('NEWLY_DOWNLOADED_FASTQS'))!=True:
    os.mkdir('NEWLY_DOWNLOADED_FASTQS')


if (os.path.exists('NEWLY_DOWNLOADED_FASTQS'))==True:

    if (os.path.exists('NEWLY_DOWNLOADED_FASTQS/PE1'))!=True:
        os.mkdir('NEWLY_DOWNLOADED_FASTQS/PE1')
        
    if (os.path.exists('NEWLY_DOWNLOADED_FASTQS/PE2'))!=True:
        os.mkdir('NEWLY_DOWNLOADED_FASTQS/PE2')
        
    if (os.path.exists('NEWLY_DOWNLOADED_FASTQS/NPE'))!=True:
        os.mkdir('NEWLY_DOWNLOADED_FASTQS/NPE')












if os.path.exists('METADATA'): ##### Should check labels of METADATA file to see if the required ones are there and in which order

    LABELS=META_DATA_FILE.readline().strip().split('\t')
    

    
    ## Read METADATA file, create Dictionary

    for LINE in META_DATA_FILE:
        
        LINE=LINE.strip().split('\t')
        
        
        SAMPLE_HERE=LINE[LABELS.index('sample_accession')]
        SPECIES_HERE=LINE[LABELS.index('scientific_name')]
        RUN_ACCESSION_HERE=LINE[LABELS.index('run_accession')]
        FASTQS_LINKS = LINE[ LABELS.index('fastq_ftp') ].split(';')
        
        
        FASTQS_HERE= [ X.split('/')[len(X.split('/'))-1] for X in FASTQS_LINKS] 
        
        
        
        
        
        
        RUNS_TO_DOWNLOAD_LINKS[RUN_ACCESSION_HERE] = FASTQS_LINKS
        
        RUNS_TO_FASTQS[RUN_ACCESSION_HERE] = FASTQS_HERE
        
        if len(FASTQS_HERE)==2:
            PAIRED_END[RUN_ACCESSION_HERE] = 'YES'
        if len(FASTQS_HERE)==1:
            PAIRED_END[RUN_ACCESSION_HERE] = 'NO'
        if len(FASTQS_HERE)==3:  #### If 3 files present, usually the first one is non Pair end and the 2 next ones are paired end, so this will grab the first one  
            PAIRED_END[RUN_ACCESSION_HERE] = 'NO'
            
            
            
        SAMPLE_TO_SPECIES_ID[ SAMPLE_HERE ] = SPECIES_HERE
        
        
        
        
        if SAMPLE_HERE in SAMPLE_TO_RUNS.keys():
            SAMPLE_TO_RUNS[SAMPLE_HERE].append(RUN_ACCESSION_HERE)
        else:
            SAMPLE_TO_RUNS[SAMPLE_HERE]=[RUN_ACCESSION_HERE]
            SAMPLES.append(SAMPLE_HERE)


##### Check which Fastqs are available locally!

for FASTQ_PATH in FASTQ_FOLDERS: #### Find which samples exist in fastq folder, make 'library' of them

    FASTQ_FILES = [f for f in listdir(FASTQ_PATH) if isfile(join(FASTQ_PATH, f))]
    FASTQ_FILES = [f for f in FASTQ_FILES if '.fastq' in f]
    for FILE in FASTQ_FILES:    
        DIR_LOC_FASTQ[FILE]=FASTQ_PATH
        




############################################################################################################################################################################################################################################################


















############################################################################################################################################################################################################################################################

####### FINAL SAMPLE SET UP
SAMPLES=list(set(SAMPLES))
SAMPLES.sort()
SAMPLES=SAMPLES[1] ############################### FOR TESTING!


###### One Rule to Rule them ALL

rule all:
    input:
        REF+'.bwt',
        # expand("{sample}.sorted", sample=FASTQ_SAMPLES),###### Initial Sorting Step
        expand("{sample}.bam", sample=SAMPLES),    ############## Combined into samples
        expand("{sample}.bam.bai", sample=SAMPLES) ############ Combined into samples

























######################################################################################################################
#### SAMPLE IDENTIFICATION - DOWNLOAD



##### DOWNLOAD LINKS

#### PE-1
rule Download_FastQ1_From_Link:
    input:
        FastQ_Link='METADATA'
    output:
        FastQ_File='NEWLY_DOWNLOADED_FASTQS/PE1/{sample}_1.fastq.gz'
    run:
        DL_LINK=''
        
        if wildcards.sample in RUNS_TO_DOWNLOAD_LINKS.keys():
            DL_LINK=RUNS_TO_DOWNLOAD_LINKS[wildcards.sample][0]
            
        shell(F"wget --continue --progress=dot:mega --tries=0 {DL_LINK} -O {output.FastQ_File}")


#### PE-2
rule Download_FastQ2_From_Link:
    input:
        FastQ_Link='METADATA'
    output:
        FastQ_File='NEWLY_DOWNLOADED_FASTQS/PE2/{sample}_2.fastq.gz'
    run:
        DL_LINK=''
        
        if wildcards.sample in RUNS_TO_DOWNLOAD_LINKS.keys():
            DL_LINK=RUNS_TO_DOWNLOAD_LINKS[wildcards.sample][1]
            
        shell(F"wget --continue --progress=dot:mega --tries=0 {DL_LINK} -O {output.FastQ_File}")


#### NON-PE
rule Download_FastQNPE_From_Link:
    input:
        FastQ_Link='METADATA'
    output:
        FastQ_File='NEWLY_DOWNLOADED_FASTQS/NPE/{sample}.fastq.gz'
    run:
        DL_LINK=''
        
        if wildcards.sample in RUNS_TO_DOWNLOAD_LINKS.keys():
            DL_LINK=RUNS_TO_DOWNLOAD_LINKS[wildcards.sample][0]
            
        shell(F"wget --continue --progress=dot:mega --tries=0 {DL_LINK} -O {output.FastQ_File}")







######################################################################################################################
#### REF PREP

rule Index_Reference:
    input:
        REF
    output:
        REF+'.bwt'
    shell:
        "bwa index {input}"
        
        
        
        
        
        
        
        


###################################################################################################################### SAMPLE PREP ############################################################################################################################
######### TRIMMING


###### PE
rule Adapter_trimming_PE:
    input:
        r1=lambda wildcards: LocationOfPairedSample1(wildcards.sample),
        r2=lambda wildcards: LocationOfPairedSample2(wildcards.sample)

    output:
        o1=temp('{sample}_1_FASTP.fastq.gz'),
        o2=temp('{sample}_2_FASTP.fastq.gz')
    threads: 16
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.o1}  -O {output.o2} --length_required 30  --thread {threads}"



######## NON PE

rule Adapter_trimming_NON_PE:
    input:
        R= lambda wildcards: LocationOfNonPairedSample(wildcards.sample)
    output:
        TR=temp("{sample}.trimmed.fastq.gz")
    threads: 16
    shell:
        "fastp -i {input.R} -o {output.TR} --length_required 30 --thread {threads}"










########################################################################################################################### ALIGNMENT ###########################################################################################################################
##### ALIGNMENT



######## PE 
rule BWA_PE_aln_1:
    input:
        r1='{sample}_1_FASTP.fastq.gz',
        R=REF,
        RI=REF+'.bwt'
    output:
        o1=temp("{sample}_1.sai")
    threads: 16
    shell:
        "bwa aln -t {threads} {input.R} {input.r1} > {output.o1}"
   


rule BWA_PE_aln_2:
    input:
        r2='{sample}_2_FASTP.fastq.gz',
        R=REF,
        RI=REF+'.bwt'
    output:
        o2=temp("{sample}_2.sai")
    threads: 16
    shell:
        "bwa aln -t {threads} {input.R} {input.r2} > {output.o2}"



rule BWA_PE_sampe:
    input:
        r1="{sample}_1.sai",
        r2="{sample}_2.sai",
        raw1='{sample}_1_FASTP.fastq.gz',
        raw2='{sample}_2_FASTP.fastq.gz',
        R=REF,
        RI=REF+'.bwt'
    output:
        SAM=temp("{sample}_PE.sam")
    threads: 16
    shell:
        "bwa sampe {input.R} {input.r1} {input.r2} {input.raw1} {input.raw2} > {output.SAM}"



######## NON PE 
### Create SAM file, maybe switch to bwa mem for modern DNA?

rule BWA_aln_NON_PE:
    input:
        TR="{sample}.trimmed.fastq.gz",
        R=REF,
        RI=REF+'.bwt'
    output:
        SAI=temp("{sample}.sai")
    threads: 16
    shell:
        "bwa aln -t {threads} {input.R} {input.TR} > {output.SAI}"
    

rule BWA_samse_NON_PE:
    input:
        SAI="{sample}.sai",
        R=REF,
        RI=REF+'.bwt',
        TR="{sample}.trimmed.fastq.gz"
    output:
        SAM=temp("{sample}_NONPE.sam")
    shell:
        "bwa samse {input.R} {input.SAI} {input.TR} > {output.SAM}"










################################################################################## BAM CONVERSION AND CLEAN-UP $#################################################################################################################
####### BAM convertion






rule samtools_convert_to_bam:
    input:
        SAM=lambda wildcards: SAM_PE_OR_NOT(wildcards.sample)
    output:
        BAM=temp("{sample}.initial")
    threads: 4
    shell:
        "samtools view -S --threads={threads} -b {input.SAM} > {output.BAM}"


rule samtools_fix_bam:
    input:
        BAM="{sample}.initial"
    output:
        FIXED_BAM=temp("{sample}.fixed")
    threads: 4
    shell:
        "samtools fixmate -@ {threads} -m {input.BAM} {output.FIXED_BAM}"


rule samtools_sort:
    input:
        FIXED_BAM="{sample}.fixed"
    output:
        SORTED_BAM="{sample}.sorted"
    threads: 4
    shell:
        "samtools sort {input.FIXED_BAM}  -@ {threads} -o {output.SORTED_BAM}"









################################################################################## LIBRARY MERGING AND FINAL FILTERING #################################################################################################################

### Merge Different BAM files if they correspond to the same sample

rule Combine_Same_Sample_Bams:
    input:
        SORTED_BAMS= lambda wildcards: GetSortedForSample(wildcards.sample)
    output:
        COMBINED=temp('{sample}.combined')
    threads: 4
    run:
        RENAMED_SORTED_BAMS=' '.join(input.SORTED_BAMS)
        shell('samtools merge --threads {threads} -o {output.COMBINED} -f {input.SORTED_BAMS}')


### Markdup BAM files for different things: Mark supplementary reads of duplicates as duplicates, remove duplicates

rule samtools_markedup_bam:
    input:
        SORTED_BAM="{sample}.combined"
    output:
        MARKED_BAM=temp("{sample}.marked")
    threads: 4
    shell:
        "samtools markdup --mode t -S -r -@ {threads} {input.SORTED_BAM} {output.MARKED_BAM}"


### Filter BAM files for MAPQ of >30 and 

rule samtools_final_filter_bam:
    input:
        MARKED_BAM="{sample}.marked"
    output:
        FILTERED_BAM="{sample}.filtered"
    threads: 4
    shell:
        "samtools view -@ {threads} -q 15 -F 0x404 -b -o {output.FILTERED_BAM} {input.MARKED_BAM}"


### Index the final bam file


rule samtools_index:
    input:
        FILTERED_BAM="{sample}.filtered"
    output:
        BAI="{sample}.filtered.bai"
    threads: 4
    shell:
        "samtools index -b {input.FILTERED_BAM} -@ {threads}"


### Rename the BAM file

rule rename_files:
    input:
        FILTERED_BAM="{sample}.filtered",
        FILTERED_BAI="{sample}.filtered.bai"
    output:
        RENAMED_BAM="{sample}.bam",
        RENAMED_BAI="{sample}.bam.bai"
    run:
        shell("mv {input.FILTERED_BAM} {output.RENAMED_BAM}")
        shell("mv {input.FILTERED_BAI} {output.RENAMED_BAI}")