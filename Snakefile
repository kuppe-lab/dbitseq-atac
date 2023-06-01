import os
import fnmatch

########################################################
# global varibales
# 1. modify global varibales here
# 2. modify cluster.json
########################################################

rawData_dir = "/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/raw" # raw sequencing data directory
processedData_dir = '/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing' # processed sequencing data directory
fastq_ext = '.fastq.gz' # fastq file extension
# Make sure read 1 has R1 in its name and read 2 R2, otherwise you need to change the fnmatch below
core = 20 # threads to use
core_cellRanger = int(core*0.8) # threads to use in cellRanger
mem = 64 # memory to use in cellRanger

bbduk_dir = "/rwthfs/rz/cluster/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing/bbmap/bbduk.sh" # bbduk.sh directory
cellranger_dir = "/rwthfs/rz/cluster/work/rwth1209/software/cellranger-atac-2.1.0-dbitseq/cellranger-atac" # cellranger-atac directory
mouse_dir = "/rwthfs/rz/cluster/hpcwork/rwth1209/reference_genome/mouse_mm10_atac_10x" # reference data directory
human_dir ="/rwthfs/rz/cluster/hpcwork/rwth1209/reference_genome/human_GRCh38_atac_10x"
########################################################

########################################################
# Do not change codes below
########################################################
#shell.executable("/usr/local_rwth/bin/zsh")
#shell.prefix("source /home/mz637064/mambaforge/etc/profile.d/conda.sh; ")
samples = os.listdir(rawData_dir)
output = [processedData_dir+'/'+sample+'/CellRanger/'+sample+'/outs/fragments.tsv.gz' for sample in samples]

def get_genome(sample):
    if "human" in sample.lower():
        return human_dir
    elif "mouse" in sample.lower():
        return mouse_dir
    else:
        raise ValueError("{} does not contain human or mouse in name".format(sample))

list_fastq={}
for sample in samples:
    list_fastq[sample]=(fnmatch.filter(os.listdir(rawData_dir+'/'+sample), '*R1*'+fastq_ext)[0].split(".")[0],
                        fnmatch.filter(os.listdir(rawData_dir+'/'+sample), '*R2*'+fastq_ext)[0].split(".")[0])
                        
# generating folders

for sample in samples:
    
    output_dir=processedData_dir+'/'+sample+'/CellRanger'
    tmp_data=output_dir+'/tmp_data'
    qc_raw_data=tmp_data+'/qc_raw_data'
    log_dir=output_dir+'/log'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(tmp_data):
        os.makedirs(tmp_data)
    if not os.path.exists(qc_raw_data):
        os.makedirs(qc_raw_data)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)


rule all:
  input:
    output


rule filter_primer:
  input:
    in1 = lambda wildcards: rawData_dir+"/{sample}/"+list_fastq[wildcards.sample][0]+fastq_ext,
    in2 = lambda wildcards: rawData_dir+"/{sample}/"+list_fastq[wildcards.sample][1]+fastq_ext
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_primer_R1.fastq.gz",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_primer_R2.fastq.gz"
  shell:
    '''
    {bbduk_dir} \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k=22 mm=f rcomp=f restrictleft=30 skipr1=t \
    hdist=2 \
    stats={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data/qc_raw_data/{wildcards.sample}_stats.primer.txt \
    threads={core} \
    literal=CAAGCGTTGGCTTCTCGCATCT
    '''    
rule filter_L1:
  input:
    in1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_primer_R1.fastq.gz",
    in2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_primer_R2.fastq.gz"
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R1.fastq.gz",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R2.fastq.gz"
  shell:
    '''
    {bbduk_dir} \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k=30 mm=f rcomp=f restrictleft=108 skipr1=t \
    hdist=3 \
    stats={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data/qc_raw_data/{wildcards.sample}_stats.linker1.txt \
    threads={core} \
    literal=GTGGCCGATGTTTCGCATCGGCGTACGACT
    '''
    
rule filter_L2:
  input:
    in1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R1.fastq.gz",
    in2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_linker1_R2.fastq.gz"
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R1.fastq.gz",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R2.fastq.gz"
  shell:
    '''
    {bbduk_dir} \
    in1={input.in1} \
    in2={input.in2} \
    outm1={output.out1} \
    outm2={output.out2} \
    k=30 mm=f rcomp=f restrictleft=70 skipr1=t \
    hdist=3 \
    stats={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data/qc_raw_data/{wildcards.sample}_stats.linker2.txt \
    threads={core} \
    literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG
    '''

rule bc_process:
  input:
    processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R2.fastq.gz"
  output:
    out1 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R1_001.fastq",
    out2 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R2_001.fastq"
#  conda:
#  		"dbitseq.yaml"
  shell:
    '''
   # conda activate /work/rwth1209/enviroments/dbitseq
    /work/rwth1209/enviroments/dbitseq/bin/python3 BC_process.py --input {input} --output_R1 {output.out1} --output_R2 {output.out2}
    '''
    
rule R1_rename:
  input:
    processedData_dir + "/{sample}/CellRanger/tmp_data/qc_raw_data/{sample}_raw_qc_R1.fastq.gz"
  output:
    processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R3_001.fastq.gz"
  shell:
    '''
    cp {input} {output}
    '''
   
    
rule cell_ranger:
  input:
    in1 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R1_001.fastq",
    in2 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R2_001.fastq",
    in3 = processedData_dir + "/{sample}/CellRanger/tmp_data/{sample}_S1_L001_R3_001.fastq.gz",
    genome = lambda wildcards: get_genome(wildcards.sample)
  output:
    processedData_dir + "/{sample}/CellRanger/{sample}/outs/fragments.tsv.gz"
  shell:
    '''
    cd {processedData_dir}/{wildcards.sample}/CellRanger/
    rm -r {wildcards.sample}
    
    {cellranger_dir} count \
    --id={wildcards.sample} \
    --reference={input.genome}\
    --fastqs={processedData_dir}/{wildcards.sample}/CellRanger/tmp_data \
    --sample={wildcards.sample} \
    --localcores={core_cellRanger} \
    --localmem={mem}
    '''
