#!/bin/sh
# properties = {"type": "single", "rule": "filter_L1", "local": false, "input": ["/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/tmp_data/qc_raw_data/mouse_Library1_raw_qc_primer_R1.fastq.gz", "/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/tmp_data/qc_raw_data/mouse_Library1_raw_qc_primer_R2.fastq.gz"], "output": ["/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/tmp_data/qc_raw_data/mouse_Library1_raw_qc_linker1_R1.fastq.gz", "/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/tmp_data/qc_raw_data/mouse_Library1_raw_qc_linker1_R2.fastq.gz"], "wildcards": {"sample": "mouse_Library1"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 89486, "mem_mib": 85341, "disk_mb": 89486, "disk_mib": 85341, "tmpdir": "<TBD>"}, "jobid": 4, "cluster": {"partition": "general", "job-name": "filter_L1_mouse_Library1", "ntasks": "1", "cpus-per-task": "20", "mem": "64g", "time": "120:00:00", "output": "/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/log/filter_L1.%j.out", "error": "/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/log/filter_L1.%j.err"}}
cd '/rwthfs/rz/cluster/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing' && /work/rwth1209/enviroments/dbitseq/bin/python3.11 -m snakemake --snakefile '/rwthfs/rz/cluster/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing/Snakefile' --target-jobs 'filter_L1:sample=mouse_Library1' --allowed-rules 'filter_L1' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=89486' 'mem_mib=85341' 'disk_mb=89486' 'disk_mib=85341' --wait-for-files '/rwthfs/rz/cluster/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing/.snakemake/tmp.hmtw7437' '/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/tmp_data/qc_raw_data/mouse_Library1_raw_qc_primer_R1.fastq.gz' '/rwthfs/rz/cluster/hpcwork/rwth1209/data/dbitseq/atac/mouse/230522/data_processing/mouse_Library1/CellRanger/tmp_data/qc_raw_data/mouse_Library1_raw_qc_primer_R2.fastq.gz' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'mtime' 'input' 'software-env' 'params' 'code' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/work/rwth1209/enviroments/dbitseq/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/rwthfs/rz/cluster/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing/.snakemake/tmp.hmtw7437/4.jobfinished' || (touch '/rwthfs/rz/cluster/work/rwth1209/software/DbiT-seq-pipeline/Data_preprocessing/.snakemake/tmp.hmtw7437/4.jobfailed'; exit 1)
