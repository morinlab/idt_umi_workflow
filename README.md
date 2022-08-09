# IDT UMI workflow
Snakemake workflow for processing CAPPSeq data sequenced using IDT's UMIs. Note that this workflow requires paired-end FASTQ files which have been prepared using IDT's UMI structure


## Configuration
1. Download the workflow: `git clone https://github.com/morinlab/idt_umi_workflow`
2. Edit the sample file under `config/samplelist.tsv` with your samples of interest
  - This will be a 3-column config file containing `sample_id`, `path_to_R1_fastq`, `path_to_R2_fastq`
  - Use "#" to comment out lines
3. Edit the workflow config file under `config/cappseq_umi_config.yaml` for your project's parameters
4. Create a conda environment containing `snakemake` version 7 or newer, and activate that environment
5. Run the workflow using `snakemake --use-conda -s cappseq_umi_workflow.smk -j <number_of_threads>`

## Workflow
![CAPPSeq_workflow](https://user-images.githubusercontent.com/22357668/183552741-0842ba1e-19b6-4a64-8aae-4651dbbaa1a8.png)

## Results
  - Final processed BAM files (collapsed and error corrected using UMIs) can be found under `99-final/`
  - Final QC report can be found under `Q9-multiqc`

