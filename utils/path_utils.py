import pandas as pd
import os


class FastQInfo(object):
    """"""

    def __init__(self, samplesheet, fastq_dirs):
        super(FastQInfo, self).__init__()
        self.samplesheet = pd.read_csv(samplesheet, sep='\t', comment='#')
        self.sample_ids = self.samplesheet['sample_id'].tolist()
        self.R1_fastsq, self.R2_fastsq = self.get_fastq_paths(fastq_dirs, self.sample_id)
        self.sampleID_to_Run = self.match_sampleID_to_Run(self.samplesheet)

    
    def get_fastq_paths(self, fastq_dirs: list, sample_ids: list) -> dict:

        # find all fastq files in all listed fastq_dirs
        all_fastqs = [file for dir in fastq_dirs for file in os.listdir(dir) if file.endswith('.fastq.gz')]
        # make dict with sample_id as key and R1 fastq path as value
        R1_paths = {sample: file for sample in sample_ids for file in all_fastqs if sample in file and '_R1_' in file}
        # make dict with sample_id as key and R2 fastq path as value
        R2_paths = {sample: file for sample in sample_ids for file in all_fastqs if sample in file and '_R2_' in file}

        return R1_paths, R2_paths

    def match_sampleID_to_Run(self, samplesheet: pd.DataFrame) -> dict:
        sampleID_to_Run = {sample: run for sample, run in zip(samplesheet['sample_id'], samplesheet['run_id'])}
        return sampleID_to_Run