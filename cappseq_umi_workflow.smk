#!/usr/bin/env snakemake

import os
import gzip
import math
import sys

# configpath = "config/cappseq_umi_config.yaml"
# configfile: configpath

sys.path.append(config["cappseq_umi_workflow"]["repo_dir"])
import utils.path_utils as pu
import utils.fastq_utils as fu

# Ensure config file is correct, and all required attributes are present
# pathkeys = {"refgenome", "captureregions", "captureregionsil", "dbsnpvcf", "samplefile", "baseoutdir"}  # config elements that are filepaths
# for ckey, attribute in config["cappseq_umi_workflow"].items():
#     if attribute == "__UPDATE__":
#         raise AttributeError(f"\'__UPDATE__\' found for \'{ckey}\' in config file \'{configpath}\'. Please ensure the config file is updated with parameters relevant for your analysis")
#     # Check that required filepaths exist
#     if ckey in pathkeys:
#         if not os.path.exists(attribute):
#             raise AttributeError(f"Unable to locate \'{attribute}\' for key \'{ckey}\' in config file \'{configpath}\': No such file or directory")

# Load samples
# This assumes the sample file has 4 columns:
# sample_id    R1_fastq_path    R2_fastq_path    run_id
# initialize sample info object
SAMPLEINFO = pu.FastQInfo(config["cappseq_umi_workflow"]["samplefile"], config["cappseq_umi_workflow"]["fastq_dirs"])

SAMPLELIST = SAMPLEINFO.sample_ids
r1_fastqs = SAMPLEINFO.R1_fastqs
r2_fastqs = SAMPLEINFO.R2_fastqs
sample_to_runid = SAMPLEINFO.sampleID_to_Run

# Process sample file
# with open(samplefile) as f:
#     i = 0  # Line counter
#     for line in f:
#         i += 1
#         line = line.rstrip("\n").rstrip("\r")
#         if line.startswith("#"):  # Ignore comment lines
#             continue
#         try:
#             cols = line.split("\t")
#             sample_id = cols[0]
#             r1_fastq = cols[1]
#             r2_fastq = cols[2]
#             run_id = cols[3]
#         except IndexError as e:
#             raise AttributeError(f"Unable to parse line {i} of sample file \'{samplefile}\'. Expecting three tab-delineated columns corresponding to \'sample_id\', \'R1_fastq\', \'R2_fastq\', and \'run_id\'") from e
#         # Check that the specified FASTQs exist
#         if not os.path.exists(r1_fastq):
#             raise AttributeError(f"Unable to parse line {i} of sample file \'{samplefile}\': Unable to locate \'{r1_fastq}\': No such file or directory")
#         if not os.path.exists(r2_fastq):
#             raise AttributeError(f"Unable to parse line {i} of sample file \'{samplefile}\': Unable to locate \'{r2_fastq}\': No such file or directory")

#         # Check that this sample doesn't already exist
#         if sample_id in r1_fastqs:
#             raise AttributeError(f"Duplicate sample ID \'{sample_id}\' in sample file \'{samplefile}")
#         samplelist.append(sample_id)
#         r1_fastqs[sample_id] = r1_fastq
#         r2_fastqs[sample_id] = r2_fastq

#         # Store the runID for this sample
#         sample_to_runid[sample_id] = run_id

outdir = config["cappseq_umi_workflow"]["baseoutdir"]

rule trim_umi:
    """
    Trim the UMI from the front (and end) of each read pair
    Also does some FASTQ cleanup (polyG tails etc)

    Note that we are using --trim_front here instead of --umi, as we do not want the UMI stored in the read name
    as then we can't use fgbio AnnotateBamWithUMIs, as the read names will be different from the FASTQ file
    (In theory we could use CopyUmiFromReadName, but that requires a "-" deliminator between the forward and reverse
    UMIs while fastp uses a "_" so)
    """
    input:
        r1 = lambda w: r1_fastqs[w.samplename],
        r2 = lambda w: r2_fastqs[w.samplename]
    output:
        r1 = temp(os.path.join(outdir, "01-trimmedfastqs", "{samplename}.R1.trimmed.fastq.gz")),
        r2 = temp(os.path.join(outdir, "01-trimmedfastqs", "{samplename}.R2.trimmed.fastq.gz")),
        fastp_report = os.path.join(outdir, "01-trimmedfastqs", "{samplename}.fastp.json")
    params:
        barcodelength = config["cappseq_umi_workflow"]["barcodelength"],
        outdir = os.path.join(outdir, "01-trimmedfastqs")
    threads: 4
    conda:
        "envs/fastp.yaml"
    log:
        os.path.join(outdir, "logs", "{samplename}.fastp.log")
    shell: """
fastp --overrepresentation_analysis --detect_adapter_for_pe --trim_front1 {params.barcodelength} \
--trim_front2 {params.barcodelength} --in1 {input.r1} --in2 {input.r2} --thread {threads} --out1 {output.r1} --out2 {output.r2} \
--report_title "fastp {wildcards.samplename}" --json {output.fastp_report} --trim_poly_x \
--qualified_quality_phred 20 2> {log}
"""


rule bwa_align_unsorted:
    input:
        r1 = rules.trim_umi.output.r1,
        r2 = rules.trim_umi.output.r2,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = temp(os.path.join(outdir , "02-BWA" , "{samplename}.bwa.unsort.bam"))
    params:
        readgroup = lambda w: fu.generate_read_group(r1_fastqs[w.samplename], w.samplename, config),
    threads:
        config["cappseq_umi_workflow"]["bwa_threads"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir, "logs" , "{samplename}.bwa_allreads.log")
    shell:
        "bwa mem -t {threads} -R \"{params.readgroup}\" {input.refgenome} {input.r1} {input.r2} 2> {log} | samtools view -b > {output.bam} 2>> {log}"


# Add UMI tag
rule fgbio_annotate_umis:
    input:
        bam = rules.bwa_align_unsorted.output.bam,
        r1 = lambda w: r1_fastqs[w.samplename],
        r2 = lambda w: r2_fastqs[w.samplename]
    output:
        bam = temp(os.path.join(outdir , "03-withumis" , "{samplename}.bwa.umi.namesort.bam"))
    params:
        umiloc = config["cappseq_umi_workflow"]["barcodelocation"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.annotateumis.log")
    shell:
        "fgbio AnnotateBamWithUmis --input {input.bam} --fastq {input.r1} --fastq {input.r2} --read-structure {params.umiloc} --output {output.bam} &> {log}"

# Group reads by UMI into families
rule fgbio_group_umis:
    input:
        bam = rules.fgbio_annotate_umis.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = os.path.join(outdir, "04-umigrouped", "{samplename}.umigrouped.sort.bam"),
        txt = os.path.join(outdir , "04-umigrouped" , "{samplename}.umigrouped.famsize.txt")
    params:
        maxedits = config["cappseq_umi_workflow"]["umiedits"],
        outdir = os.path.join(outdir ,"04-umigrouped")
    threads:
        config["cappseq_umi_workflow"]["samtools_sort_threads"]
    resources:
        mem_mb = "2G"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs", "{samplename}.groupumis.log")
    shell:
        "samtools sort -m {resources.mem_mb} -@ {threads} -n {input.bam} | fgbio SetMateInformation --ref {input.refgenome} 2> {log} | fgbio GroupReadsByUmi --edits {params.maxedits} --family-size-histogram {output.txt} --strategy paired > {output.bam} 2>> {log}"

# Generate a consensus of these families
rule fgbio_duplex_consensus:
    input:
        bam = rules.fgbio_group_umis.output.bam
    output:
        bam = temp(os.path.join(outdir , "05-duplexconsensus" , "{samplename}.consensus.unmapped.bam"))
    params:
        minreads = config["cappseq_umi_workflow"]["minreads"],
        sampleid = "{samplename}"
    threads:
        config["cappseq_umi_workflow"]["duplexconsensus_threads"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.duplexconsensus.log")
    shell:
        "fgbio CallDuplexConsensusReads --input {input.bam} --output {output.bam} --threads {threads} --read-group-id {params.sampleid} --min-reads {params.minreads} &> {log}"


# Because CallDuplexConsensusReads recalculates base qualities, and those numbers can be above those supported by certain tools, set a upper
# limit to base qualities.
# Also, remove all the extra space-taking tags
rule sanitize_bam:
    input:
        bam = rules.fgbio_duplex_consensus.output.bam
    output:
        bam = temp(os.path.join(outdir , "06-sanitizebam" , "{samplename}.consensus.unmapped.capqual.bam"))
    params:
        max_base_qual = int(config["cappseq_umi_workflow"]["max_base_qual"]),  # Bases with quality scores above this are capped at this
        tagstoremove = config["cappseq_umi_workflow"]["tags_to_remove"],
        min_base_qual = int(config["cappseq_umi_workflow"]["min_base_qual"])  # Bases with quality scores below this are masked
    conda:
        "envs/pysam.yaml"
    threads:
        config["cappseq_umi_workflow"]["basequal_threads"]
    shell:
        """python utils/sanitize_bam.py --in_bam {input.bam} --out_bam {output.bam}
        --tagstoremove {params.tagstoremove} --max_base_qual {params.max_base_qual} 
        --threads {threads}"""

# rule sanitize_bam:
#     input:
#         bam = rules.fgbio_duplex_consensus.output.bam
#     output:
#         bam = temp(outdir + os.sep + "06-sanitizebam" + os.sep + "{samplename}.consensus.unmapped.capqual.bam")
#     params:
#         max_base_qual = int(config["cappseq_umi_workflow"]["max_base_qual"]),  # Bases with quality scores above this are capped at this
#         tagstoremove = config["cappseq_umi_workflow"]["tags_to_remove"],
#         min_base_qual = int(config["cappseq_umi_workflow"]["min_base_qual"])  # Bases with quality scores below this are masked
#     conda:
#         "envs/pysam.yaml"
#     threads:
#         config["cappseq_umi_workflow"]["basequal_threads"]
#     run:
#         inFile = pysam.AlignmentFile(input.bam, check_sq=False, mode = "rb", threads=threads)  # We have to provide check_sq=False in case this is an unaligned BAM
#         tagstoremove = set(params.tagstoremove)
#         # Remove the entries for these tags from the BAM header
#         inHeader = inFile.header.to_dict()  # Returns a multi-level dictionary
#         outHeader = {}
#         for level, attributes in inHeader.items():
#             # If we are removing the RG tag, remove this level
#             if "RG" in tagstoremove and level == "RG":
#                 continue
#             outHeader[level] = attributes

#         outFile = pysam.AlignmentFile(output.bam, header=outHeader, mode = "wb")

#         # Process reads
#         for read in inFile.fetch(until_eof=True):
#             # Cap the upper limit of base qualities
#             outqual = list(qual if qual <= params.max_base_qual else params.max_base_qual for qual in read.query_qualities)
#             # Mask bases with quality scores below the specified threshold
#             masked_seq = "".join(read.query_sequence[i] if read.query_qualities[i] >= 20 else "N" for i in range(0, read.query_length))
#             read.query_sequence = masked_seq
#             read.query_qualities = outqual

#             # Remove the unwanted tags
#             # pysam returns read tags as a list of tuples
#             # ex. [(NM, 2), (RG, "GJP00TM04")]
#             outtags = list(tag for tag in read.get_tags() if tag[0] not in tagstoremove)
#             read.set_tags(outtags)

#             outFile.write(read)

#         # For safety, manually close files
#         inFile.close()
#         outFile.close()


# Covert unaligned BAM back to FASTQ and map reads
rule bwa_realign_bam:
    input:
        bam = rules.sanitize_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = temp(os.path.join(outdir, "07-consensus_aligned", "{samplename}.consensus.mapped.namesort.bam"))
    threads:
        config["cappseq_umi_workflow"]["bwa_threads"]
    params:
        readgroup = lambda w: fu.generate_read_group(r1_fastqs[w.samplename], w.samplename, config)
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs", "{samplename}.bwa_realign.log")
    shell:
        "samtools fastq -@ {threads} -N {input.bam} 2> {log} | bwa mem -R \"{params.readgroup}\" -p -t {threads} {input.refgenome} - 2>> {log} | samtools view -b | "
        "picard SortSam -SO queryname -I /dev/stdin -O /dev/stdout 2>> {log} | fgbio SetMateInformation --ref {input.refgenome} --output {output.bam} &> {log}"

# Add back in family information from the unaligned consensus BAM
rule picard_annotate_bam:
    input:
        unaligned_bam = rules.sanitize_bam.output.bam,
        aligned_bam = rules.bwa_realign_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = os.path.join(outdir, "99-final", "{samplename}.consensus.mapped.annot.bam")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.picardannotate.log")
    shell:
        "picard SortSam -I {input.unaligned_bam} -O /dev/stdout -SO queryname --REFERENCE_SEQUENCE {input.refgenome} | picard MergeBamAlignment --ALIGNED_BAM {input.aligned_bam} --UNMAPPED_BAM /dev/stdin --REFERENCE_SEQUENCE {input.refgenome} -O {output.bam} &> {log} && "
        "samtools index {output.bam}"

# Output sentinel confirming that the final BAMs are valid
rule picard_validate_sam:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = os.path.join(outdir , "09-validoutput" , "{samplename}.consensus.mapped.ValidateSamFile.is_valid")
    params:
        outdir = os.path.join(outdir , "09-validoutput")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir, "logs" , "{samplename}.picardvalidatesam.log")
    shell:
        "picard ValidateSamFile -I {input.bam} -R {input.refgenome} > {output.txt} 2> {log}"


### QUALITY CONTROL RULES ###
### Currently supported tools:
### FASTQC
### PICARD

rule qc_fastqc:
    input:
        bam = rules.bwa_align_unsorted.output.bam
    output:
        qc = os.path.join(outdir , "Q1-fastqc" , "{samplename}.bwa.unsort_fastqc.html")
    params:
        outdir = os.path.join(outdir ,"Q1-fastqc")
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.fastqc.log")
    shell:
        "fastqc -o {params.outdir} --nogroup -f bam {input.bam} 2> {log}"

rule qc_picard_hsmetrics:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        hsmet = os.path.join(outdir , "Q2-hs_metrics" , "{samplename}.hs_metrics.txt"),
        tarcov = os.path.join(outdir , "Q2-hs_metrics" , "{samplename}.target_coverage.txt")
    params:
        capture_reg_il = config["cappseq_umi_workflow"]["captureregionsil"],
        outdir = outdir + os.sep + "Q2-hs_metrics",
        max_ram_records = "5000000",
        cov_cap_sens = "20000"
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.picard_hsmet.log")
    shell:
        "picard CollectHsMetrics -R {input.refgenome} -TI {params.capture_reg_il} -BI {params.capture_reg_il} -I {input.bam} -O {output.hsmet} --PER_TARGET_COVERAGE {output.tarcov} --MAX_RECORDS_IN_RAM {params.max_ram_records} --COVERAGE_CAP {params.cov_cap_sens} 2> {log}"

rule qc_picard_oxog:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = os.path.join(outdir, "Q3-oxog_metrics", "{samplename}.oxoG_metrics.txt")
    params:
        outdir = os.path.join(outdir , "Q3-oxog_metrics")
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(outdir , "logs", "{samplename}.picard_oxoG.log")
    shell:
        "picard CollectOxoGMetrics -I {input.bam} -R {input.refgenome} -O {output.txt} 2> {log}"

rule qc_picard_insertsize:
    input:
        bam = rules.picard_annotate_bam.output.bam
    output:
        txt = os.path.join(outdir, "Q4-insert_size", "{samplename}.insert_size_metrics.txt"),
        pdf = os.path.join(outdir , "Q4-insert_size" , "{samplename}.insert_size_histogram.pdf")
    params:
        outdir = os.path.join(outdir , "Q4-insert_size")
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.picard_insertsize.log")
    shell:
        "picard CollectInsertSizeMetrics -I {input.bam} -O {output.txt} -H {output.pdf} 2> {log}"

rule qc_fgbio_errorrate:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = os.path.join(outdir , "Q5-error_rate" , "{samplename}.error_rate_by_read_position.txt")
    params:
        outprefix = outdir + os.sep + "Q5-error_rate" + os.sep + "{samplename}",
        outdir = outdir + os.sep + "Q5-error_rate",
        capture_reg_il = config["cappseq_umi_workflow"]["captureregionsil"],
        dbsnpvcf = config["cappseq_umi_workflow"]["dbsnpvcf"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.error_rate_by_position.log")
    shell:
        "fgbio ErrorRateByReadPosition -i {input.bam} -r {input.refgenome} -v {params.dbsnpvcf} -l {params.capture_reg_il} --collapse -o {params.outprefix} 2> {log}"

rule qc_calc_dupl:
    input:
        collapsed_bam = rules.picard_annotate_bam.output.bam,
        all_reads_bam = rules.bwa_align_unsorted.output.bam
    output:
        txt = os.path.join(outdir , "Q6-dupl_rate" , "{samplename}.dup_metrics.txt")
    params:
        outdir = os.path.join(outdir, "Q6-dupl_rate")
    conda:
        "envs/pysam.yaml"
    shell:
        """python utils/qc_calc_dupl.py --collapsed_bam {input.collapsed_bam} --all_reads_bam {input.all_reads_bam}
        --outdir {params.outdir} --samplename {wildcards.samplename}"""

# rule qc_calc_dupl:
#     input:
#         collapsed_bam = rules.picard_annotate_bam.output.bam,
#         all_reads_bam = rules.bwa_align_unsorted.output.bam
#     output:
#         txt = outdir + os.sep + "Q6-dupl_rate" + os.sep + "{samplename}.dup_metrics.txt"
#     params:
#         outdir = outdir + os.sep + "Q6-dupl_rate"
#     conda:
#         "envs/pysam.yaml"
#     run:
#         # This definitely isn't the most efficient way to calculate duplicate rate, but it works
#         # We are going to generate a Picard-style output file
#         # First, parse the final (collapsed) BAM, and calculate the total number of reads
#         col_bam = pysam.AlignmentFile(input.collapsed_bam)
#         consensus_reads = col_bam.count(until_eof=True, read_callback=lambda x: not x.is_duplicate and x.is_mapped and x.is_paired and not x.is_supplementary and not x.is_secondary)
#         col_read_pairs = consensus_reads / 2
#         col_unpaired_reads = col_bam.count(until_eof=True, read_callback=lambda x: not x.is_paired and not x.is_supplementary and not x.is_secondary)

#         # Now, parse the original (non-consensus) BAM, and calculate the total number
#         # of read pairs, unmapped reads, and unpaired reads
#         orig_bam = pysam.AlignmentFile(input.all_reads_bam, require_index=False)
#         orig_reads = orig_bam.count(until_eof=True, read_callback=lambda x: not x.is_duplicate and x.is_mapped and x.is_paired and not x.is_supplementary and not x.is_secondary)
#         orig_pairs = orig_reads / 2
#         unmapped_reads = orig_bam.count(until_eof=True, read_callback=lambda x: x.is_unmapped and not x.is_supplementary and not x.is_secondary)
#         unpaired_reads = orig_bam.count(until_eof=True, read_callback=lambda x: not x.is_paired and not x.is_supplementary and not x.is_secondary)
#         secondary_reads = orig_bam.count(until_eof=True, read_callback=lambda x: x.is_supplementary or x.is_secondary)

#         # Now, calculate the output stats
#         library = wildcards.samplename
#         read_pairs_examined = int(orig_pairs)
#         secondary_or_supplemental = secondary_reads
#         unpaired_dups = int(unpaired_reads - col_unpaired_reads)
#         read_pair_dups = int(orig_pairs - col_read_pairs)
#         optical_dup = 0  # In this consensus approach (via UMIs) we lose info on which are optical duplicates
#         per_dupl = read_pair_dups / orig_pairs
#         # Custom added by Chris. Calculate the total size of this library, and how much we have sequenced
#         estimated_library_size = int(fu.estimateLibrarySize(orig_pairs, col_read_pairs))  # Use MultiQC's function to calculate the library size
#         prop_library_seq = col_read_pairs / estimated_library_size

#         # Close input
#         col_bam.close()
#         orig_bam.close()

#         # Write output
#         with open(output.txt, "w") as o:
#             # To trick MultiQC into thinking this is a Picard output
#             o.write("##picard.sam.markduplicates.MarkDuplicates BUT NOT ACTUALLY PICARD")
#             o.write(os.linesep)
#             header = ["LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
#                       "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE", "PERCENT_LIBRARY_SEQUENCED"]
#             o.write("\t".join(header))
#             o.write(os.linesep)
#             # Write out stats for this sample
#             out_values = [library, str(unpaired_reads), str(read_pairs_examined), str(secondary_or_supplemental), str(unmapped_reads),
#                           str(unpaired_dups), str(read_pair_dups), str(optical_dup), str(per_dupl), str(estimated_library_size), str(prop_library_seq)]
#             o.write("\t".join(out_values))
#             o.write(os.linesep)

# Output sentinel confirming that the final BAMs are valid
rule qc_validate_sam:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = os.path.join(outdir , "Q7-validatesam" , "{samplename}.consensus.mapped.ValidateSamFile.is_valid")
    params:
        outdir = os.path.join(outdir , "Q7-validatesam")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(outdir , "logs" , "{samplename}.picardvalidatesam.log")
    shell:
        "picard ValidateSamFile -I {input.bam} -R {input.refgenome} > {output.txt} 2> {log}"

# Merge QC results via multiqc
checkpoint qc_multiqc:
    input:
        # Run multiqc once per run, and merge all samples from that run
        dupl = lambda w: list(os.path.join(outdir, "Q6-dupl_rate", x + ".dup_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        errorate = lambda w: list(os.path.join(outdir, "Q5-error_rate", x + ".error_rate_by_read_position.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        insertsize = lambda w: list(os.path.join(outdir,"Q4-insert_size", x + ".insert_size_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        oxog = lambda w: list(os.path.join(outdir, "Q3-oxog_metrics", x + ".oxoG_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        hsmet = lambda w: list(os.path.join(outdir, "Q2-hs_metrics", x + ".hs_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        fastp = lambda w: list(os.path.join(outdir, "01-trimmedfastqs", x + ".fastp.json") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        fastqc = lambda w: list(os.path.join(outdir, "Q1-fastqc", x + ".bwa.unsort_fastqc.html") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        validatesam = lambda w: list(os.path.join(outdir, "Q7-validatesam", x + ".consensus.mapped.ValidateSamFile.is_valid") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        famsizehist = lambda w: list(os.path.join(outdir, "04-umigrouped", x + ".umigrouped.famsize.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid)
    output:
        html = os.path.join(outdir , "Q9-multiqc" , "multiqc_report.{runid}.html"),
    params:
        outdir = os.path.join(outdir , "Q9-multiqc"),
        outname = lambda w: "multiqc_report." + w.runid + ".html",
        ignoresamples = lambda w: "*\' --ignore-samples \'".join(x for x in SAMPLELIST if sample_to_runid[x] != w.runid),
        modules = "-m picard -m fastqc -m fgbio -m fastp",  # Should start with -m flag
        config = os.path.join(config["cappseq_umi_workflow"]["repo_dir"], "config/multiqc_config.yaml"),
        dupl_dir = rules.qc_calc_dupl.params.outdir,
        errorrate_dir = rules.qc_fgbio_errorrate.params.outdir,
        insertsize_dir = rules.qc_picard_insertsize.params.outdir,
        oxog_dir = rules.qc_picard_oxog.params.outdir,
        hsmet_dir = rules.qc_picard_hsmetrics.params.outdir,
        fastp_dir = rules.trim_umi.params.outdir,
        fastqc_dir = rules.qc_fastqc.params.outdir,
        validsam_dir = rules.qc_validate_sam.params.outdir,
        famsize_dir = rules.fgbio_group_umis.params.outdir
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(outdir, "logs" , "multiqc_{runid}.log")
    shell:
        "multiqc --no-data-dir --interactive --config {params.config} --outdir {params.outdir} --filename {params.outname} --force {params.modules} {params.dupl_dir} {params.errorrate_dir} {params.insertsize_dir} {params.oxog_dir} {params.hsmet_dir} {params.fastqc_dir} {params.validsam_dir} {params.famsize_dir} {params.fastp_dir} --ignore-samples \'{params.ignoresamples}*\' > {log}"

rule all:
    input:
        expand([str(rules.picard_annotate_bam.output.bam),
                str(rules.qc_multiqc.output.html)],
            samplename=SAMPLELIST,
            runid=set(sample_to_runid.values()))
    default_target: True

