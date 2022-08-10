#!/usr/bin/env snakemake

import os
import gzip
import datetime
import pysam
import math

configpath = "config/cappseq_umi_config.yaml"
configfile: configpath

# Ensure config file is correct, and all required attributes are present
pathkeys = {"refgenome", "captureregions", "captureregionsil", "dbsnpvcf", "samplefile", "baseoutdir"}  # config elements that are filepaths
for ckey, attribute in config["cappseq_umi_workflow"].items():
    if attribute == "__UPDATE__":
        raise AttributeError(f"\'__UPDATE__\' found for \'{ckey}\' in config file \'{configpath}\'. Please ensure the config file is updated with parameters relevant for your analysis")
    # Check that required filepaths exist
    if ckey in pathkeys:
        if not os.path.exists(attribute):
            raise AttributeError(f"Unable to locate \'{attribute}\' for key \'{ckey}\' in config file \'{configpath}\': No such file or directory")

# Load samples
# This assumes the sample file has three columns:
# sample_id    R1_fastq_path    R2_fastq_path
samplefile = config["cappseq_umi_workflow"]["samplefile"]
samplelist = []
r1_fastqs = {}
r2_fastqs = {}

with open(samplefile) as f:
    i = 0  # Line counter
    for line in f:
        i += 1
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith("#"):  # Ignore comment lines
            continue
        try:
            cols = line.split("\t")
            sample_id = cols[0]
            r1_fastq = cols[1]
            r2_fastq = cols[2]
        except IndexError as e:
            raise AttributeError(f"Unable to parse line {i} of sample file \'{samplefile}\'. Expecting three tab-delineated columns corresponding to \'sample_id\', \'R1_fastq\', and \'R2_fastq\'") from e
        # Check that the specified FASTQs exist
        if not os.path.exists(r1_fastq):
            raise AttributeError(f"Unable to parse line {i} of sample file \'{samplefile}\': Unable to locate \'{r1_fastq}\': No such file or directory")
        if not os.path.exists(r2_fastq):
            raise AttributeError(f"Unable to parse line {i} of sample file \'{samplefile}\': Unable to locate \'{r2_fastq}\': No such file or directory")

        # Check that this sample doesn't already exist
        if sample_id in r1_fastqs:
            raise AttributeError(f"Duplicate sample ID \'{sample_id}\' in sample file \'{samplefile}")
        samplelist.append(sample_id)
        r1_fastqs[sample_id] = r1_fastq
        r2_fastqs[sample_id] = r2_fastq

outdir = config["cappseq_umi_workflow"]["baseoutdir"]


def is_gzipped(filepath):
    with open(filepath, "rb") as f:
        magicnum = f.read(2)
        return magicnum == b'\x1f\x8b'  # Magic number of gzipped files


def generate_read_group(fastq, sample):
    # Parses flowcell, lane, and barcode information from FASTQ read names
    # Uses this information (and config file info) to generate a read group line
    if is_gzipped(fastq):
        readname = gzip.open(fastq, "rt").readline()
    else:
        readname = open(fastq, "r").readline()

    # Parse out the attributes for the read group from the read name
    readname = readname.rstrip("\n").rstrip("\r")
    cols = readname.split(":")

    sequencer = cols[0]
    sequencer = sequencer.replace("@","")
    flowcell = cols[2]
    flowcell = flowcell.split("-")[-1]
    lane = cols[3]
    barcode = cols[-1]

    # From the config (generic and should be consistent between runs)
    description = config["cappseq_umi_workflow"]["readgroup"]["description"]
    centre = config["cappseq_umi_workflow"]["readgroup"]["centre"]
    platform = config["cappseq_umi_workflow"]["readgroup"]["platformunit"]
    platformmodel = config["cappseq_umi_workflow"]["readgroup"]["platformmodel"]

    date = datetime.datetime.now().date().isoformat()  # I KNOW ITS UGLY BUT I JUST WANT THE DATE IN ISO FORMAT
    platformunit = flowcell + "-" + lane + ":" + barcode

    readgroup = f"@RG\\tID:{sample}\\tBC:{barcode}\\tCN:{centre}\\tDS:\'{description}\'\\tDT:{date}\\tLB:{sample}\\tPL:{platform}\\tPM:{platformmodel}\\tPU:{platformunit}\\tSM:{sample}"
    return readgroup


# Adapted from the multiQC module
# Used to calculate the estimated total number of molecules in a library
def estimateLibrarySize(readPairs, uniqueReadPairs):
    """
    Picard calculation to estimate library size
    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L153-L164
    Note: Optical duplicates are contained in duplicates and therefore do not enter the calculation here.
    See also the computation of READ_PAIR_NOT_OPTICAL_DUPLICATES.
     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs observed.
     * <p>
     * Based on the Lander-Waterman equation that states:
     * C/X = 1 - exp( -N/X )
     * where
     * X = number of distinct molecules in library
     * N = number of read pairs
     * C = number of distinct fragments observed in read pairs
    """

    readPairDuplicates = readPairs - uniqueReadPairs

    if readPairs > 0 and readPairDuplicates > 0:

        m = 1.0
        M = 100.0

        if uniqueReadPairs >= readPairs or f(m * uniqueReadPairs, uniqueReadPairs, readPairs) < 0:
            logging.warning("Picard recalculation of ESTIMATED_LIBRARY_SIZE skipped - metrics look wrong")
            return None

        # find value of M, large enough to act as other side for bisection method
        while f(M * uniqueReadPairs, uniqueReadPairs, readPairs) > 0:
            M *= 10.0

        # use bisection method (no more than 40 times) to find solution
        for i in range(40):
            r = (m + M) / 2.0
            u = f(r * uniqueReadPairs, uniqueReadPairs, readPairs)
            if u == 0:
                break
            elif u > 0:
                m = r
            elif u < 0:
                M = r

        return uniqueReadPairs * (m + M) / 2.0
    else:
        return None


def f(x, c, n):
    """
    Picard calculation used when estimating library size
    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L172-L177
    * Method that is used in the computation of estimated library size.
    """
    return c / x - 1 + math.exp(-n / x)


rule trim_umi:
    input:
        r1 = lambda w: r1_fastqs[w.samplename],
        r2 = lambda w: r2_fastqs[w.samplename]
    output:
        r1 = temp(outdir + os.sep + "01-trimmedfastqs" + os.sep + "{samplename}.r1.trimmed.fastq.gz"),
        r2 = temp(outdir + os.sep + "01-trimmedfastqs" + os.sep + "{samplename}.r2.trimmed.fastq.gz")
    params:
        barcodelength = config["cappseq_umi_workflow"]["barcodelength"]
    run:
        # Sanity check barcode
        barcodelength = int(params.barcodelength)
        if barcodelength <= 0:
            raise AttributeError("\'barcodelength\' must be greater than 0")
        # Trim read1
        open_func = lambda x: gzip.open(x, "rt") if is_gzipped(x) else open(x)
        write_func = lambda x: gzip.open(x, "wt") if x.endswith(".gz") else open(x, "w")
        with open_func(input.r1) as f:
            with write_func(output.r1) as o:
                # WARNING: WILL NOT WORK WITH READS SHORTER THAN 8BP. WILL BREAK.
                i = 0
                for line in f:
                    i += 1
                    # Trim the sequence and quality lines for each read
                    if i % 2 == 0:
                        line = line[barcodelength:]
                    o.write(line)

        # Trim read2
        with open_func(input.r2) as f:
            with write_func(output.r2) as o:
                # WARNING: WILL NOT WORK WITH READS SHORTER THAN 8BP. WILL BREAK.
                i = 0
                for line in f:
                    i += 1
                    # Trim the sequence and quality lines for each read
                    if i % 2 == 0:
                        line = line[barcodelength:]
                    o.write(line)


rule bwa_align_unsorted:
    input:
        r1 = rules.trim_umi.output.r1,
        r2 = rules.trim_umi.output.r2,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = outdir + os.sep + "02-BWA" + os.sep + "{samplename}.bwa.unsort.bam"
    params:
        readgroup = lambda w: generate_read_group(r1_fastqs[w.samplename], w.samplename),
    threads:
        config["cappseq_umi_workflow"]["bwa_threads"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.bwa_allreads.log"
    shell:
        "bwa mem -t {threads} -R \"{params.readgroup}\" {input.refgenome} {input.r1} {input.r2} 2> {log} | samtools view -b > {output.bam} 2> {log}"


# Add UMI tag
rule fgbio_annotate_umis:
    input:
        bam = rules.bwa_align_unsorted.output.bam,
        r1 = lambda w: r1_fastqs[w.samplename],
        r2 = lambda w: r2_fastqs[w.samplename]
    output:
        bam = temp(outdir + os.sep + "03-withumis" + os.sep + "{samplename}.bwa.umi.namesort.bam")
    params:
        umiloc = config["cappseq_umi_workflow"]["barcodelocation"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.annotateumis.log"
    shell:
        "fgbio AnnotateBamWithUmis --input {input.bam} --fastq {input.r1} --fastq {input.r2} --sorted true --read-structure {params.umiloc} --output {output.bam} &> {log}"

# Group reads by UMI into families
rule fgbio_group_umis:
    input:
        bam = rules.fgbio_annotate_umis.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = temp(outdir + os.sep + "04-umigrouped" + os.sep + "{samplename}.umigrouped.sort.bam"),
        txt = outdir + os.sep + "04-umigrouped" + os.sep + "{samplename}.umigrouped.famsize.txt"
    params:
        maxedits = config["cappseq_umi_workflow"]["umiedits"],
        outdir = outdir + os.sep + "04-umigrouped"
    threads:
        config["cappseq_umi_workflow"]["samtools_sort_threads"]
    resources:
        mem_mb = "2G"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.groupumis.log"
    shell:
        "samtools sort -m {resources.mem_mb} -@ {threads} -n {input.bam} | fgbio SetMateInformation --ref {input.refgenome} 2> {log} | fgbio GroupReadsByUmi --edits {params.maxedits} --family-size-histogram {output.txt} --strategy paired > {output.bam} 2>> {log}"

# Generate a consensus of these families
rule fgbio_duplex_consensus:
    input:
        bam = rules.fgbio_group_umis.output.bam
    output:
        bam = temp(outdir + os.sep + "05-duplexconsensus" + os.sep + "{samplename}.consensus.unmapped.bam")
    params:
        minreads = config["cappseq_umi_workflow"]["minreads"],
        sampleid = "{samplename}"
    threads:
        config["cappseq_umi_workflow"]["duplexconsensus_threads"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.duplexconsensus.log"
    shell:
        "fgbio CallDuplexConsensusReads --input {input.bam} --output {output.bam} --threads {threads} --read-group-id {params.sampleid} --min-reads {params.minreads} &> {log}"


# Because CallDuplexConsensusReads recalculates base qualities, and those numbers can be above those supported by certain tools, set a upper
# limit to base qualities.
# Also, remove all the extra space-taking tags
rule cap_base_qual:
    input:
        bam = rules.fgbio_duplex_consensus.output.bam
    output:
        bam = outdir + os.sep + "06-sanitizebam" + os.sep + "{samplename}.consensus.unmapped.capqual.bam"
    params:
        max_base_qual = int(config["cappseq_umi_workflow"]["max_base_qual"]),
        tagstoremove = config["cappseq_umi_workflow"]["tags_to_remove"]
    threads:
        config["cappseq_umi_workflow"]["basequal_threads"]
    run:
        inFile = pysam.AlignmentFile(input.bam, check_sq=False, mode = "rb", threads=threads)  # We have to provide check_sq=False in case this is an unaligned BAM
        tagstoremove = set(params.tagstoremove)
        # Remove the entries for these tags from the BAM header
        inHeader = inFile.header.to_dict()  # Returns a multi-level dictionary
        outHeader = {}
        for level, attributes in inHeader.items():
            # If we are removing the RG tag, remove this level
            if "RG" in tagstoremove and level == "RG":
                continue
            outHeader[level] = attributes

        outFile = pysam.AlignmentFile(output.bam, header=outHeader, mode = "wb")

        # Process reads
        for read in inFile.fetch(until_eof=True):
            # Cap the upper limit of base qualities
            outqual = list(qual if qual <= params.max_base_qual else params.max_base_qual for qual in read.query_qualities)
            read.query_qualities = outqual

            # Remove the unwanted tags
            # pysam returns read tags as a list of tuples
            # ex. [(NM, 2), (RG, "GJP00TM04")]
            outtags = list(tag for tag in read.get_tags() if tag[0] not in tagstoremove)
            read.set_tags(outtags)

            outFile.write(read)

        # For safety, manually close files
        inFile.close()
        outFile.close()


# Covert unaligned BAM back to FASTQ and map reads
rule bwa_realign_bam:
    input:
        bam = rules.cap_base_qual.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = temp(outdir + os.sep + "07-consensus_aligned" + os.sep + "{samplename}.consensus.mapped.namesort.bam")
    threads:
        config["cappseq_umi_workflow"]["bwa_threads"]
    params:
        readgroup = lambda w: generate_read_group(r1_fastqs[w.samplename], w.samplename)
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.bwa_realign.log"
    shell:
        "samtools fastq -@ {threads} -N {input.bam} | bwa mem -R \"{params.readgroup}\" -p -t {threads} {input.refgenome} - | samtools view -b | "
        "picard SortSam -SO queryname -I /dev/stdin -O /dev/stdout | fgbio SetMateInformation --ref {input.refgenome} --output {output.bam} &> {log}"

# Add back in family information from the unaligned consensus BAM
rule picard_annotate_bam:
    input:
        unaligned_bam = rules.cap_base_qual.output.bam,
        aligned_bam = rules.bwa_realign_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        bam = outdir + os.sep + "99-final" + os.sep + "{samplename}.consensus.mapped.annot.bam"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picardannotate.log"
    shell:
        "picard SortSam -I {input.unaligned_bam} -O /dev/stdout -SO queryname --REFERENCE_SEQUENCE {input.refgenome} | picard MergeBamAlignment --ALIGNED_BAM {input.aligned_bam} --UNMAPPED_BAM /dev/stdin --REFERENCE_SEQUENCE {input.refgenome} -O {output.bam} &> {log} && "
        "samtools index {output.bam}"

# Output sentinel confirming that the final BAMs are valid
rule picard_validate_sam:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = outdir + os.sep + "09-validoutput" + os.sep + "{samplename}.consensus.mapped.ValidateSamFile.is_valid"
    params:
        outdir = outdir + os.sep + "09-validoutput"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picardvalidatesam.log"
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
        qc = outdir + os.sep + "Q1-fastqc" + os.sep + "{samplename}.bwa.unsort_fastqc.html"
    params:
        outdir = outdir + os.sep + "Q1-fastqc"
    conda:
        "envs/picard_fastqc.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.fastqc.log"
    shell:
        "fastqc -o {params.outdir} --nogroup -f bam {input.bam} 2> {log}"

rule qc_picard_hsmetrics:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        hsmet = outdir + os.sep + "Q2-hs_metrics" + os.sep + "{samplename}.hs_metrics.txt",
        tarcov = outdir + os.sep + "Q2-hs_metrics" + os.sep + "{samplename}.target_coverage.txt"
    params:
        capture_reg_il = config["cappseq_umi_workflow"]["captureregionsil"],
        outdir = outdir + os.sep + "Q2-hs_metrics",
        max_ram_records = "5000000",
        cov_cap_sens = "20000"
    conda:
        "envs/picard_fastqc.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picard_hsmet.log"
    shell:
        "picard CollectHsMetrics -R {input.refgenome} -TI {params.capture_reg_il} -BI {params.capture_reg_il} -I {input.bam} -O {output.hsmet} --PER_TARGET_COVERAGE {output.tarcov} --MAX_RECORDS_IN_RAM {params.max_ram_records} --COVERAGE_CAP {params.cov_cap_sens} 2> {log}"

rule qc_picard_oxog:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = outdir + os.sep + "Q3-oxog_metrics" + os.sep + "{samplename}.oxoG_metrics.txt"
    params:
        outdir = outdir + os.sep + "Q3-oxog_metrics"
    conda:
        "envs/picard_fastqc.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picard_oxoG.log"
    shell:
        "picard CollectOxoGMetrics -I {input.bam} -R {input.refgenome} -O {output.txt} 2> {log}"

rule qc_picard_insertsize:
    input:
        bam = rules.picard_annotate_bam.output.bam
    output:
        txt = outdir + os.sep + "Q4-insert_size" + os.sep + "{samplename}.insert_size_metrics.txt",
        pdf = outdir + os.sep + "Q4-insert_size" + os.sep + "{samplename}.insert_size_histogram.pdf"
    params:
        outdir = outdir + os.sep + "Q4-insert_size"
    conda:
        "envs/picard_fastqc.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picard_insertsize.log"
    shell:
        "picard CollectInsertSizeMetrics -I {input.bam} -O {output.txt} -H {output.pdf} 2> {log}"

rule qc_fgbio_errorrate:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = outdir + os.sep + "Q5-error_rate" + os.sep + "{samplename}.error_rate_by_read_position.txt"
    params:
        outprefix = outdir + os.sep + "Q5-error_rate" + os.sep + "{samplename}",
        outdir = outdir + os.sep + "Q5-error_rate",
        capture_reg_il = config["cappseq_umi_workflow"]["captureregionsil"],
        dbsnpvcf = config["cappseq_umi_workflow"]["dbsnpvcf"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.error_rate_by_position.log"
    shell:
        "fgbio ErrorRateByReadPosition -i {input.bam} -r {input.refgenome} -v {params.dbsnpvcf} -l {params.capture_reg_il} --collapse -o {params.outprefix} 2> {log}"

rule qc_calc_dupl:
    input:
        collapsed_bam = rules.picard_annotate_bam.output.bam,
        all_reads_bam = rules.bwa_align_unsorted.output.bam
    output:
        txt = outdir + os.sep + "Q6-dupl_rate" + os.sep + "{samplename}.dup_metrics.txt"
    params:
        outdir = outdir + os.sep + "Q6-dupl_rate"
    run:
        # This definitely isn't the most efficient way to calculate duplicate rate, but it works
        # We are going to generate a Picard-style output file
        # First, parse the final (collapsed) BAM, and calculate the total number of reads
        col_bam = pysam.AlignmentFile(input.collapsed_bam)
        consensus_reads = col_bam.count(until_eof=True, read_callback=lambda x: not x.is_duplicate and x.is_mapped and x.is_paired and not x.is_supplementary and not x.is_secondary)
        col_read_pairs = consensus_reads / 2
        col_unpaired_reads = col_bam.count(until_eof=True, read_callback=lambda x: not x.is_paired and not x.is_supplementary and not x.is_secondary)

        # Now, parse the original (non-consensus) BAM, and calculate the total number
        # of read pairs, unmapped reads, and unpaired reads
        orig_bam = pysam.AlignmentFile(input.all_reads_bam, require_index=False)
        orig_reads = orig_bam.count(until_eof=True, read_callback=lambda x: not x.is_duplicate and x.is_mapped and x.is_paired and not x.is_supplementary and not x.is_secondary)
        orig_pairs = orig_reads / 2
        unmapped_reads = orig_bam.count(until_eof=True, read_callback=lambda x: x.is_unmapped and not x.is_supplementary and not x.is_secondary)
        unpaired_reads = orig_bam.count(until_eof=True, read_callback=lambda x: not x.is_paired and not x.is_supplementary and not x.is_secondary)
        secondary_reads = orig_bam.count(until_eof=True, read_callback=lambda x: x.is_supplementary or x.is_secondary)

        # Now, calculate the output stats
        library = wildcards.samplename
        read_pairs_examined = int(orig_pairs)
        secondary_or_supplemental = secondary_reads
        unpaired_dups = int(unpaired_reads - col_unpaired_reads)
        read_pair_dups = int(orig_pairs - col_read_pairs)
        optical_dup = 0  # In this consensus approach (via UMIs) we lose info on which are optical duplicates
        per_dupl = read_pair_dups / orig_pairs
        # Custom added by Chris. Calculate the total size of this library, and how much we have sequenced
        estimated_library_size = int(estimateLibrarySize(orig_pairs, col_read_pairs))  # Use MultiQC's function to calculate the library size
        prop_library_seq = col_read_pairs / estimated_library_size

        # Close input
        col_bam.close()
        orig_bam.close()

        # Write output
        with open(output.txt, "w") as o:
            # To trick MultiQC into thinking this is a Picard output
            o.write("##picard.sam.markduplicates.MarkDuplicates BUT NOT ACTUALLY PICARD")
            o.write(os.linesep)
            header = ["LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
                      "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE", "PERCENT_LIBRARY_SEQUENCED"]
            o.write("\t".join(header))
            o.write(os.linesep)
            # Write out stats for this sample
            out_values = [library, str(unpaired_reads), str(read_pairs_examined), str(secondary_or_supplemental), str(unmapped_reads),
                          str(unpaired_dups), str(read_pair_dups), str(optical_dup), str(per_dupl), str(estimated_library_size), str(prop_library_seq)]
            o.write("\t".join(out_values))
            o.write(os.linesep)

# Output sentinel confirming that the final BAMs are valid
rule qc_validate_sam:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["cappseq_umi_workflow"]["refgenome"]
    output:
        txt = outdir + os.sep + "Q7-validatesam" + os.sep + "{samplename}.consensus.mapped.ValidateSamFile.is_valid"
    params:
        outdir = outdir + os.sep + "Q7-validatesam"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picardvalidatesam.log"
    shell:
        "picard ValidateSamFile -I {input.bam} -R {input.refgenome} > {output.txt} 2> {log}"

# Merge QC results via multiqc
checkpoint qc_multiqc:
    input:
        dupl = expand(outdir + os.sep + "Q6-dupl_rate" + os.sep + "{samplename}.dup_metrics.txt", samplename=samplelist),
        errorate = expand(outdir + os.sep + "Q5-error_rate" + os.sep + "{samplename}.error_rate_by_read_position.txt", samplename=samplelist),
        insertsize = expand(outdir + os.sep + "Q4-insert_size" + os.sep + "{samplename}.insert_size_metrics.txt", samplename=samplelist),
        oxog = expand(outdir + os.sep + "Q3-oxog_metrics" + os.sep + "{samplename}.oxoG_metrics.txt", samplename=samplelist),
        hsmet = expand(outdir + os.sep + "Q2-hs_metrics" + os.sep + "{samplename}.hs_metrics.txt", samplename=samplelist),
        fastqc = expand(outdir + os.sep + "Q1-fastqc" + os.sep + "{samplename}.bwa.unsort_fastqc.html", samplename=samplelist),
        validatesam = expand(outdir + os.sep + "Q7-validatesam" + os.sep + "{samplename}.consensus.mapped.ValidateSamFile.is_valid", samplename=samplelist),
        famsizehist = expand(outdir + os.sep + "04-umigrouped" + os.sep + "{samplename}.umigrouped.famsize.txt", samplename=samplelist)
    output:
        html = outdir + os.sep + "Q9-multiqc" + os.sep + "multiqc_report.html",
    params:
        outdir = outdir + os.sep + "Q9-multiqc",
        modules = "-m picard -m fastqc -m fgbio",  # Should start with -m flag
        config = config["cappseq_umi_workflow"]["multiqc_config"],
        dupl_dir = rules.qc_calc_dupl.params.outdir,
        errorrate_dir = rules.qc_fgbio_errorrate.params.outdir,
        insertsize_dir = rules.qc_picard_insertsize.params.outdir,
        oxog_dir = rules.qc_picard_oxog.params.outdir,
        hsmet_dir = rules.qc_picard_hsmetrics.params.outdir,
        fastqc_dir = rules.qc_fastqc.params.outdir,
        validsam_dir = rules.qc_validate_sam.params.outdir,
        famsize_dir = rules.fgbio_group_umis.params.outdir
    conda:
        "envs/picard_fastqc.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "multiqc_all.log"
    shell:
        "multiqc --interactive --config {params.config} --outdir {params.outdir} --force {params.modules} {params.dupl_dir} {params.errorrate_dir} {params.insertsize_dir} {params.oxog_dir} {params.hsmet_dir} {params.fastqc_dir} {params.validsam_dir} {params.famsize_dir} > {log}"

rule all:
    input:
        expand([str(rules.picard_annotate_bam.output.bam),
                str(rules.qc_multiqc.output.html)],
            samplename= samplelist)
    default_target: True

