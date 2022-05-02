#!/usr/bin/env snakemake

import os
import gzip
import datetime
import glob
import pysam

configfile: "config/cappseq_umi_config.yaml"


# Load samples
samplelist = []
with open(config["cappseq_umi_workflow"]["samplefile"]) as f:
    for line in f:
        line = line.rstrip("\n").rstrip("\r")
        samplelist.append(line)

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

    # Parse out the ttributes for the read group from the read name
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


rule trim_umi:
    input:
        r1 = lambda w: glob.glob(config["cappseq_umi_workflow"]["datadir"] + os.sep + w.samplename + "*R1_001.fastq.gz"),
        r2 = lambda w: glob.glob(config["cappseq_umi_workflow"]["datadir"] + os.sep + w.samplename + "*R2_001.fastq.gz")
    output:
        r1 = temp(outdir + os.sep + "01-trimmedfastqs" + os.sep + "{samplename}.r1.trimmed.fastq.gz"),
        r2 = temp(outdir + os.sep + "01-trimmedfastqs" + os.sep + "{samplename}.r2.trimmed.fastq.gz")
    params:
        barcodelength = config["cappseq_umi_workflow"]["barcodelength"]
    run:
        # Because globs are weird
        input.r1 = input.r1[0]
        input.r2 = input.r2[0]
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
        bam = temp(outdir + os.sep + "02-BWA" + os.sep + "{samplename}.bwa.unsort.bam")
    params:
        readgroup = lambda w: generate_read_group(glob.glob(config["cappseq_umi_workflow"]["datadir"] + os.sep + w.samplename + "*R1_001.fastq.gz")[0], w.samplename),
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
        r1 = lambda w: glob.glob(config["cappseq_umi_workflow"]["datadir"] + os.sep + w.samplename + "*R1_001.fastq.gz"),
        r2 = lambda w: glob.glob(config["cappseq_umi_workflow"]["datadir"] + os.sep + w.samplename + "*R2_001.fastq.gz"),
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
        bam = temp(outdir + os.sep + "04-umigrouped" + os.sep + "{samplename}.umigrouped.sort.bam")
    params:
        maxedits = config["cappseq_umi_workflow"]["umiedits"]
    threads: 
        config["cappseq_umi_workflow"]["samtools_sort_threads"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.groupumis.log"
    shell:
        "samtools sort -@ {threads} -n {input.bam} | fgbio SetMateInformation --ref {input.refgenome} | fgbio GroupReadsByUmi --edits {params.maxedits} --strategy paired > {output.bam} 2> {log}"

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
        readgroup = lambda w: generate_read_group(glob.glob(config["cappseq_umi_workflow"]["datadir"] + os.sep + w.samplename + "*R1_001.fastq.gz")[0], w.samplename)
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
        txt = outdir + os.sep + "09-validoutput" + os.sep + "{samplename}.consensus.mapped.annot.bam.is_valid"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        outdir + os.sep + "logs" + os.sep + "{samplename}.picardvalidatesam.log"
    shell:
        "picard ValidateSamFile -I {input.bam} -R {input.refgenome} > {output.txt} 2> {log}"

# Calculate duplicate rate
rule samtools_calculate_duplicate_rate:
    input:
        collapsed_bam = rules.fgbio_duplex_consensus.output.bam,
        all_reads_bam = rules.fgbio_group_umis.output.bam
    output:
        txt = outdir + os.sep + "08-duplicaterate" + os.sep + "{samplename}.duplicaterate.txt"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    shell:
        'seq 1 | awk -v totalreads=$(samtools flagstat {input.all_reads_bam} | grep read1 | cut -d " " -f 1) -v collapsedreads=$(samtools flagstat {input.collapsed_bam} | grep read1 | cut -d " " -f 1) \'BEGIN{{print((totalreads - collapsedreads)/totalreads)}}\' > {output.txt}'

rule all:
    input:
        expand([str(rules.picard_annotate_bam.output.bam),
                str(rules.samtools_calculate_duplicate_rate.output.txt),
                str(rules.picard_validate_sam.output.txt)], 
            samplename= samplelist)
