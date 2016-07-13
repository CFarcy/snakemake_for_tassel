#Snakefile to launch Tassel 5 GBS pipeline v2.5.18

#Loading configuration from YAML file format
configfile: "config.yaml"

ruleorder:
    GBSSeqToTagDB > TagExportToFastq > bowtie2 > SAMToGBSdb > DiscoverySNPCaller > SNPQualityProfiler > ProductionSNPCaller > targets

rule targets:
    input:
        "01_TagDataBase/" + config['db_name'],
        "02_TagFastq/tagsForAlign.fa.gz",
        "03_SamFile/tagsForAlignFullvs.sam",
        "04_snpStatistics/snp_stats.txt",
        "05_snpCalling/" + config['vcf']



rule GBSSeqToTagDB:
    input:
        fastq_dir = config['fastq_dir'],
        keyfile = config['keyfile']

    output:
        database = "01_TagDataBase/" + config['db_name']

    params:
        java_memory = config['java_memory'][1],
        enzyme = config['enzyme'] ,
        kmer_length = config['GBSSeqToTagDBPlugin']['kmerLength'],
        minimum_kmer_length = config['GBSSeqToTagDBPlugin']['minKmerL'],
        batch_size = config['GBSSeqToTagDBPlugin']['batchSize'],
        minimum_quality_score = config['GBSSeqToTagDBPlugin']['mnQS'],
        maximum_number_of_kmer = config['GBSSeqToTagDBPlugin']['mxKmerNum'],
        min_kmer_count = config['GBSSeqToTagDBPlugin']['c'],




    shell:
        "run_pipeline.pl -Xmx{params.java_memory} -fork1 -GBSSeqToTagDBPlugin "
        "-i {input.fastq_dir} "
        "-k {input.keyfile} "
        "-db {output.database} "
        "-e {params.enzyme} "
        "-c {params.min_kmer_count} "
        "-mnQS {params.minimum_quality_score} "
        "-kmerLength {params.kmer_length} "
        "-minKmerL {params.minimum_kmer_length} "
        "-mxKmerNum {params.maximum_number_of_kmer} "
        "-batchSize {params.batch_size} "
        "-endPlugin -runfork1"


rule TagExportToFastq:
    input:
         database = rules.GBSSeqToTagDB.output

    output:
        fastq = "02_TagFastq/tagsForAlign.fa.gz"

    params:
        java_memory = config['java_memory'][0],
        min_count = config['TagExportToFastqPlugin']['c']


    shell:
        "run_pipeline.pl -Xmx{params.java_memory} -fork1 -TagExportToFastqPlugin "
        "-db {input.database} "
        "-o {output.fastq} "
        "-c {params.min_count} "
        "-endPlugin -runfork1"


rule bowtie2:
    input:
        fastq = rules.TagExportToFastq.output,
        genome = config['genome']

    output:
        sam = "03_SamFile/tagsForAlignFullvs.sam"


    shell:
        "bowtie2 -p 15 --very-sensitive -x {input.genome} -U {input.fastq} -S {output.sam}"

rule SAMToGBSdb :
    input:
        sam = rules.bowtie2.output,
        database = rules.GBSSeqToTagDB.output

    output:
        "03_SamFile/sam_to_gbs.out"

    params:
        java_memory = config['java_memory'][0],
        aLen = config['SAMToGBSdbPlugin']['aLen'],
        aProp = config['SAMToGBSdbPlugin']['aProp']


    shell:
        "run_pipeline.pl -Xmx{params.java_memory} -fork1 -SAMToGBSdbPlugin "
        "-i {input.sam} "
        "-db {input.database} "
        "-aLen {params.aLen} "
        "-aProp {params.aProp} "
        "-endPlugin -runfork1;"
        "echo 'SAMToGBSdb done!' > {output};"
        "touch 02_TagFastq/*"

rule DiscoverySNPCaller:
    input:
        database = rules.GBSSeqToTagDB.output,
        temp_file = rules.SAMToGBSdb.output

    output:
        "03_SamFile/snp_caller.out"

    params:
        java_memory = config['java_memory'][0],
        start_chromosome = config['DiscoverySNPCallerPluginV2']['sC'],
        stop_chromosome = config['DiscoverySNPCallerPluginV2']['eC'],
        min_minor_allel_freq = config['DiscoverySNPCallerPluginV2']['mnMAF'],
        min_locus_coverage = config['DiscoverySNPCallerPluginV2']['mnLCov'],
        max_tags_cut_site = config['DiscoverySNPCallerPluginV2']['maxTagsCutSite']



    shell:
        "run_pipeline.pl -Xmx{params.java_memory} -fork1 -DiscoverySNPCallerPluginV2 "
        "-db {input.database} "
        "-sC {params.start_chromosome} "
        "-eC {params.stop_chromosome} "
        "-mnLCov {params.min_locus_coverage} "
        "-mnMAF {params.min_minor_allel_freq} "
        "-maxTagsCutSite {params.max_tags_cut_site} "
        "-endPlugin -runfork1;"
        "echo 'DiscoverySNPCaller done!' > {output};"
        "touch 02_TagFastq/*; touch 03_SamFile/*"

rule SNPQualityProfiler:
    input:
        database = rules.GBSSeqToTagDB.output,
        temp_file = rules.DiscoverySNPCaller.output

    output:
        stat_file = "04_snpStatistics/snp_stats.txt"

    params:
        java_memory = config['java_memory'][0],
        tname = config['SNPQualityProfilerPlugin']['tname']


    shell:
        "run_pipeline.pl -debug -Xmx{params.java_memory} -fork1 -SNPQualityProfilerPlugin "
        "-db {input.database} "
        "-tname  {params.tname} "
        "-statFile {output.stat_file} "
        "-endPlugin -runfork1"

rule ProductionSNPCaller:
    input:
        keyfile = config['keyfile'],
        fastq_dir = config['fastq_dir'],
        stat_file = "04_snpStatistics/snp_stats.txt"

    output:
        vcf = "05_snpCalling/" + config['vcf']

    params:
        database = rules.GBSSeqToTagDB.output,
        java_memory = config['java_memory'][1],
        enzyme = config['enzyme'],
        kmer_length = config['ProductionSNPCallerPluginV2']['kmerLength']


    shell:
        "run_pipeline.pl -Xmx{params.java_memory} -fork1 -ProductionSNPCallerPluginV2 "
        "-db {params.database} "
        "-i {input.fastq_dir} "
        "-k {input.keyfile} "
        "-e {params.enzyme} "
        "-o {output.vcf} "
        "-kmerLength {params.kmer_length} "
        "-endPlugin -runfork1"
