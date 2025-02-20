include { QualityControl } from './binning_qc'
include { MGEAnalysis    } from './mge_analysis'

process ProcessGFA {
    cpus 2
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(gfa_file)

    output:
    path("segments.fna"), emit: segments
    path("circular.ids"), emit: circular_ids

    """
    gfa_utils dump-segments $gfa_file segments.fna
    gfa_utils isolates --circular $gfa_file | sed 1d | cut -d, -f1 > circular.ids
    """
}

process MappabilityMasking {
    cpus 8
    memory '64 GB'
    conda params.conda.genmap
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(segments)

    output:
    path('genmap.*')
    path("segments.masked.fna"), emit: masked

    shell:
    '''
    genmap index -F !{segments} -I genmap.index && \
        genmap map -T !{task.cpus} -I genmap.index -O genmap -bg -K !{params.genmap.kmer} && \
        gawk -F '\t' 'BEGIN {OFS="\t"} {if ($4 < !{params.genmap.minscore}) print $1, $2, $3, ".", $4; }' genmap.bedgraph > genmap.bed && \
        bedtools maskfasta -fi !{segments} -fo segments.masked.fna -bed genmap.bed
    '''
}

process IndexAndMap {
    cpus 25
    memory '32 GB'
    publishDir params.outdir, mode: 'copy'
    scratch params.scratch_dir
    //conda params.conda.bin3c

    input:
    path(asm_fasta)
    path(hic_r1)
    path(hic_r2)

    output:
    path("hic2ctg.bam"), emit: bam_file

    """
    bwa index $asm_fasta
    bwa mem -5SP -t ${task.cpus} $asm_fasta $hic_r1 $hic_r2 | \
        samtools view -F 0x904 -@4 -uS | \
        samtools sort -@ ${task.cpus} -n -o "hic2ctg.bam"
    """
}

process MakeContactMap {
    cpus 4
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(asm_fasta)
    path(bam_file)

    output:
    path("map_out")

    script:
    // single enzyme digest
    if (params.bin3c.enz2 == null) {
        """
        bin3C mkmap --threads ${task.cpus} \
            --min-insert ${params.bin3c.min_insert} \
            --min-extent ${params.bin3c.min_extent} \
            --min-mapq ${params.bin3c.min_mapq} \
            --max-edist ${params.bin3c.max_edist} \
            --bin-size ${params.bin3c.bin_size} \
            --min-alen ${params.bin3c.min_alen} \
            -e ${params.bin3c.enz1} \
            $asm_fasta $bam_file map_out
        """
    }
    // double enzyme digest
    else {
        """
        bin3C mkmap --threads ${task.cpus} \
            --min-insert ${params.bin3c.min_insert} \
            --min-extent ${params.bin3c.min_extent} \
            --min-mapq ${params.bin3c.min_mapq} \
            --max-edist ${params.bin3c.max_edist} \
            --bin-size ${params.bin3c.bin_size} \
            --min-alen ${params.bin3c.min_alen} \
            -e ${params.bin3c.enz1} -e ${params.bin3c.enz2} \
            $asm_fasta $bam_file map_out
        """
    }
}

process ClusterMetagenome {
    cpus 2
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(asm_fasta)
    path(map_out)
    path(excludes)

    output:
    path("cluster_out")

    """
    bin3C cluster --verbose -s ${params.bin3c.seed} \
        --markov-scale ${params.bin3c.markov_scale} \
        --min-reflen ${params.bin3c.min_reflen} \
        --assembler ${params.bin3c.assembler} \
        --min-sig ${params.bin3c.min_sig} \
        --fdr-alpha ${params.bin3c.fdr_alpha} \
        --from-extent --norm-method ${params.bin3c.norm_method} \
        --n-iter ${params.bin3c.n_iter} \
        --min-extent ${params.bin3c.min_extent} \
        --plot-contrast ${params.bin3c.plot_contrast} \
        --exclude-from $excludes \
        --fasta $asm_fasta \
        "${map_out}/contact_map.p.gz" cluster_out
    """
}

process CAT {
    cpus 32
    memory '256 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "contigs/${fn}"}
    conda params.conda.catbat
    scratch params.scratch_dir

    input:
    path(asm_fasta)

    output:
    path("cat_out")
    path('cat_out/CAT.contig2classification.txt'), emit: result_table

    shell:
    '''
    mkdir cat_out
    !{params.cat.exe} contigs \
        --compress \
        --tmpdir . --force \
        -n !{task.cpus} \
        -t !{params.cat.tax_db_path} \
        -d !{params.cat.data_db_path} \
        -c !{asm_fasta} -o "cat_out/CAT"
    '''
}

// Workflow 1: Preprocessing
workflow Preprocessing {
    take:
    gfa_file
    hic_r1
    hic_r2

    main:
    ProcessGFA(gfa_file)
    MappabilityMasking(ProcessGFA.out.segments)
    IndexAndMap(
            MappabilityMasking.out.masked,
            hic_r1,
            hic_r2)

    emit:
    segments = ProcessGFA.out.segments
    circular_ids = ProcessGFA.out.circular_ids
    bam_file = IndexAndMap.out.bam_file
}

// Workflow 2: ContactMapAndClustering
workflow ContactMapAndClustering {
    take:
    segments
    bam_file
    excludes

    main:
    MakeContactMap(segments, bam_file)
    ClusterMetagenome(
            segments,
            MakeContactMap.out,
            excludes)

    emit:
    ClusterMetagenome.out
}


// Main workflow to chain everything together
workflow {

    // Run Preprocessing
    Preprocessing(
            params.gfa_file,
            params.hic_r1,
            params.hic_r2)

    // Taxonomic classification of segments
    CAT(Preprocessing.out.segments)

    // Predict virus/plasmid sequences, return confident ids
    MGEAnalysis(
            Preprocessing.out.segments,
            Preprocessing.out.circular_ids,
            CAT.out.result_table)

    // Run ContactMapAndClustering
    ContactMapAndClustering(
            Preprocessing.out.segments,
            Preprocessing.out.bam_file,
            MGEAnalysis.out.confident_mges)

    // Run QualityControl
    QualityControl(
            CAT.out.result_table,
            ContactMapAndClustering.out)
}
