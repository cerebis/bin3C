process AnalyseContacts {
    cpus 2
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(map_dir)
    path(cluster_dir)
    path(coverage)
    path(genmap_text)

    output:
    path("sig3c_out")

    """
    sig3C analyse --clobber --seed ${params.bin3c.seed} \
        "${map_dir}/contact_map.p.gz" \
        "${cluster_dir}/clustering.p.gz" \
        $coverage \
        $genmap_text \
        sig3c_out
    """
}

process ComputeEmbeddings {
    queue 'gpuq'
    cpus 8
    memory '64 GB'
    accelerator 1
    conda params.conda.seq_embed
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(segments)

    output:
    path('embeddings.p.gz')

    """
    embed.py \
        --seed ${params.bin3c.seed} \
        --device-id ${params.seq_embed.device_id} \
        --chunk-size ${params.seq_embed.chunk_size} \
        --batch-size ${params.seq_embed.batch_size} \
        $segments embeddings.p.gz
    """
}

process LabelTrainingData {
    cpus 2
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(embeddings)
    path(cluster_dir)
    path(fasta_index)
    path(qc_dir)
    path(sig3c_dir)

    output:
    path(sig3c_dir)

    """
    sig3C labeller --clobber --seed ${params.bin3c.seed} \
        $embeddings \
        "${cluster_dir}/clustering.p.gz" \
        ${fasta_index} \
        "${qc_dir}/qc_summary.csv" \
        ${sig3c_dir}
    """
}

process ClassifyContacts {
    cpus 2
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(sig3c_dir)

    output:
    path(sig3c_dir)

    """
    sig3C classify --verbose --clobber --seed ${params.bin3c.seed} \
        --learning-rate ${params.sig3c.learning_rate} \
        --n-epochs ${params.sig3c.n_epochs} \
        --n-nodes ${params.sig3c.n_nodes} \
        --n-layers ${params.sig3c.n_layers} \
        --early-stopping \
        ${sig3c_dir}
    """
}

workflow PredictIntracellularContacts {
    take:
    map_dir
    cluster_dir
    coverage
    genmap_text
    segments
    fasta_index
    qc_dir

    main:
    AnalyseContacts(map_dir,
                    cluster_dir,
                    coverage,
                    genmap_text)

    ComputeEmbeddings(segments)

    LabelTrainingData(ComputeEmbeddings.out,
                      cluster_dir,
                      fasta_index,
                      qc_dir,
                      AnalyseContacts.out)

    ClassifyContacts(LabelTrainingData.out)

    emit:
    classify_out = ClassifyContacts.out
}
