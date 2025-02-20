include { RNA_report } from './rna_report'

process CheckM {
    cpus 24
    memory '128 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "binning_qc/${fn}"}
    conda params.conda.checkm1
    scratch params.scratch_dir

    input:
    path(cluster_dir)

    output:
    path("checkm_out")

    """
    checkm lineage_wf --threads ${task.cpus} --extension ${params.bin_suffix} \
        --tab_table --file "checkm_out/quality.tsv" "${cluster_dir}/fasta" "checkm_out"
    """
}

process CheckM2 {
    cpus 24
    memory '128 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "binning_qc/${fn}"}
    conda params.conda.checkm2
    scratch params.scratch_dir

    input:
    path(cluster_dir)

    output:
    path("checkm2_out")

    """
    checkm2 predict -x ${params.bin_suffix} --threads ${task.cpus} --tmpdir . \
        --input "${cluster_dir}/fasta" --output-directory "checkm2_out"
    """
}

process CocoPye {
    cpus 8
    memory '64 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "binning_qc/${fn}"}
    conda params.conda.cocopye
    scratch params.scratch_dir

    input:
    path(cluster_dir)

    output:
    path("cocopye_out.csv")

    """
    cocopye run -t ${task.cpus} -i "${cluster_dir}/fasta" -o cocopye_out.csv
    """
}

process GTDBtk {
    cpus 32
    memory '400 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "binning_qc/${fn}"}
    conda params.conda.gtdbtk
    scratch params.scratch_dir

    input:
    path(cluster_dir)

    output:
    path("gtdb_out")

    shell:
    '''
    gtdbtk classify_wf -x !{params.bin_suffix} --cpus !{task.cpus} \
        --genome_dir "!{cluster_dir}/fasta" \
        --tmpdir . --out_dir "gtdb_out" \
        --mash_db !{params.gtdbtk.mash_db} 
    '''
}

process CollateResults {
    cpus 1
    memory '8 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "binning_qc/${fn}"}
    conda params.conda.rnareport

    input:
    path(cluster_dir)
    path(rna_report)
    path(checkm_out)
    path(checkm2_out)
    path(cocopye_out)
    path(gtdb_out)

    output:
    path("qc_collated")

    """
    collate_binqc.py \
        qc_collated \
        ${cluster_dir}/cluster_report.csv \
        $rna_report \
        ${checkm_out}/quality.tsv \
        ${checkm2_out}/quality_report.tsv \
        $cocopye_out \
        ${gtdb_out}/gtdbtk.bac120.summary.tsv
    """
}

process SegmentTopology {
    cpus 2
    memory '32 GB'
    conda params.conda.bin3c
    scratch params.scratch_dir
    publishDir params.outdir, mode: 'copy'

    input:
    path(gfa_file)

    output:
    path('circular.csv')
    path('circular.ids'), emit: id_list

    """
    gfa_utils isolates --circular $gfa_file circular.csv
    sed 1d circular.csv | cut -d, -f1 > circular.ids
    """
}

workflow QualityControl {
    take:
    cat_table
    cluster_dir

    main:
    CheckM(cluster_dir)
    CheckM2(cluster_dir)
    CocoPye(cluster_dir)
    GTDBtk(cluster_dir)
    RNA_report(cluster_dir)

    CollateResults(
        cluster_dir,
        RNA_report.out,
        CheckM.out, 
        CheckM2.out, 
        CocoPye.out, 
        GTDBtk.out)

}
