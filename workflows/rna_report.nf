process Aragorn {
    cpus 1
    memory '16 GB'
    conda params.conda.rnareport

    input:
    path(fasta_files)

    output:
    path('aragorn_results.csv')

    shell:
    '''
    for fn in !{fasta_files.join(" ")};
    do
        echo $(basename $fn .fna),$(aragorn -w -l -t $fn | aragorn_parser.py --total-only)
    done > aragorn_results.csv
    '''
}

process Barrnap {
    cpus 4
    memory '32 GB'
    conda params.conda.rnareport

    input:
    path(fasta_files)

    output:
    path('barrnap_results.csv')

    shell:
    '''
    for fn in !{fasta_files.join(" ")};
    do
        barrnap --threads !{task.cpus} $fn 2>/dev/null | barrnap_parser.py --name $(basename $fn .fna)
    done | awk -F, \'!unique[$8]++\' > barrnap_results.csv
    '''
}

workflow CountRNA {
    take:
    cluster_dir

    main:

    x = cluster_dir.map{ dir ->
        def bin_dir = dir / 'fasta/'
        file(file(bin_dir).resolve("*${params.bin_suffix}"))
    }.flatten()

    Aragorn(x | buffer(size: params.rna_report.batch_size, remainder: true))
        .collectFile(name: 'aragorn_report.csv', sort: true)
        .set{ tRNA }

    Barrnap(x | buffer(size: params.rna_report.batch_size, remainder: true))
        .collectFile(name: 'barrnap_report.csv', sort: true, keepHeader: true, skip: 1)
        .set{ rRNA }

    emit:
    tRNA = tRNA
    rRNA = rRNA
}

process CombineReports {
    cpus 1
    memory '8 GB'
    conda params.conda.rnareport
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "${publish_subdir}/rna/${fn}"}

    input:
    path(tRNA_report)
    path(rRNA_report)
    val(publish_subdir)

    output:
    path('combined.csv')

    """
    #!/usr/bin/env python
    import pandas as pd
    df1 = pd.read_csv("$tRNA_report", names=["bin","tRNA_count"]).set_index("bin")
    df2 = pd.read_csv("$rRNA_report").rename(columns={"name": "bin"}).set_index("bin")
    df1.join(df2, how="inner").reset_index().to_csv("combined.csv", index=False)
    """
}

workflow RNA_report {
    take:
    cluster_dir
    publish_subdir

    main:
    CountRNA(cluster_dir)
    CombineReports(CountRNA.out.tRNA,
                   CountRNA.out.rRNA,
                   publish_subdir)

    emit:
    CombineReports.out
}
