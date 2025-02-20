include { SegmentTopology } from './binning_qc'

process VirSorter {
    cpus 24
    memory '100 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/virus/${fn}"}
    conda params.conda.virsorter
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('virsorter_out')
    path('virsorter_out/final-viral-score.tsv'), emit: prediction

    """
    virsorter run --tmpdir ./local_tmp --rm-tmpdir -j ${task.cpus} \
        --include-groups ${params.virsorter.groups} --min-length 2500 \
        --provirus-off -w virsorter_out -i $contigs
    """
}

process Vibrant {
    cpus 24
    memory '100 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/virus/${fn}"}
    conda params.conda.vibrant
    scratch params.scratch_dir

    input:
    path('seqs.fasta')

    output:
    path('vibrant_out')
    path('vibrant_out/VIBRANT_seqs/VIBRANT_phages_seqs/seqs.phages_combined.txt'), emit: prediction

    """
    ${params.vibrant.exe} -t ${task.cpus} -f nucl \
        -l 2500 -i seqs.fasta -folder vibrant_out
    """
}

process Marvel {
    cpus 24
    memory '64 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/virus/${fn}"}
    conda params.conda.marvel
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('marvel_out/all_predictions.tsv'), emit: prediction

    """
    ${params.marvel.exe} ${task.cpus/2} $contigs marvel_out
    """
}

process DeepVirFinder {
    cpus 16
    memory '64 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/virus/${fn}"}
    conda params.conda.dvf
    scratch params.scratch_dir

    input:
    path('seqs.fasta')

    output:
    path('dvf_out')
    path('dvf_out/seqs.fasta_gt2500bp_lt1500000bp_dvfpred.txt'), emit: prediction

    """
    dvf.py -l 2500 -L 1500000 -b 10 -c ${task.cpus} \
        -i seqs.fasta -o dvf_out 
    """
}

process Phamer {
    cpus 16
    memory '128 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/virus/${fn}"}
    conda params.conda.phabox
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('phamer_out')
    path('phamer_out/final_prediction/phamer_prediction.tsv'), emit: prediction

    """
    phabox2 --task phamer --threads ${task.cpus} \
        -d ${params.phabox.db} \
        --contigs $contigs --outpth phamer_out        
    """
}

process PlasClass {
    cpus 24
    memory '128 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/plasmid/${fn}"}
    conda params.conda.plasclass
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('plasclass_out/plasclass.txt'), emit: prediction

    shell:
    '''
    mkdir plasclass_out
    export PATH=!{params.plasclass.path}:$PATH
    classify_fasta.py -f !{contigs} -o plasclass_out/plasclass.txt -p !{task.cpus}
    '''
}

process PlasFlow {
    cpus 24
    memory '128 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/plasmid/${fn}"}
    conda params.conda.plasflow
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('plasflow_out')
    path('plasflow_out/plasflow_0.9.txt'), emit: prediction

    """
    mkdir plasflow_out
    PlasFlow.py --input $contigs --output plasflow_out/plasflow.txt --batch_size 1000
    plasflow_selector.py plasflow_out/plasflow.txt > plasflow_out/plasflow_0.9.txt    
    """
}

process PlasForest {
    cpus 16
    memory '128 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/plasmid/${fn}"}
    conda params.conda.plasforest
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('plasforest_out/plasforest.csv'), emit: prediction

    """
    mkdir plasforest_out
    python3 ${params.plasforest.exe} --threads ${task.cpus} \
        -i $contigs -o plasforest_out/plasforest.csv
    """
}

process PlasmidHunter {
    cpus 16
    memory '128 GB'
    publishDir params.outdir, mode: 'link', saveAs: {fn -> "mge_analysis/plasmid/${fn}"}
    conda params.conda.plasmidhunter
    scratch params.scratch_dir

    input:
    path(contigs)

    output:
    path('plasmidhunter_out/predictions.tsv'), emit: prediction

    """
    mkdir plasmidhunter_out
    plasmidhunter -c ${task.cpus} \
        -i $contigs -o plasmidhunter_out
    """
}

process PrepFasta {
    cpus 1
    memory '16 GB'

    input:
    path(contigs)

    output:
    path('tmp_contigs.fasta')

    script:
    if (contigs.endsWith('.gz')) {
        """
        gzip -d -c $contigs > tmp_contigs.fasta
        """
    } else {
        """
        cp $contigs tmp_contigs.fasta
        """
    }
}

process CollateResults {
    cpus 16
    memory '128 GB'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "mge_analysis/${fn}"}
    conda params.conda.bin3c

    input:
    path('seqs.fasta')
    path('circular.ids')
    path('cat_result.tsv')

    // virus finders
    path(VirSorter)
    path(Vibrant)
    path(Marvel)
    path(DeepVirFinder)
    path(Phamer)

    // plasmid finders
    path(PlasClass)
    path(PlasFlow)
    path(PlasForest)
    path(PlasmidHunter)

    output:
    path("mge_collated")
    path('mge_collated/seqs_mge_2way.ids'), emit: mge_ids

    """
    collate_mge.py --no_partial \
        mge_collated seqs.fasta circular.ids cat_result.tsv \
        $VirSorter $Vibrant $Marvel $DeepVirFinder $Phamer \
        $PlasClass $PlasFlow $PlasForest $PlasmidHunter  
    """
}

workflow MGEAnalysis {

    take:
    segments
    circular_ids
    cat_result

    main:
    // Viruses finders
    VirSorter(segments)
    Vibrant(segments)
    DeepVirFinder(segments)
    Marvel(segments)
    Phamer(segments)

    // Plasmids finders
    PlasClass(segments)
    PlasFlow(segments)
    PlasForest(segments)
    PlasmidHunter(segments)

    CollateResults(
            segments,
            circular_ids,
            cat_result,
            VirSorter.out.prediction,
            Vibrant.out.prediction,
            DeepVirFinder.out.prediction,
            Marvel.out.prediction,
            Phamer.out.prediction,
            PlasFlow.out.prediction,
            PlasClass.out.prediction,
            PlasForest.out.prediction,
            PlasmidHunter.out.prediction)

    emit:
    confident_mges = CollateResults.out.mge_ids
}
