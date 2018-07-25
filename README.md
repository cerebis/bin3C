# bin3C
Extract metagenome-assembled genomes (MAGs) from metagenomic data using Hi-C.

## Introduction
bin3C is a tool which attempts to extract so called metagenome-assembled genones or MAGs from metagenomic shotgun sequencing experiments. The major prerequisite for achieving this is an accompanying Hi-C sequencing data-set.

In short, the shotgun sequencing is assembled and the Hi-C read-set is used to cluster the assembly contigs/scaffolds into genome bins. In our testing on simulated validation datasets, bin3C genome binning solutions have high precision (>0.95), in that extracted genomes possess low contamination in comparison to non-HiC methods and good recall (0.65) for a realistic experimental scenario.

Currently bin3C has only been tested with short-read sequencing data.

## Prerequisites
1. Metagenomic shotgun
  - bin3C was developed using Illuminina short-read data, but other forms resulting in assemblies with high-basecalling accuracy should work as well.
2. Hi-C read-sets from the same sample.
  - For simplicity and good quality data, we reccommend the use of a commercial Hi-C library prep kit.
  - If not using a Hi-C library prep kit, we encourage users to carry-out two digestions with different 4-cutter enzymes. This will the uniformity of coverage over individual contigs in the metagenome.
  - A high quality clean-up of proximity ligation fragments (biotin pulldown) is also very important to maximise signal.

## Typical Workflow
1. Assemble metagenome.
  - We have been using MetaSPAdes (v3.11.1) in our work, but there is no reason another assembler could be used.
2. Map Hi-C reads to assembly contigs or scaffolds.
  - We have been using recent BWA MEM (0.7.15-r1142-dirty) with the following command line: `bwa mem -5SP [contigs] [readset]`
3. Sort the resulting "Hi-C to contigs" BAM file in name order.
4. Create a contact map for analysis
  - Run `python2 ./bin3C mkmap [assembly fasta] [hic bam] [output dir]`
5. Cluster the resulting contact map into genome bins
  - Run `python2 ./bin3C cluster [contact_map] [output dir]`
