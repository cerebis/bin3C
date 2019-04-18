# bin3C
Extract metagenome-assembled genomes (MAGs) from metagenomic data using Hi-C.

## NEW - docker/singularity based analysis environment

To simplify the setup of a bin3C computational environment for metagenomic analysis, we have recently published docker images which include all the neccessary tools to take raw shotgun and Hi-C reads through to metagenome-assembled genomes [cerebis/bin3c](https://cloud.docker.com/u/cerebis/repository/docker/cerebis/bin3c). These images include the primary tools we ourselves have chosen to take raw shotgun and Hi-C reads through to genome bins. 

In the near future, this image will also include a complete workflow.

Two images are available, which differ in the age of the supported kernel.

**Current distributions: cerebis/bin3c:latest**

This image should support the majority of users and is based upon Fedora release 29.

**Older distributions: cerebis/bin3c:centos6**

This image supports unforutnately souls, whose computational environments may be getting out of date. Built upon Centos Release 6 and the 2.6 Linux kernel.

In either case, users should be able to deploy these images either using Docker or Singularity. With Singularity, we reccomend that users specify a contained environment (`-c`) and a home directory (`-H`).

With Docker
```
docker pull cerebis/bin3c:latest
docker run cerebis/bin3c:latest bin3C -h
```

With Singularity
```
singularity build bin3c.img docker://cerebis/bin3c:latest
singularity run -c -H $PWD:/home/user bin3c.img bin3C -h
```

## Introduction
bin3C is a tool which attempts to extract so called metagenome-assembled genones or MAGs from metagenomic shotgun sequencing experiments. The major prerequisite for achieving this is an accompanying Hi-C sequencing data-set.

In short, the shotgun sequencing is assembled and the Hi-C read-set is used to cluster the assembly contigs/scaffolds into genome bins. In our testing on simulated validation datasets, bin3C genome binning solutions have high precision (>0.95), in that extracted genomes possess low contamination in comparison to non-HiC methods and good recall (0.65) for a realistic experimental scenario.

Currently bin3C has only been tested with short-read sequencing data.

## Data requirements and recommendations
* A metagenomic shotgun sequencing read-set.
    - bin3C was developed using Illuminina short-read data, but other forms resulting in assemblies with high-basecalling accuracy should work as well.
* A Hi-C read-set derived from the same sample.
    - For simplicity and good quality data, we reccommend the use of a commercial Hi-C library prep kit.
    - If not using a Hi-C library prep kit, we encourage users to carry-out two digestions with different 4-cutter enzymes. This will the uniformity of coverage over individual contigs in the metagenome.
    - A high quality clean-up of proximity ligation fragments (biotin pulldown) is also very important to maximise signal.

## External software requirements
- Metagenomic assembly software
    - We have been using SPAdes http://cab.spbu.ru/software/spades/
- Short-read aligner
    - We have been using BWA https://github.com/lh3/bwa
- Samtools or other SAM/BAM handling package
    - https://github.com/samtools/samtools

## Runtime and setup
- Python 2.7
- Pip (>=v19) for installation dependency resolution
- GCC is required by some dependencies (eg. llvmlite)

### Obtaining bin3C

The bin3C development branch can be installed easily using pip, but we strongly encourage users to make use of a virtual python environment to avoid dependency version conflicts. 

Detailed below are the steps involved in setting up a working application. After pip completes successfully, the bin3C executable entry-point can be found in your python environment's bin directory.

***Note**: bin3C has a number of Python dependencies and a few external binary executables, these will all be resolved automatically when using pip, however NumPy and Cython must be installed first.*

Eg. To install the development branch.
```bash
# 1. create a virtual environment
mkdir bin3C_env
virtualenv -p python2.7 bin3C_env
cd bin3C_env

# 2. make sure pip is up to date
bin/pip install -U pip

# 3. install numpy and cython
bin/pip install "numpy<1.15" cython

# 4. install bin3C
bin/pip install git+https://github.com/cerebis/bin3C@pgtk

# 5. get bin3C's help
bin/bin3C -h
```
## Typical Workflow

### Initial data preparation

**1. Assemble the metagenome**

We have been using MetaSPAdes (v3.11.1) in our work, but there is no reason that an alternative assembler with a low missassembly rate could be used. One minor caveat at this time however is that while producing the per-cluster report, we cheat a fraction and extract depth of coverage from SPAdes' sequence names rather than calculate this statistic ourselves. Therefore, with other assemblers this information with not be reported at the present time.

Assuming your read-set is in separate Forward/Reverse form, you could create the assembly as so:

```bash
spades.py --meta -1 shotgun_R1.fastq.gz -2 shotgun_R2.fastq.gz -o spades_out
```

_**Optional step. Split references into fragments**_

Prior to mapping Hi-C reads to assembly contigs, the contigs can be broken into smaller pieces using the tool `split_ref.py`. 

Our experiments have shown that splitting may improve both Precision and Recall (i.e. improve the reconstruction of genomes), compensating, we surmise, for assembly errors. It, however, should be noted that the evidence was not unamimiously favourable. Further, `split_ref.py` currently takes only a simplistic approach and makes no attempt to identify the points at which these suspected errors occur (or other features which contradict internal assumptions).

If you choose to perform this step, we suggest a target fragment size no smaller than 5kb, with the default being 10kb. For a given target size, any sufficiently long contig will be split into smaller pieces. As majority of contigs will not evenly divide by the chosen target size, the remainder (or shortfall) is distributed across the fragments. Using this approach, fragments are uniform in size within a contig and close to the target size but between contigs can very slightly.

The resulting split contigs will retain their original identifiers, but each will be appended with the coordinates of the split to retain uniqueness within the multi-fasta file.

e.g. Splitting a 37kb contig `NODE_887_length_37685_cov_9.578421` with a target size of 10kb, the fragments would be as follows:

  - `NODE_887_length_37685_cov_9.578421.0_9421`
  - `NODE_887_length_37685_cov_9.578421.9421_18842`
  - `NODE_887_length_37685_cov_9.578421.28263_37685`

**2. Map the Hi-C reads to assembly contigs or scaffolds**

We have been using recent versions of BWA MEM, which have the option `-5` (such as v0.7.15-r1142-dirty). Assuming your Hi-C reads are interleaved, with BWA MEM we recommend users ignore mate-pairing steps (`-SP`) and importantly require 5-prime alignments are primary (`-5`).

```bash
bwa mem -5SP contigs.fasta.gz hic_paired.fastq.gz | \
    samtools view -bS - > hic2ctg_unsorted.bam
```

**3. Sort the resulting "Hi-C to contigs" BAM file in name order**

Downstream, bin3C expects read-pairs to come sequentially which requires the BAM file be sorted in name order. In future, we may internalise this step.

```bash
samtools sort -o hic2ctg.bam -n hic2ctg_unsorted.bam
```

Steps **2** and **3** can be combined with some pipes on a single line, where we can also filter out alignments which will not contribute to the map and save processing time downstream. These being read-pairs where either is flagged secondary, supplementary and unmapped. Depending on your environment, you may wish to add concurrency with BWA (`-t`) and Samtools command (`-@`).

```bash
bwa mem -5SP contigs.fasta.gz hic_paired.fastq.gz | \
    samtools view -F 0x904 -bS - | \
    samtools sort -o hic2ctg.bam -
```

### bin3C analysis

We have split bin3C into two primary stages. 
1) mkmap: Create a contact map
2) cluster: Clustering the contact map (genome binning)

**4. Create a contact map for analysis**

The contact map is the primary data object on which the genome binning is performed. The map is derived from the assembly sequences chosen (scaffolds or contigs) and the corresponding Hi-C alignments BAM file.

Assuming Pip-style invocation, a map for a library constructed using the MluCI restriction enzyme and verbose output could be made as follows:

```bash
# inside bin3C's virtual env
bin/bin3C mkmap -e MluCI -v contigs.fasta.gz hic2ctg.bam bin3c_out
```

While running, the user will be presented with progress information. 

All files from this stage will be stored in the output directory, including a detailed log. In particular, the compressed contact map is named: `contact_map.p.gz`.

#### Output directory contents after stage 1.

|#| Filename         | Description                      |
|-|------------------|----------------------------------|
|1| bin3C.log        | Appending log file               |
|2| contact_map.p.gz | contact map results from stage 1 |



**5. Cluster the resulting contact map into genome bins**

After a map has been created, the second stage of analysis is clustering the map from stage one. 

Using defaults aside from verbose output, and storing the results from the clustering stage in a directory called `bin3c_clust`, the clustering is performed as follows:

```bash
# inside bin3C's virtual env
bin/bin3C cluster -v bin3c_out/contact_map.p.gz bin3c_clust
```

#### Output directory contents after stage 2

Assuming the same directory was used in both stages, all analysis files will be located together.

|# | Filename           | Description                                       |
|--|--------------------|---------------------------------------------------|
| 1| bin3C.log          | Appending log file                                |
| 2| clustering.mcl     | Clustering solution in MCL format                 |
| 3| clustering.p.gz    | Clustering solution as a picked python dictionary |
| 4| cluster_plot.png   | Heatmap of the contact map after clustering       |
| 5| cluster_report.csv | A per-cluster report of various statistics        | 
| 6| cm_graph.edges     | The graph used in clustering in edge list format  |
| 7| cm_graph.tree      | Infomap clustering output                         |
| 8| fasta              | Per-cluster multi-fasta sequences                 |
| 9| infomap.log        | Infomap runtime log                               |

## Inspecting results

#### Fasta for each genome bin 
After stage 2 has been successfully completed, binned assembly contigs can be found in the `fasta/` subdirectory. We suggest analysing these with CheckM, BUSCO for inferring completeness and contamination. GTDBtk can be used for a more thorough taxonomic analysis.

#### Clustering result
The most portable format of the clustering result is `clustering.mcl` which follows the format established by MCL (surprisingly). In this format, each line pertains to a cluster (counting from 1), while the contents of a line is a space-separated list of member sequences.

E.g. For a MCL file of three lines with the following contents:

```
L1: seq1 seq4 seq5
L2: seq2
L3: seq2
```

There are 3 clusters, with memberships as follows:

- Cluster1: seq1, seq4, seq5
- Cluster2: seq3
- Cluster3: seq2

#### Qualitative inspection of result

Users are recommended to inspect the heatmap, to qualitatively appraise the result (see below), where heatmap pixel intensity is the logarithm of the Hi-C interaction strength. The heatmap is downsampled to a maximum dimension 4000x4000 due to the need to create a density matrix representation. Using a higher resolution is possible if memory is available.

When inspecting a heatmap, signal should be concentrated down the diagonal in crisp blocks of varying size. Relatively intense unattached off-diagonal blocks *might* indicate bin splitting errors. For data with good signal to noise, the overall appearance of the off-diagonal field should be dark or substantially much less intense than the diagonal. 

Good clustering and signal to noise. ![good map](https://drive.google.com/uc?id=1MZNRmU4PwTwkdI4WXU_qG9d9kIMGmoQl )

Potential clustering errors and poor signal to noise.
*( add bad heatmap )*
