# bin3C
Extract metagenome-assembled genomes (MAGs) from metagenomic data using Hi-C.

## Note on latest codebase: 

At present, users are encouraged checkout the developmental `extraction` or `gothic` branches. These extensions will eventually be merged into the master branch.

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
- Pip or Pipenv for dependency resolution
    - Pipenv and Conda do not place nicely together.
- Git if you wish to clone from our repository: https://github.com/cerebis/bin3C
- GCC is required by some modules within the dependency hierarchy (llvmlite)

### Obtaining bin3C

Installation of bin3C is currently only via Github. Users can either clone the repository with git or obtain an archive of the tarball using github URL syntax.

Clone the repo using git.

```bash
git clone --recursive https://github.com/cerebis/bin3C
```

Pull down the repo using curl.

```bash
mkdir bin3C && curl -L https://api.github.com/repos/cerebis/bin3C/tarball/master | tar -C bin3C --strip-components 1 -xvz 
```
### Dependency resolution

bin3C has a number of Python dependencies, which can be resolved using either Pip or Pipenv.

#### Using Pip

For Pip, we suggest using a virtual environment to resolve dependencies and avoid clashes with system-wide modules or other applications setup by yourself.

This is easiest to set up with the Python module `virtualenv`, which can be obtained from Pip and installed in userspace as follows:

```base
pip install --user virtualenv
```

Once virtualenv is installed, you can install bin3C in its own environment.

```bash
# enter the bin3C folder
cd bin3C
# make a new ve in the bin3C folder
virtualenv .
# using the pip of the ve, installed requirements
bin/pip2 install -r requirements.txt
```

Once complete and assuming you are in the repository folder, bin3C is invoked as follows to make certain we use the local python2 we just set up.

```bash
# run bin3C using the ve python
bin/python2 ./bin3C.py --help
```

#### Using Pipenv

Pipenv can be installed with Pip in userspace as follows:

```bash
pip install --user pipenv
```

Once you have Pipenv, the following will create a virtual environment and install the dependencies.

```bash
# enter the bin3C folder
cd bin3C
# use pipenv to create your ve
pipenv --python 2.7
# now ask pipenv to resolve dependencies
pipenv install
```

To invoke bin3C using the Pipenv environment, things are slightly different. First, we launch a shell which will be preconfigured by Pipenv, then we can confidently invoke bin3C directly.

```bash
# initialise the ve
pipenv shell
# now run bin3C directly
./bin3C.py --help
```

or you can invoke bin3C without going into a deeper shell as follows.

```bash
pipenv run ./bin3C.py --help
```

## Optional - building Infomap

### Requirements
- make
- g++

A statically built Infomap executable has been included in the repository, but this can still cause issues if your runtime environment is particularly old, or you wish to run bin3C on something other than Linux. There is an alternative, which requires a few extra steps.

The informap repository is a submodule of bin3C. If you've checked out the repo with recursion, you should already have the source code to Infomap. In that case, assuming you have `make` and `g++` installed, you can build your own executable very easily from the bin3C root folder.

```bash
# go into the bin3C root folde
cd bin3C
# build infomap for your local environment
make -f Makefile.infomap
```

This will build infomap for source and replace the executable that came with bin3C.

If you do not have the infomap source in your bin3C/external directory, you may have forgotten to clone bin3C with recursion or perhaps you downloaded a tarball. No worries, you can get the submodule without having to pull down the entire bin3C repository again.

```bash
# go into the bin3C root folder
cd bin3C
# get the submodules -- there's only one for now
git submodule update --init --recursive
```

## Typical Workflow

### Initial data preparation

**1. Assemble the metagenome**

We have been using MetaSPAdes (v3.11.1) in our work, but there is no reason that an alternative assembler with a low missassembly rate could be used. One minor caveat at this time however is that while producing the per-cluster report, we cheat a fraction and extract depth of coverage from SPAdes' sequence names rather than calculate this statistic ourselves. Therefore, with other assemblers this information with not be reported at the present time.

Assuming your read-set is in separate Forward/Reverse form, you could create the assembly as so:

```bash
> spades.py --meta -1 shotgun_R1.fastq.gz -2 shotgun_R2.fastq.gz -o spades_out
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
python2 ./bin3C mkmap -e MluCI -v contigs.fasta.gz hic2ctg.bam bin3c_out
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
python2 ./bin3C cluster -v bin3c_out/contact_map.p.gz bin3c_clust
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

