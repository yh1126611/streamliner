# streamliner

Streamline from assembly, annotation and raw reads to modification probability (MP) score information on transcription start & termination sites (TSS/TES). Output may be used as data for visualization in R.

## Installation
### Conda
streamliner can be installed via conda
```
conda install yh1126::streamliner
```
### Conda recipe
Alternatively, streamliner can also be built locally from conda recipe following [conda package building manual](https://docs.conda.io/projects/conda-build/en/latest/user-guide/tutorials/build-pkgs.html)

* If not already installed, install conda-build:

        conda install conda-build

* Download [streamliner.sh](conda/streamliner.sh), [build.sh](conda/build.sh) and [meta.yaml](conda/meta.yaml) and move them all under a same directory.
* After moving to the same directory, build package:

        conda build .

* Install locally:

        conda install --use-local streamliner

## Input preparation (BURST)
![Species information table format](img/burst.png?raw=true "BURST")

Input for streamliner is __BURST__, a text file representing table of all data required streamliner, created as following:

1. Create a table comprising RefSeq accession ID from [NCBI Database](https://www.ncbi.nlm.nih.gov/) and 42basepairs download links for their raw PacBio HiFi reads (BAM) from [GenomeArk](https://www.genomeark.org/) for every species you wish to analyze, parsed into a two-column table-format file (e.g. by MS Excel). Only BAM files with kinetics tags (`fi`, `fp`, `ri`, `rp`) generated from PacificBioscience's [5-base sequencing](https://www.pacb.com/wp-content/uploads/application-brief-measuring-dna-methylation-with-5-base-hifi-sequencing.pdf) feature may be included in input. Every row of the table corresponds to one species and two columns correspond to accession ID and download links, respectively. Multiple BAM files may be associated with one species, in which case, download links should be separated by newlines (\n) inside a single cell. ![Diagram of input preparation](img/input_prep_1.png?raw=true "Diagram of input preparation")  
  Schematic format of a species information table: ![Species information table format](img/species_information_table_format.png?raw=true "Species information table format")  

2. Export table into a tab-delimited text format. Columns are separated by tabs (\t) and rows are separated by newlines (\n). A multi-line cell with multiple links which are also separated by newlines (\n) are enclosed by double quotes (").  
   Example arrangement of a species information table in text format:
  
        <Accession1>	"<Link1_1>
        <Link1_2>
        <Link1_3>
        <Link1_n>"
        <Accession2>	"<Link2_1>
        <Link2_2>
        <Link2_3>
        <Link2_n>"  
   Schematic diagram of BURST in text format: ![Schematic diagram of a species table converted to text format](img/species_information_table_txt_format.png?raw=true "Schematic diagram of a species table converted to text format")  

## Sample input
[This BURST](input/BURST_sample.txt) is a functional input to streamliner comprising GenomeArk download links to all raw PacBio HiFi reads for every VGP-sequenced species for which BAM files with kinetics tags exist (78).

## Usage

```
streamliner [OPTIONS] <species_information>
```

## Options

* `-t`: Set feature to use as locational indicator for TSS/TES. Possible choices are: gene [Default], transcript and mRNA.

* `-m`: Set maximum memory allocated for sorting process. Suffixes K/M/G are accepted. Maximum is 99G. [Default 99G]

* `-i`: Set furthest distance from TSS/TES to be analyzed for MP (i.e. interval ÷ 2) (bp) [Default 10,000]

* `-w`: Set size of window inside interval to calculate GC and BpB contents (bp) [Default 100]  
     Note: `-w` (window size) must be smaller than `-i` (half of interval size).

* `-p`: Run streamliner partially by entering two unique values to specify stages to start and end the pipeline. Default is 0 10 (Complete pipeline). The unique values for each stage are as follows:

     ０- Downloading input (i.e. complete pipeline) [Default start]  
     １- Merging (Engine: [samtools merge](https://www.htslib.org/doc/samtools-merge.html))  
     ２- Extracting HiFi (Engine: [pbtk](https://github.com/PacificBiosciences/pbtk))  
     ３- Modification-calling (Engine: [pbjasmine](https://github.com/PacificBiosciences/jasmine))  
     ４- Alignment (Engine: [pbmm2](https://github.com/PacificBiosciences/pbmm2))  
     ５- Sorting (Engine: [samtools sort](https://www.htslib.org/doc/samtools-sort.html))  
     ６- Indexing (Engine: [samtools index](https://www.htslib.org/doc/samtools-index.html))  
     ７- Modification probability (MP) computation (Engine: [pb-cpg-tools](https://github.com/PacificBiosciences/pb-CpG-tools))  
     ８- Extraction of MP at transcription start sites (TSS) & transcription termination sites (TES) (Engine: [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html))  
     ９- GC content calculation (Engine: [streamGC](https://github.com/yh1126611/streamGC), [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html))  
     10 - Calculation of default BpB (CpC, CpG, GpC, GpG) contents (Engine: [streamBpB](https://github.com/yh1126611/streamBpB), [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html)) [Default end] 

     Note: For a partial run, intermediate files from all preceeding stages for every constituent in the input species information must exist in respective folders, appropriately-named.

## Notes
* 42basepairs should be used as opposed to s3 for download links in species information input.
* Setting the interval size by `-i` will only affect GC and BpB calculations and not MP calculation.
* Sorting part of the pipeline comprises dividing the BAM file into smaller files of size permitted by maximum memory (maxMem) allocated to process. A smaller maxMem will require dividing the BAM file into more number of files and too many files may breach the system’s limit for number of open files. Therefore, maxMem to process must be sufficiently large to allow the BAM file to be split into a number less than the maximum number of open files (preferably over 20G) while being smaller than available memory. Check your system-permitted capacity for number of open files by `ulimit -n` and available memory by `free -mh` and allocate maxMem to streamliner accordingly using `-m` option. The default is set to maximum value `99G` which can be larger than available memory of some systems depending on number of simultaneous users and processes in which case should be reduced.
