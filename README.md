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

## Input preparation

1. Create a __species information table__ comprising every species you wish to analyze by obtaining their RefSeq accession ID from [NCBI Database](https://www.ncbi.nlm.nih.gov/) and 42basepairs download links for their raw PacBio HiFi reads (BAM) from [GenomeArk](https://www.genomeark.org/) and parsing them into a two-column table-format file (e.g. excel sheet). Only BAM files generated by [5-base sequencing](https://www.pacb.com/wp-content/uploads/application-brief-measuring-dna-methylation-with-5-base-hifi-sequencing.pdf) thus including kinetics tags (`fi`, `fp`, `ri`, `rp`) may be included. Every row of the table corresponds to one species and two columns correspond to accession ID and download links, respectively. Multiple BAM files may be associated with one species, in which case, download links should be separated by newlines (\n) inside a single cell. ![Diagram of input preparation](img/input_prep_1.png?raw=true "Diagram of input preparation")  
  Schematic format of a species information table: ![Species information table format](img/species_information_table_format.png?raw=true "Species information table format")  

2. Export the species information table into a tab-delimited text format. Columns are separated by tabs (\t) and rows are separated by newlines (\n). A multi-line cell with multiple links which are also separated by newlines (\n) are enclosed by double quotes (").  
   Example arrangement of a species information table in text format:
  
        <Accession1>  "<Link1_1>
        <Link1_2>
        <Link1_3>
        <Link1_n>"
        <Accession2> "<Link2_1>
        <Link2_2>
        <Link2_3>
        <Link2_n>"  
   Schematic diagram of a species information table in text format: ![Schematic diagram of a species table converted to text format](img/species_information_table_txt_format.png?raw=true "Schematic diagram of a species table converted to text format")  

## Sample input
[This text-converted species information table](input/species_information.txt) is a functional input to streamliner comprising GenomeArk download links to all raw PacBio HiFi reads for every VGP-sequenced species for which BAM files with kinetics tags exist (78).

## Usage

```
streamliner [OPTIONS] <species_information>
```

## Options

* `-t`: Set feature to use as locational indicator for TSS/TES [Default "gene"]

* `-m`: Set maximum memory allocated for sorting process. Suffixes K/M/G are accepted [Default 99G]

* `-i`: Set furthest distance from TSS/TES to be analyzed for MP (i.e. interval ÷ 2) (bp) [Default 10,000]

* `-w`: Set size of window inside interval to calculate GC and BpB contents (bp) [Default 100]  
     Note: Window size (`-w`) must be smaller than half of interval size (`-i`).

* `-p`: Run streamliner partially by entering a unique value to specify a stage to start the pipeline. The unique values for each stage are as follows:

     ０- Downloading input (i.e. complete pipeline) [Default]  
     １- Merging (Engine: [samtools merge](https://www.htslib.org/doc/samtools-merge.html))  
     ２- Extracting HiFi (Engine: [pbtk](https://github.com/PacificBiosciences/pbtk))  
     ３- Modification-calling (Engine: [pbjasmine](https://github.com/PacificBiosciences/jasmine))  
     ４- Alignment (Engine: [pbmm2](https://github.com/PacificBiosciences/pbmm2))  
     ５- Sorting (Engine: [samtools sort](https://www.htslib.org/doc/samtools-sort.html))  
     ６- Indexing (Engine: [samtools index](https://www.htslib.org/doc/samtools-index.html))  
     ７- Modification probability (MP) computation (Engine: [pb-cpg-tools](https://github.com/PacificBiosciences/pb-CpG-tools))  
     ８- Extraction of MP at transcription start sites (TSS) & transcription termination sites (TES)  
     ９- GC content calculation (Engine: [streamGC](https://github.com/yh1126611/streamGC))  
     10 - Calculation of default BpB (CpC, CpG, GpC, GpG) contents (Engine: [streamBpB](https://github.com/yh1126611/streamBpB))  

     Note: For a partial run, intermediate files of all preceeding stages for every constituent in the input species information must exist in respective folders.

## Notes
* Sorting part of the pipeline comprises dividing the BAM file into smaller files of size permitted by maximum memory (maxMem) allocated to process. A smaller maxMem will require dividing the BAM file into more number of files. Too many files may breach the system’s limit for number of open files. Therefore, maxMem to process must be smaller than available memory but large enough to allow the BAM file to be split into a number less than the maximum number of open files. Check your system-permitted capacity for number of open files by `ulimit -n` and available memory by `free -mh` and allocate maxMem to streamliner accordingly using `-m` option. The default is set to maximum value `99G` which can be larger than available memory of some systems depending on number of simultaneous users and processes.
