# streamliner

Streamline from assembly, annotation and raw reads to modification probability score information on transcription start/end sites (TSS/TES). Resulting output may be used as data for visualization in R.

## Installation
### Conda
streamliner can be installed via conda
```
conda install yh1126::streamliner
```
### Conda recipe
Alternatively, streamliner can also be installed built locally using conda recipe following [conda package building manual](https://docs.conda.io/projects/conda-build/en/latest/user-guide/tutorials/build-pkgs.html)

* If not already installed, install conda-build:

        conda install conda-build

* Download [streamliner.sh](streamliner.sh), [build.sh](build.sh) and [meta.yaml](meta.yaml) and move them all under a same directory.
* After moving to the same directory, build package:

        conda build .

* Install locally:

        conda install --use-local streamliner

## Input preparation

1. For every species you wish to analyze, create a species information table by obtaining RefSeq accession ID from NCBI and 42basepairs download links for raw PacBio HiFi reads (BAM) from GenomeArk (https://www.genomeark.org/) and parsing them into a two-column table-format file (e.g. excel sheet). Every row pertains to one species and two columns each pertain to accession ID and download links, respectively. Multiple BAM files (thus, multiple download links) may be associated with one species, in which case, they should be separated by newlines (\n) inside a single cell. ![Diagram of input preparation](input_prep_1.png?raw=true "Diagram of input preparation")  
  Schematic format of a species information table: ![Species information table format](species_information_table_format.png?raw=true "Species information table format")  
  Example of a species information table: ![Species information table example](species_information_table_example.png?raw=true "Example of a species information table")

2. Export the table consisting of information for your species to be analyzed into a tab-delimited text format. Columns are separated by tabs (\t) and rows are separated by newlines (\n). A multi-line cell with multiple links which are also separated by newlines (\n) are enclosed by double quotes (").  
  
  
        <Accession1>  "<Link1_1>
        <Link1_2>
        <Link1_3>
        <Link1_n>"
        <Accession2> "<Link2_1>
        <Link2_2>
        <Link2_3>
        <Link2_n>"  
  
  Schematic diagram of a species information table in text format: ![Schematic diagram of a species table converted to text format](species_information_table_txt_format.png?raw=true "Schematic diagram of a species table converted to text format")  
Example of a species information table converted to a text format: ![Example of a species information table in text format](species_information_table_txt_example.png?raw=true "Example of a species information table in text format")

## Notes

* Sorting part of the pipeline comprises dividing the BAM file into smaller files of size permitted by maximum memory (maxMem) allocated to process. A smaller maxMem will require dividing the BAM file into more number of files. Too many files may breach the system’s limit for number of open files. Therefore, maxMem to process must be smaller than available memory but large enough to allow the BAM file to be split into a number less than the maximum number of open files. Check your system-permitted capacity for number of open files by `ulimit -n` and available memory by `free -mh` and allocate maxMem to streamliner accordingly using `-m` option. The default is set to maximum value `99G` which can be larger than available memory of some systems depending on number of simultaneous users and processes.
