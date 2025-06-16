#!/bin/bash

if [ "$#" -eq 0 ]; then
    "$0" --help
    exit 0
fi

if [[ "$1" == "--version" ]]; then
    echo "Version 1.2.0"
    exit 0
fi

if [[ "$1" == "--help" ]]; then
    echo ""
    echo "Program:        streamliner"
    echo "Version:        1.2.0"
    echo "Code:           https://github.com/yh1126611/streamliner"
    echo ""
    echo "Description:"
    echo "Streamline from assembly, annotation and raw reads to modification"
    echo "probability score information on transcription start/end sites (TSS/TES)."
    echo "Resulting output may be used as data for visualization in R."
    echo ""
    echo "Usage:"
    echo "streamliner [OPTIONS] <BURST.txt|tsv>"
    echo ""
    echo "<BURST.txt|tsv>         A two-column tab-delimited file with following"
    echo "                        columns:"
    echo "                        1. NCBI RefSeq accession ID"
    echo "                        2. Links to PacBio HiFi sequencing reads on"
    echo "                           GenomeArk (https://www.genomeark.org/)"
    echo ""                                     
    echo "                        The reads are BAM files and must contain"
    echo "                        kinetics tags (fi, fp, ri, rp)."
    echo "                        In case there are multiple read files for one"
    echo "                        species, a multi-line cell is allowed where"
    echo "                        links are separated by newline characters(\n)"
    echo "                        and all links for one species are enclosed by"
    echo "                        double quotes(\")."
    echo ""
    echo "Optional:"
    echo "-t, --type      STRING  Feature to use as locational indicator for"
    echo "                        TSS/TES [Default gene]"
    echo "-m              INT     Maximum memory allocated for sorting process."
    echo "                        Suffix K/M/G allowed. Cannot exceed 99G."
    echo "                        [Default 99G]"
    echo "-i, --interval  INT     Furthest distance from coordinate to be"
    echo "                        analyzed (bp) [Default 10,000]"
    echo "-w, --window    INT     Size of each window inside interval (bp)."
    echo "                        [Default 100]"
    echo "-p, --partial   INT     Two values between 0 and 10 indicating which"
    echo "                        stages to start and end streamliner,"
    echo "                        respectively if pipeline is to be run"
    echo "                        partially [Default 0, 10 (Complete pipeline)]"
    echo ""
    exit 0
fi

# Initialize defaults
maxmem="99G"
type_string="gene"
window_size=100
interval_size=10000
partial_start=0
partial_end=10

# Manual option parsing to handle both short and long options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m)
            maxmem="$2"
            shift 2
            ;;
        -t|--type)
            type_string="$2"
            shift 2
            ;;
        -w|--window)
            window_size="$2"
            shift 2
            ;;
        -i|--interval)
            interval_size="$2"
            shift 2
            ;;
	-p|--partial)
            partial_start="$2"
            partial_end="$3"
            shift 3
            ;;
        -*)
            echo "Invalid option: $1" >&2
            exit 1
            ;;
        *)
            # Break when non-option arguments are found
            break
            ;;
    esac
done

read_links=$1

if [ -z "$read_links" ]; then
    echo "Error: Missing input file" >&2
    exit 1
fi

if [ ! -f "$read_links" ]; then
    echo "Error: File $read_links not found" >&2
    exit 1
fi

if [ $# -ne 1 ]; then
    echo "Usage: streamliner [OPTIONS] <species_information>"
    exit 1
fi

for var_name in interval_size window_size; do
  value="${!var_name}"
  if [[ -n "$value" && ! "$value" =~ ^-?[0-9]+$ ]]; then
    echo "Error: Option -${var_name:0:1} must be an integer." >&2
    exit 1
  fi
done

if [[ -n "$maxmem" && ! "$maxmem" =~ ^-?[0-9]+[KMG]?$ ]]; then
	echo "Error: Option -m must be an integer or an integer followed by K, M, or G (e.g. 100M)." >&2
	exit 1
fi

# Check if interval is bigger than window
if [ "$interval_size" -le "$window_size" ]; then
    echo "Error: Interval size ($interval_size) must be larger than window size ($window_size)" >&2
	exit 1
fi

if [[ -n "$type_string" && ! "$type_string" =~ ^(gene|transcript|mRNA)$ ]]; then
	echo "Error: Option -t must be one of 'gene', 'transcript', or 'mRNA'." >&2
	exit 1
fi

if ! [[ "$partial_start" =~ ^[0-9]+$ && "$partial_end" =~ ^[0-9]+$ ]]; then
	echo "Error: -p|--partial requires two integer arguments between 0 and 10." >&2
	exit 1
fi

if [ "$partial_start" -le 0 ] && [ "$partial_end" -ge 0 ]; then
	# Extract accession numbers
	awk -v FS=" " -v OFS="\t" '$1~/^GCF/{print $1}' "$read_links" | while read accession; do
		echo $accession >> accession.txt
	done

	# Determine assembly name registered in NCBI
	awk -v FS=" " -v OFS="\t" '$1~/^GCF/{print $1}' "$read_links" | while read accession; do
		datasets summary genome accession $accession | grep assembly_name | awk '{gsub(/"/, "", $0); gsub(/ /, "_", $0); print $0}' | awk -v FS=",|:" '{for (i=1; i<=NF; i++) if ($i~/^assembly_name$/) {print $(i+1)}}' | awk -v FS="\t" -v OFS="\t" '{print $0}' >> assembly_name.txt
	done

	# Write head of FTP download link
	awk -v FS=" " -v OFS="\t" '$1~/^GCF/{print $1}' "$read_links" | while read accession; do
		echo $accession | awk -v FS="_|\." '{print $2}' | awk -v FS="" -v OFS="" '{print "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/", $1, $2, $3, "/", $4, $5, $6, "/", $7, $8, $9, "/"}' >> ftp.txt
	done

	# Write full FTP download link (less assembly identifier (_genomic.fna.gz) or annotation identifier (_genomic.gtf.gz))
	paste -d "_" <(paste -d "/" <(paste -d "_" <(paste -d "" ftp.txt accession.txt) assembly_name.txt) accession.txt) assembly_name.txt > ftp_url.txt

	# Create folders by assembly name and download assembly and annotation to respective folder
	paste <(awk -v FS="/|\t" -v OFS="\t" '$1~/^GCF/{print $10}' "$read_links") ftp_url.txt | while read assembly_name ftp_url; do
		mkdir -p "$assembly_name"
		cd "$assembly_name"
		wget "${ftp_url}_genomic.fna.gz"
		wget "${ftp_url}_genomic.gtf.gz"
		gzip -df *.gz
		for file in *.fna; do
			mv "$file" "${file%.fna}.fasta"
		done
		cd ..
	done

	# Download reads
	awk -F'\t' '{gsub(/"/, "", $NF); gsub(/\r$|\n$/, "", $NF); print $NF}' "$read_links" | awk -v FS="/|\t" -v OFS="\t" '{print $9, $0}' | while read -r col1 col2; do
		cd "$col1"
		wget "$col2"
		cd ..
	done
fi

# Pipeline (modification score calling -> TSS/TES profile identification -> GC/BpB content calculation)
awk -v FS="\t|/" -v OFS="\t" '$1~/^GCF/{print $10}' "$read_links" | while read -r col1; do
	cd "$col1"
	if [ "$partial_start" -le 1 ] && [ "$partial_end" -ge 1 ]; then
		samtools merge "$col1".bam *.bam
	fi
	if [ "$partial_start" -le 2 ] && [ "$partial_end" -ge 2 ]; then
		[ -f "$col1".bam ] || { echo "Error: Merged BAM file for ${col1} not found"; exit 1; }
		extracthifi "$col1".bam "$col1"_hifi.bam
	fi
	if [ "$partial_start" -le 3 ] && [ "$partial_end" -ge 3 ]; then
		[ -f "$col1"_hifi.bam ] || { echo "Error: HiFi-extracted BAM file for ${col1} not found"; exit 1; }
		jasmine "$col1"_hifi.bam "$col1"_hifi_methCalled.bam
	fi
	if [ "$partial_start" -le 4 ] && [ "$partial_end" -ge 4 ]; then
		[ -f "$col1"_hifi_methCalled.bam ] || { echo "Error: Modification-called BAM file for ${col1} not found"; exit 1; }
		pbmm2 align *.fasta "$col1"_hifi_methCalled.bam > "$col1"_hifi_methCalled_aligned.bam
	fi
	if [ "$partial_start" -le 5 ] && [ "$partial_end" -ge 5 ]; then
		[ -f "$col1"_hifi_methCalled_aligned.bam ] || { echo "Error: Aligned BAM file for ${col1} not found"; exit 1; }
		samtools sort -m "$maxmem" "$col1"_hifi_methCalled_aligned.bam -o "$col1"_hifi_methCalled_aligned_sorted.bam
	fi
	if [ "$partial_start" -le 6 ] && [ "$partial_end" -ge 6 ]; then
		[ -f "$col1"_hifi_methCalled_aligned_sorted.bam ] || { echo "Error: Sorted BAM file for ${col1} not found"; exit 1; }
		samtools index "$col1"_hifi_methCalled_aligned_sorted.bam "$col1"_hifi_methCalled_aligned_sorted.bam.bai
	fi
	if [ "$partial_start" -le 7 ] && [ "$partial_end" -ge 7 ]; then
		[ -f "$col1"_hifi_methCalled_aligned_sorted.bam ] || { echo "Error: Sorted BAM file for ${col1} not found"; exit 1; }
		aligned_bam_to_cpg_scores --bam "$col1"_hifi_methCalled_aligned_sorted.bam --output-prefix "$col1" --modsites-mode reference --ref *.fasta
		gzip -df *.bed.gz
		awk -v FS="\t" -v OFS="\t" '$1!~"^#"{print $1, $2, $3, $4}' *combined.bed > CG.mp
	fi
	if [ "$partial_start" -le 8 ] && [ "$partial_end" -ge 8 ]; then
		bedtools intersect -wo -a CG.mp -b <(awk -v type_string="$type_string" -v FS="\t" -v OFS="\t" '{if($7~/^\+$/&&$3==type_string){print $1, $4, $7} else if($7~/^\-$/&&$3==type_string){print $1, $5, $7}}' *.gtf | sort -n | uniq | awk -v FS="\t" -v OFS="\t" '{print $1, $2-10000, $2+10000, $3}' | awk -v FS="\t" -v OFS="\t" '{if($2<0){print $1, 0, $3, $2+10000, $4}else{print $1, $2, $3, $2+10000, $4}}') | awk -v FS="\t" -v OFS="" '{print $5, ".", $8, "\t", $2-$8, "\t", $4, "\t", $9}' > MP_TSS_"$col1".tsv
		bedtools intersect -wo -a CG.mp -b <(awk -v type_string="$type_string" -v FS="\t" -v OFS="\t" '{if($7~/^\+$/&&$3==type_string){print $1, $5, $7} else if($7~/^\-$/&&$3==type_string){print $1, $4, $7}}' *.gtf | sort -n | uniq | awk -v FS="\t" -v OFS="\t" '{print $1, $2-10000, $2+10000, $3}' | awk -v FS="\t" -v OFS="\t" '{if($2<0){print $1, 0, $3, $2+10000, $4}else{print $1, $2, $3, $2+10000, $4}}') | awk -v FS="\t" -v OFS="" '{print $5, ".", $8, "\t", $2-$8, "\t", $4, "\t", $9}' > MP_TES_"$col1".tsv
		sed -i '/^NC/!d' MP_TSS_"$col1".tsv
		sed -i '/^NC/!d' MP_TES_"$col1".tsv
	fi
	if [ "$partial_start" -le 9 ] && [ "$partial_end" -ge 9 ]; then
		awk -v type_string="$type_string" -v FS="\t" -v OFS="\t" '{if($7~/^\+$/&&$3==type_string){print $1, $4, $7} else if($7~/^\-$/&&$3==type_string){print $1, $5, $7}}' *.gtf | sort | uniq > TSS.txt
		awk -v type_string="$type_string" -v FS="\t" -v OFS="\t" '{if($7~/^\+$/&&$3==type_string){print $1, $5, $7} else if($7~/^\-$/&&$3==type_string){print $1, $4, $7}}' *.gtf | sort | uniq > TES.txt
		sed -i '/^NC/!d' TSS.txt
		sed -i '/^NC/!d' TES.txt
		streamgc *.fasta TSS.txt GC_TSS_"$col1".tsv
		streamgc *.fasta TES.txt GC_TES_"$col1".tsv
	fi
	if [ "$partial_start" -le 10 ] && [ "$partial_end" -ge 10 ]; then
		streambpb -b CC -w "$window_size" -i "$interval_size" *.fasta TSS.txt CpC_TSS_"$col1".tsv
		streambpb -b CC -w "$window_size" -i "$interval_size" *.fasta TES.txt CpC_TES_"$col1".tsv
		streambpb -b CG -w "$window_size" -i "$interval_size" *.fasta TSS.txt CpG_TSS_"$col1".tsv
		streambpb -b CG -w "$window_size" -i "$interval_size" *.fasta TES.txt CpG_TES_"$col1".tsv
		streambpb -b GC -w "$window_size" -i "$interval_size" *.fasta TSS.txt GpC_TSS_"$col1".tsv
		streambpb -b GC -w "$window_size" -i "$interval_size" *.fasta TES.txt GpC_TES_"$col1".tsv
		streambpb -b GG -w "$window_size" -i "$interval_size" *.fasta TSS.txt GpG_TSS_"$col1".tsv
		streambpb -b GG -w "$window_size" -i "$interval_size" *.fasta TES.txt GpG_TES_"$col1".tsv
	fi
 	cd ..
done
