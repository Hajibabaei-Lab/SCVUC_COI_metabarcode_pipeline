# README

This repository outlines how COI metabarcodes are processed by Teresita M. Porter. **SCVUC** refers to the programs, algorithms, and reference datasets used in this data flow: **S**EQPREP, **C**UTADAPT, **V**SEARCH, **U**NOISE, **C**OI classifier. 

The pipeline begins with raw Illumina MiSeq fastq.gz files with paired-end reads. Reads are paired. Primers are trimmed. All the samples are pooled for a global analysis. Reads are dereplicated and denoised producing a reference set of exact sequence variants (ESVs). ESVs are filtered again by retaining the longest coding sequence so long as they exceed a minimum length cutoff.  These coding sequences are taxonomically assigned using the COI mtDNA reference set available from https://github.com/terrimporter/CO1Classifier and is used with the RDP Classifier (Wang et al., 2007) available from https://sourceforge.net/projects/rdp-classifier/ .

This data flow has been developed using a conda environment and snakemake pipeline for improved reproducibility. It will be updated on a regular basis so check for the latest version at https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline/releases .

## How to cite

If you use this dataflow or any of the provided scripts, consider citing the CO1 classifier paper (Porter & Hajibabaei, 2018 Sci Rep) if you use it and provide a link to this page https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline in your publication.

## Outline

[Standard pipeline](#standard-pipeline)  

[Implementation notes](#implementation-notes)  

[References](#references)  

[Acknowledgements](#acknowledgements)  

## Standard pipeline

### Overview of the standard pipeline

If you are comfortable reading code, read through the snakefile to see how the pipeline runs, and which programs and versions are used.

#### A brief overview:

Raw paired-end reads are merged using SEQPREP v1.3.2 from bioconda (St. John, 2016).  This step looks for a minimum Phred quality score of 20 in the overlap region, requires at least 25bp overlap.

Primers are trimmed in two steps using CUTADAPT v2.4 from bioconda (Martin, 2011).  This step looks for a minimum Phred quality score of at least 20 at the ends, forward primer is trimmed first, no more than 3 N's allowed, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.  The output from the first step, is used as in put for the second step.  This step looks for a minimum Phred quality score of at least 20 at the ends, the reverse primer is trimmed, no more than 3 N's allowed, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.

Files are reformatted and samples are combined for a global analysis.

Reads are dereplicated (only unique sequences are retained) using VSEARCH v2.13.6 from bioconda (Rognes et al., 2016).

Denoised exact sequence variants (ESVs) are generated using USEARCH v11.0.667 with the unoise3 algorithm (Edgar, 2016).  This step removes any PhiX contamination, putative chimeric sequences, sequences with predicted errors, and rare sequences.  This step produces zero-radius OTUs (Zotus) also referred to commonly as amplicon sequence variants (ASVs), ESVs, or 100% operational taxonomic unit (OTU) clusters.  Here, we define rare sequences to be sequence clusters containing only one or two reads (singletons and doubletons) and these are removed as 'noise'.

The ESVs are translated into every possible open reading frame and the longest coding sequence is retained so long as they exceed a minimum length cutoff.  The cutoff is determined empirically to screen out the most obvious pseudogenenes that may have a shorter than expected length due to indels and frameshifts (see [Implementation notes](#implementation-notes) below).

An ESV table that tracks read number for each coding sequence in each sample is generated with VSEARCH.

COI mtDNA taxonomic assignments are made using the Ribosomal Database classifier v2.12 (RDP classifier) available from https://sourceforge.net/projects/rdp-classifier/ (Wang et al., 2007) using the COI classifier v3.2 reference dataset (Porter and Hajibabaei, 2018 Sci Rep).

The final output is reformatted to add read numbers for each sample and column headers to improve readability.

Read and ESV statistics are provided for various steps of the program are also provided.

### Prepare your environment to run the pipeline

1. This pipeline includes a conda environment that provides most of the programs needed to run this pipeline (SNAKEMAKE, SEQPREP, CUTADAPT, VSEARCH, etc.).

```linux
# Create the environment from the provided environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate myenv
```
2. The pipeline requires commercial software for the denoising step.  A free 32-bit version of USEARCH v11.0.667 can be obtained from https://drive5.com/usearch/download.html .  Be sure to put the program in your PATH, ex. ~/bin .  Make it executable and rename it to simply usearch11.

```linux
mv usearch11.0.667_i86linux32 ~/bin/.
cd ~/bin
chmod 755 usearch11.0.667_i86linux32
mv usearch11.0.667_i86linux32 usearch11
```

3. The pipeline requires ORFfinder 0.4.3 available from the NCBI at ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ .  This program should be downloaded, made executable, and put in your path.

```linux
# download
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz

# decompress
gunzip ORFfinder.gz

# make executable
chmod 755 ORFfinder

# put in your PATH (ex. ~/bin)
mv ORFfinder ~/bin/.
```

4. The pipeline also requires the RDP classifier for the taxonomic assignment step.  Although the RDP classifier v2.2 is available through conda, a newer v2.12 is available form SourceForge at https://sourceforge.net/projects/rdp-classifier/ .  Download it and take note of where the classifier.jar file is as this needs to be added to config.yaml .

The RDP classifier comes with the training sets to classify 16S, fungal ITS and fungal LSU rDNA sequences.  To classify COI mtDNA sequences, obtain the COI classifier v4 reference set from GitHub 
https://github.com/terrimporter/CO1Classifier .  Take note of where the rRNAclassifier.properties file is as this needs to be added to the config.yaml .

```linux
RDP:
    jar: "/path/to/rdp_classifier_2.12/dist/classifier.jar"
    t: "/path/to/CO1Classifier/v4/NCBI_BOLD_merged/mydata/mydata_trained/rRNAClassifier.properties"
```

4. In most cases, your raw paired-end Illumina reads can go into a directory called 'data' which should be placed in the same directory as the other files that come with this pipeline.

```linux
# Create a new directory to hold your raw data
mkdir data
```

5. Please go through the config.yaml file and edit directory names, filename patterns, etc. as necessary to work with your filenames.

6. Be sure to edit the first line of each Perl script (shebang) in the perl_scripts directory to point to where Perl is installed.

```linux
# The usual shebang if you already have Perl installed
#!/usr/bin/perl

# Alternate shebang if you want to run perl using the conda environment (edit this)
#!/path/to/miniconda3/envs/myenv/bin/perl
```

### Run the standard pipeline

Run snakemake by indicating the number of jobs or cores that are available to run the whole pipeline.  

```linux
# Choose and/or edit the appropriate configuration file (ex. config_BR5.yaml, config_F230R.yaml)
snakemake --jobs 24 --snakefile snakefile --configfile config_BR5.yaml
```

You can view read number and length (min, max, mean, median, mode) statistics for each sample at steps of the bioinformatic pipeline.  A simple report can be generated like so, modify to summarize reports for different bioinformatic steps (raw reads, paired reads, primer trimmed reads):

```linux
# Generate a report for raw R1 reads
cd stats/raw/R1
cat *.stats > R1.stats
```

When you are done, deactivate the conda environment:

```linux
conda deactivate
```

## Implementation notes

### Changing the minimum coding sequence cutoff value

The minimum coding sequence cutoff value chosen dependson the amplicon being analyzed and is chosen empirically.  It is a good idea to test a suite of cutoff values and select a cutoff based on the results.  

1. You can easily change this value in the config.yaml file.

```linux
min_config_length: 309
```

2. Edit the snakefile to stop after generating the cds.fasta file by removing the hash (#) symbol before cds_out and adding a hash symbol to the rdp_csv2 line:

```linux
rule all:
    input:
        ...
		# 5_Get CDS, filter CDS, keep longest
		cds_out
        ...
        # 8_Taxonomic assignment (edit ESV id's to include amplicon name) [Final output file]
#		rdp_csv2
```

3. Run the edited snakefile:

```linux
snakemake --jobs 24 --snakefile snakefile --configfile config.yaml
```

4. Rename the the default outfile:

```linux
mv BR5/cds.fasta BR5/cds.fasta.309
```

Repeat steps 1-4 with new cutoff values.  Compare the number of retained coding sequences in the cds.fasta files (ex. plot cutoff versus number of cds and choose a value before a steep dropoff in the number of retained cds).  Alternatively, use the average/mode length of the primer trimmed sequences as the cutoff value, plot a histogram of the resulting cds lengths, choose a cutoff based on the results.

### Batch renaming of files

Sometimes it necessary to rename raw data files in batches.  I use Perl-rename (Gergely, 2018) that is available at https://github.com/subogero/rename not linux rename.  I prefer the Perl implementation so that you can easily use regular expressions.  I first run the command with the -n flag so you can review the changes without making any actual changes.  If you're happy with the results, re-run without the -n flag.

```linux
rename -n 's/PATTERN/NEW PATTERN/g' *.gz
```

### Symbolic links

Symbolic links are like shortcuts or aliases that can also be placed in your ~/bin directory that point to files or programs that reside elsewhere on your system.  So long as those scripts are executable (e.x. chmod 755 script.plx) then the shortcut will also be executable without having to type out the complete path or copy and pasting the script into the current directory.

```linux
ln -s /path/to/target/directory shortcutName
ln -s /path/to/target/directory fileName
ln -s /path/to/script/script.sh commandName
```

## References

Edgar, R. C. (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. BioRxiv. doi:10.1101/081257  
Gergely, S. (2018, January). Perl-rename. Retrieved from https://github.com/subogero/rename  
Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. Journal, 17(1), pp–10.  
Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.  
Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. doi:10.7717/peerj.2584  
St. John, J. (2016, Downloaded). SeqPrep. Retrieved from https://github.com/jstjohn/SeqPrep/releases  
Tange, O. (2011). GNU Parallel - The Command-Line Power Tool. ;;Login: The USENIX Magazine, February, 42–47.  
Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07  

## Acknowledgements

I would like to acknowedge funding from the Canadian government through the Genomics Research and Development Initiative (GRDI) EcoBiomics project.

Last updated: Sept. 6, 2019
