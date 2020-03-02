# README

This repository contains a conda environemnt and Snakemake pipeline to processes Illumina paired-end COI metabarcodes (Koster and Rahmann, 2012). **SCVUC** refers to the programs, algorithms, and reference datasets used in this data flow: **S**EQPREP, **C**UTADAPT, **V**SEARCH, **U**NOISE, **C**OI classifier. 

## Overview

The pipeline begins with raw paired-end Illumina MiSeq fastq.gz files.  Reads are paired.  Primers are trimmed.  All the samples are pooled for a global analysis.  Reads are dereplicated and denoised producing a reference set of exact sequence variants (ESVs). The ESVs are taxonomically assigned using the COI mtDNA reference set available from https://github.com/terrimporter/CO1Classifier and is used with the RDP Classifier (Wang et al., 2007) available from https://sourceforge.net/projects/rdp-classifier/ .  At this point, the final output file is prepared OR one can proceed to filter the dataset to focus on arthropods and attempt to remove apparent pseudogenes.  The final output file contains taxonomic assignments and number of reads for each ESV and sample.  If pseudogenes were filtered, then the final output file will also contain the sequence of the longest retained open reading frame for all Arthropoda.

This data flow will be updated on a regular basis so check for the latest version at https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline/releases .

## How to cite

If you use this dataflow or any of the provided scripts, consider citing the CO1 classifier paper (Porter & Hajibabaei, 2018 Sci Rep) if you use it and provide a link to this page https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline in your publication.

## Outline

[Pipeline details](#pipeline-details)  

[Prepare your environment to run the pipeline](#prepare-your-environment-to-run-the-pipeline)  

[Run the pipeline](#run-the-pipeline)  

[Implementation notes](#implementation-notes)  

[References](#references)  

[Acknowledgements](#acknowledgements)  

## Pipeline details

If you are comfortable reading code, read through the snakefile to see how the pipeline runs as well as which programs and versions are used.  Otherwise you can just list all the programs in the conda environment, see [Implementation notes](#implementation-notes).  

Raw paired-end reads are merged using SEQPREP v1.3.2 from bioconda (St. John, 2016).  This step looks for a minimum Phred quality score of 20 in the overlap region, requires at least a 25bp overlap.

Primers are trimmed in two steps using CUTADAPT v2.6 from bioconda (Martin, 2011).  This step looks for a minimum Phred quality score of at least 20 at the ends, the forward primer is trimmed first based on its sequence, no more than 3 N's allowed, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.  The output from the first step, is used as input for the second step.  This step looks for a minimum Phred quality score of at least 20 at the ends, the reverse primer is trimmed based on its sequence, no more than 3 N's allowed, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.

Files are reformatted and samples are combined for a global analysis.

Reads are dereplicated (only unique sequences are retained) using VSEARCH v2.14.1 from bioconda (Rognes et al., 2016).

Denoised exact sequence variants (ESVs) are generated using VSEARCH with the unoise3 algorithm (Edgar, 2016).  This step removes any PhiX contamination, sequences with predicted errors, and rare sequences.  This step also produces zero-radius OTUs (Zotus) also referred to commonly as amplicon sequence variants (ASVs), ESVs, or 100% operational taxonomic unit (OTU) clusters.  Here, we define rare sequences to be sequence clusters containing only one or two reads (singletons and doubletons) and these are removed as 'noise'.  Putative chimeric sequences are then removed using the uchime3_denovo algorithm in VSEARCH.

An ESV x sample table that tracks read number for each ESV (longest ORF) is generated with VSEARCH.

COI mtDNA taxonomic assignments are made using the Ribosomal Database classifier v2.12 (RDP classifier) available from https://sourceforge.net/projects/rdp-classifier/ (Wang et al., 2007) using the COI classifier v4 reference dataset available from https://github.com/terrimporter/CO1Classifier (Porter and Hajibabaei, 2018 Sci Rep).

If you use the pipeline that attempts to remove arthropod pseudogenes, then arthropod ESVs are translated into every possible open reading frame (ORF) on the plus strand.  The longest ORFs are retained for each ESV and ESVs with outlier lengths are removed, i.e., ORFs that are too short or too long are removed.  Outliers are identified as ORFs with lengths outside the range of the 25th percentile - 1.5\*IQR and the 75th percentile + 1.5\*IQR (IQR, inter quartile range).  This method should help to screen out the most obvious pseudogenes that may have a shorter than expected length due to premature stop codons introduced through frameshifts caused by sequence errors, indels, or longer than expected lengths due to insertions.  There is no guarantee that genuine coding sequences are not erroneously removed during this step.  If your dataset contains taxa with coding sequences known to be unusually shorter or longer than usual, then the alternate pipeline should be used.

The final output is reformatted to add read numbers for each sample and column headers to improve readability.  If you ran the pipeline that attempts to remove arthropod pseudogenes, then the sequene for the longest retained open reading frame is also provided.

Statistics and log files are also provided for each major bioinformatic step.

## Prepare your environment to run the pipeline

1. This pipeline includes a conda environment that provides most of the programs needed to run this pipeline (SNAKEMAKE, SEQPREP, CUTADAPT, VSEARCH, etc.).

```linux
# Create the environment from the provided environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate SCVUCv4.3

# On the GPSC activate using source
source ~/miniconda/bin/activate SCVUCv4.3
```

2. The pipeline requires ORFfinder 0.4.3 available from the NCBI at ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ .  This program should be downloaded, made executable, and put in your path.

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

Run the program to test that it works:

```linux
ORFfinder
```

If you get an error that requries a GLIBC_2.14 libc.so.6 library, then follow the instructions at [Use conda's libc library for NCBI's ORFfinder](#use-condas-libc-library-for-ncbis-orffinder).

3. The pipeline also requires the RDP classifier for the taxonomic assignment step.  Although the RDP classifier v2.2 is available through conda, a newer v2.12 is available form SourceForge at https://sourceforge.net/projects/rdp-classifier/ .  Download it and take note of where the classifier.jar file is as this needs to be added to config.yaml .

The RDP classifier comes with the training sets to classify 16S, fungal ITS, and fungal LSU rDNA sequences.  To classify COI mtDNA sequences, obtain the COI classifier v4 reference set from GitHub 
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
#!/path/to/miniconda3/envs/SCVUCv4.3/bin/perl
```

## Run the pipeline

Run Snakemake by indicating the number of jobs or cores that are available to run the whole pipeline.  

```linux
# command to run the pipeline without attempting to remove any pseudogenes
# be sure to choose the appropriate config file for your amplicon, or edit to match your primers (eg. BR5 shown here)
snakemake --jobs 24 --snakefile snakefile_withoutPseudogeneFiltering --config_BR5.yaml

# command to run the pipeline to remove arthropod pseudogenes
snakemake --jobs 24 --snakefile snakefile_withArthropodaPseudogeneFiltering --config_BR5.yaml
```

You can view read number and length (min, max, mean, median, mode) statistics for each sample at steps of the bioinformatic pipeline.  

When you are done, deactivate the conda environment:

```linux
conda deactivate
```

## Implementation notes

### Installing Conda and Snakemake

Conda is an open source package and environment management system.  Miniconda is a lightweight version of conda that only contains conda, python, and their dependencies.  Using conda and the environment.yml file provided here can help get all the necessary programs in one place to run this pipeline.  Snakemake is a Python-based workflow management tool meant to define the rules for running this bioinformatic pipeline.  There is usually no need to edit the snakefile file directly.  Changes to select parameters can be made in the config.yaml file.  If you install conda and activate the environment provided, then you will also get the correct versions of the open source programs used in this pipeline including Snakemake.

Install miniconda as follows:

```linux
# Download miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install miniconda3 and initialize
sh Miniconda3-latest-Linux-x86_64.sh

# Add conda to your PATH, ex. to ~/bin
cd ~/bin
ln -s ~/miniconda3/bin/conda conda

# Activate conda method 1 (working in a container)
source ~/miniconda3/bin/activate SCVUCv4.3

# Activate conda method 2
conda activate SCVUCv4.3
```

### Check program versions

Ensure the program versions in the environment are being used.

```linux
# create conda environment from file
conda env create -f environment.yml

# activate the environment
conda activate SCVUCv4.3

# list all programs available in the environment at once
conda list > programs.list

# you can check that key programs in the conda environment are being used (not locally installed versions)
which SeqPrep
which cutadapt
which vsearch
which perl

# you can also check their version numbers one at a time instead of running 'conda list'
cutadapt --version
vsearch --version
```

Version numbers are also tracked in the snakefile.

### Use conda's libc library for NCBI's ORFfinder

The glibc 2.14 library is already available in the SCVUCv4.3 environment.  The LD_LIBRARY_PATH environment variable will need to be activated (and deactivated) by adding the following scripts as follows:

Create the shell script file LD_PATH.sh in the following location to set the environment variable: ~/miniconda3/envs/SCVUCv4.3/etc/conda/activate.d/LD_PATH.sh

Put the following text in the LD_PATH.sh file:

```linux
export LD_LIBRARY_PATH_CONDA_BACKUP=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib64:$LD_LIBRARY_PATH
```

Create the file LD_PATH.sh in the following location to unset the environment variable: ~/miniconda3/envs/SCVUCv4.3/etc/conda/deactivate.d/LD_PATH.sh

Put the following text in the LC_PATH.sh file:

```linux
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CONDA_BACKUP
```

Deactivate then reactivate the environment.

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
Koster, J., Rahmann, S. (2012) Snakemake - a scalable bioinformatics workflow engine.  Bioinformatics, 29(19): 2520-2522.
Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. Journal, 17(1), pp–10.  
Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.  
Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. doi:10.7717/peerj.2584  
St. John, J. (2016, Downloaded). SeqPrep. Retrieved from https://github.com/jstjohn/SeqPrep/releases  
Tange, O. (2011). GNU Parallel - The Command-Line Power Tool. ;;Login: The USENIX Magazine, February, 42–47.  
Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07  

## Acknowledgements

I would like to acknowedge funding from the Canadian government through the Genomics Research and Development Initiative (GRDI) EcoBiomics project.

Last updated: March 2, 2020
