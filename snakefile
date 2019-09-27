""" SCVUC v3 COI mtDNA metabarcode pipeline

Metabarcode pipeline to process Illumina MiSeq reads as follows:

1_Pair reads
2_Trim reads
3_Concatenate samples for global analysis
4_Dereplicate reads
5_Get CDS
6_Create denoised ESVs
7_Create ESV table
8_Taxonomic assignment

Note:
If you get a memory error at the USEARCH unoise3 denoising step, then use the COI alternate metabarcode pipeline instead.

"""

#######################################################################
# Read in vars from config.yaml
# Generate new vars, dirs, and filenames
# All code not in rules (this section) is run first before any rules (next section onwards)

# 1_Pair reads
raw_sample_read_wildcards=config["raw_sample_read_wildcards"]
print(raw_sample_read_wildcards)

SAMPLES,READS = glob_wildcards(raw_sample_read_wildcards)
SAMPLES_UNIQUE = list(set(SAMPLES))

raw_sample_forward_wildcard=config["raw_sample_forward_wildcard"]

dir=config["dir"]
stats_out="{sample}.stats"
raw1_stats_out=dir+"/stats/R1/"+stats_out
cat_raw1_stats_out=dir+"/stats/R1.stats"

raw_sample_reverse_wildcard=config["raw_sample_reverse_wildcard"]
raw2_stats_out=dir+"/stats/R2/"+stats_out
cat_raw2_stats_out=dir+"/stats/R2.stats"

#dir=config["dir"]
fastq_gz="{sample}.fastq.gz"
seqprep_out=dir+"/paired/"+fastq_gz
paired_stats_out=dir+"/stats/paired/"+stats_out
cat_paired_stats_out=dir+"/stats/paired.stats"

# 2_Trim reads
cutadapt_f_out=dir+"/Fprimer_trimmed/"+fastq_gz
Ftrimmed_stats_out=dir+"/stats/Ftrimmed/"+stats_out
cat_Ftrimmed_stats_out=dir+"/stats/Ftrimmed.stats"

fasta_gz="{sample}.fasta.gz"
cutadapt_r_out=dir+"/Rprimer_trimmed/"+fasta_gz
Rtrimmed_stats_out=dir+"/stats/Rtrimmed/"+stats_out
cat_Rtrimmed_stats_out=dir+"/stats/Rtrimmed.stats"

# 3_Concatenate samples for global analysis
fasta="{sample}.fasta"
concatenate_pattern=dir+"/"+fasta
output=dir+"/cat.fasta.tmp"
output2=dir+"/cat.fasta"
gzip_out=dir+"/cat.fasta.gz"

# 4_Dereplicate reads
vsearch_out=dir+"/cat.uniques"
vsearch_log=dir+"/dereplication.log"

# 5_Create denoised ESVs
usearch_out=dir+"/cat.denoised"
usearch_log=dir+"/usearch.log"

# 6_Get CDS
orf_out=dir+"/cds.fasta.tmp"
cds_out=dir+"/cds.fasta.tmp2"
range=dir+"/limits.csv"
plot=dir+"/cds.pdf"
cds_out2=dir+"/cds.fasta"

# 7_Create ESV table
vsearch_out2=dir+"/ESV.table"
vsearch_log2=dir+"/table.log"

# 8_Taxonomic assignment
rdp_out=dir+"/rdp.out.tmp"
rdp_csv=dir+"/rdp.csv.tmp"
rdp_csv2=dir+"/rdp.csv"

#######################################################################
# This rule defines the final target file that is needed from the pipeline
# Hashed lines (begin with # symbol) are not executed
# Target files need to be separated with a comma
# Snakemake looks at all input and output files to determine which rules need to be run to create the target file

rule all:
	input:
		# Rule testing [In order of execution]:
		# Calculate raw stats
#		expand(raw1_stats_out, sample=SAMPLES_UNIQUE),
#		expand(raw2_stats_out, sample=SAMPLES_UNIQUE),
		cat_raw1_stats_out,
		cat_raw2_stats_out,
		# 1_Pair reads
#		expand(seqprep_out, sample=SAMPLES)
		# Calculate paired stats
		expand(paired_stats_out, sample=SAMPLES_UNIQUE),
		cat_paired_stats_out,
		# 2_Trim reads (forward)
#		expand(cutadapt_f_out, sample=SAMPLES)
		# Calculate forward trimmed stats
		expand(Ftrimmed_stats_out, sample=SAMPLES_UNIQUE),
		cat_Ftrimmed_stats_out,
		# 2_Trim reads (reverse)
#		expand(cutadapt_r_out, sample=SAMPLES)
		# Calculate reverse trimmed stats
		expand(Rtrimmed_stats_out, sample=SAMPLES_UNIQUE),
		cat_Rtrimmed_stats_out,
		# 2_Trim reads (edit fasta header)
#		expand(concatenate_pattern, sample=SAMPLES)
		# 3_Concatenate samples for global analysis
#		output
		# 3_Concatenate samples for global analysis (edit fasta header)
#		output2
		# 3_Concatenate samples for global analysis (compress file)
#		gzip_out
		# 4_Dereplicate reads
#		vsearch_out
		# 5_Get CDS, translate ESVs into every possible ORF
#		orf_out
		# 5_Get CDS, filter CDS, keep longest
#		cds_out
		# 5_Get CDS, calculate range of CDS lengths
#		plot,
#		range,
		# 5_Get CDS, filter CDS, use range cutoffs
#		cds_out2
		# 6_Denoise reads
#		usearch_out
		# 7_Create ESV table
#		vsearch_out2
		# 8_Taxonomic assignment
#		rdp_out
		# 8_Taxonomic assignment (add read numbers)
#		rdp_csv
		# 8_Taxonomic assignment (edit ESV id's to include amplicon name) [Final output file]
		rdp_csv2
	 
#######################################################################
# Count raw forward reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > R1.stats'

rule raw1_stats:
	input: 
		raw_sample_forward_wildcard
	output:
		raw1_stats_out
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate R1 reads

rule cat_raw1_stats:
	input:
		expand(raw1_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_raw1_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Count raw reverse reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > R2.stats'

rule raw2_stats:
	input: 
		raw_sample_reverse_wildcard
	output:
		raw2_stats_out
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate R2 reads

rule cat_raw2_stats:
	input:
		expand(raw2_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_raw2_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Pair forward and reverse reads with SeqPrep

rule pair_reads:
	version: "1.3.2"
	input:
		f=raw_sample_forward_wildcard,
		r=raw_sample_reverse_wildcard
	output:
		X1=temp("{sample}_R1.out"),
		X2=temp("{sample}_R2.out"),
		s=seqprep_out
	shell:
		"SeqPrep -f {input.f} -r {input.r} -1 {output.X1} -2 {output.X2} -q {config[SEQPREP][q]} -s {output.s} -o {config[SEQPREP][o]}"

#######################################################################
# Count paired reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > paired.stats'

rule paired_stats:
	input: 
		seqprep_out
	output:
		paired_stats_out
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate paired reads

rule cat_paired_stats:
	input:
		expand(paired_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_paired_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Trim forward primer with CUTADAPT

rule trim_forward_primer:
	version: "2.4"
	input:
		seqprep_out
	output:
		cutadapt_f_out
	shell:
		"cutadapt -g {config[CUTADAPT_FWD][g]} -m {config[CUTADAPT_FWD][m]} -q {config[CUTADAPT_FWD][q]} --max-n={config[CUTADAPT_FWD][mn]} --discard-untrimmed -o {output} {input}"

#######################################################################
# Count forward primer trimmed reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > Ftrimmed.stats'

rule Ftrimmed_stats:
	input: 
		cutadapt_f_out
	output:
		Ftrimmed_stats_out
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate Ftrimmed reads

rule cat_Ftrimmed_stats:
	input:
		expand(Ftrimmed_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_Ftrimmed_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Trim reverse primer with CUTADAPT

rule trim_reverse_primer:
	version: "2.4"
	input:
		cutadapt_f_out
	output:
		cutadapt_r_out
	shell:
		"cutadapt -a {config[CUTADAPT_REV][a]} -m {config[CUTADAPT_REV][m]} -q {config[CUTADAPT_REV][q]} --max-n={config[CUTADAPT_REV][mn]} --discard-untrimmed -o {output} {input}"

#######################################################################
# Count reverse primer trimmed reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > Rtrimmed.stats'

rule Rtrimmed_stats:
	input: 
		cutadapt_r_out
	output:
		Rtrimmed_stats_out
	shell:
		"perl perl_scripts/fasta_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate Rtrimmed reads

rule cat_Rtrimmed_stats:
	input:
		expand(Rtrimmed_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_Rtrimmed_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Edit fasta header with Perl script

rule edit_fasta_header1:
	input:
		cutadapt_r_out
	output:
		concatenate_pattern
	threads: 1
	shell:
		"perl perl_scripts/rename_fasta_gzip.plx {input} > {output}"

#######################################################################
# Concatenate all samples into single file for global analysis using Linux cat
# Only use a single thread to work right

rule concatenate:
	input:
		expand(concatenate_pattern, sample=SAMPLES)
	output:
		temp(output)
	threads: 1
	shell:
		"cat {input} >> {output}"

#######################################################################
# Edit fasta header again using GNU sed

rule edit_fasta_header2:
	version: "4.7"
	input:
		output
	output:
		temp(output2)
	shell:
		"sed -e 's/-/_/g' {input} >> {output}"

#######################################################################
# Compress into smaller file using Linux gzip

rule compress:
	version: "1.7"
	input:
		output2
	output:
		gzip_out
	shell:
		"gzip -c {input} > {output}"

#######################################################################
# Dereplicate and track reads with VSEARCH

rule dereplicate:
	version: "2.13.6"
	input:
		gzip_out
	output:
		vsearch_out
	log: vsearch_log
	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizein --sizeout --log {log}"

#######################################################################
# Denoise with USEARCH
# Make sure this is installed locally and in your PATH
# I changed the default name of the program to 'usearch11' to be more concise
# redirect stderr to stdout(&1) so it all gets printed out into a log file

rule denoise:
	version: "11.0.667_i86linux32"
	input:
		vsearch_out
	output:
		usearch_out
	log: usearch_log
	shell:
		"usearch11 -unoise3 {input} -zotus {output} -minsize {config[USEARCH][minsize]} > {log} 2>&1"


#######################################################################
# Get CDS, first translate ESVs into every possible ORF

rule get_orfs:
	version: "0.4.3"
	input:
		usearch_out
	output:
		orf_out
	shell:
		"ORFfinder -in {input} -g {config[ORFFINDER][g]} -s {config[ORFFINDER][s]} -out {output} -outfmt {config[ORFFINDER][outfmt]}"

#######################################################################
# Get CDS, get longest CDS

rule filter_cds:
	input:
		orf_out
	output:
		cds_out
	shell:
		"perl perl_scripts/get_longest_orf_nt2.plx {input} > {output}"

#######################################################################
# Calculate outliers

rule find_outliers:
	input:
		cds_out
	output:
		out1=plot,
		out2=range
	shell:
		"Rscript R_scripts/length_hist.R {input} {output.out1} {output.out2}"

#######################################################################
# Get CDS, only keep CDS with lengths within 1.5*IQR

rule filter_cds2:
	input:
		orfs=orf_out,
		range=range
	output:
		cds_out2
	shell:
		"perl perl_scripts/get_longest_orf_nt3.plx {input.orfs} {input.range} > {output}"

#######################################################################
# Create ESV table

rule create_ESV_table:
	version: "2.13.16"
	input:
		usearch_global=gzip_out,
		db=cds_out2
	output:
		vsearch_out2
	threads: config["VSEARCH"]["t"]
	log: vsearch_log2
	shell:
		"vsearch --threads {threads} --usearch_global {input.usearch_global} --db {input.db} --id 1.0 --otutabout {output} --log {log}"

#######################################################################
# Taxonomic assignment

# Use the RDP classifier for taxonomic assignment
# Make sure this is installed locally (do not use the old 2.2 version available from conda), download the v2.12 from sourceforge, indicate the path to the classifeir in the config file
# Be sure to indicate the location of the trained reference database properties file in the config file

rule taxonomic_assignment:
	version: "2.12"
	input:
		cds_out2
	output:
		temp(rdp_out)
	shell:
		"java -Xmx8g -jar {config[RDP][jar]} classify -t {config[RDP][t]} -o {output} {input}"
	
#######################################################################
# Map read number

# Edit path in top line of Perl script to match path to perl in conda environment

rule map_read_number:
	input:
		table=vsearch_out2,
		rdp=rdp_out
	output:
		temp(rdp_csv)
	shell:
		"perl perl_scripts/add_abundance_to_rdp_out4.plx {input.table} {input.rdp} > {output}"

#######################################################################
# Edit csv to add amplicon name to Zotu
# This is needed when data from different amplicons are combined, so that the Zotu's remain unique

rule edit_csv:
	input:
		rdp_csv
	output:
		rdp_csv2
	shell:
		"sed -e {config[SED]} {input} > {output}"


# End of Metabarcode pipeline
