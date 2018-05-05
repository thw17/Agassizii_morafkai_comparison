configfile: "gopherus_config.json"

angsd_path = "/home/thwebste/Tools/angsd/angsd"
bbduksh_path = "bbduk.sh"
bedtools_path = "bedtools"
bwa_path = "bwa"
fastqc_path = "fastqc"
gatk_path = "gatk"
multiqc_path = "multiqc"
ngscovar_path = "/home/thwebste/Tools/ngsTools/ngsPopGen/ngsCovar"
samblaster_path = "samblaster"
samtools_path = "samtools"

temp_directory = "temp"
fastq_directory = "/scratch/thwebste/Gopherus_popgen/fastqs"

gmor = config["morafkai"]
gaga = [x for x in config["samples"] if x not in config["morafkai"]]

sample_dict = {"gmor": gmor, "gaga": gaga, "all": config["samples"]}

fastq_prefixes = [
	config[x]["fq1"][:-9] for x in config["samples"]] + [
		config[x]["fq2"][:-9] for x in config["samples"]]

gaga_chunks = {
	"agassizii_1": gaga[:50],
	"agassizii_2": gaga[50:100],
	"agassizii_3": gaga[100:150],
	"agassizii_4": gaga[150:200],
	"agassizii_5": gaga[200:250],
	"agassizii_6": gaga[250:]}

gaga_chunks_list = [
	"agassizii_1", "agassizii_2",
	"agassizii_3", "agassizii_4",
	"agassizii_5", "agassizii_6"]

scaffolds_to_analyze = ["ScCC6lQ_16796", "ScCC6lQ_161298"]

# Reference chunks
num_chunks = 15
chunk_range = [x for x in range(1, num_chunks + 1)]

rule all:
	input:
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html",
		expand(
			"stats/{sample}.{genome}.mkdup.sorted.bam.stats",
			sample=config["samples"], genome=["gopaga20"]),
		expand(
			"gvcf_databases/{comparison}-{genome}-{chrom}",
			comparison=["gmor", "gaga", "all"],
			genome=["gopaga20"],
			chrom=config["scaffolds_no_semi_colon"])

rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "new_reference/{assembly}.fasta",
		fai = "new_reference/{assembly}.fasta.fai",
		amb = "new_reference/{assembly}.fasta.amb",
		dict = "new_reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		shell(
			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell(
			"{params.samtools} faidx {output.new}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {output.new}")
		# bwa
		shell(
			"{params.bwa} index {output.new}")

rule chunk_reference:
	input:
		fai = "new_reference/{assembly}.fasta.fai"
	output:
		expand("new_reference/{{assembly}}_split_chunk{num}.bed", num=chunk_range)
	params:
		chunks = num_chunks,
		out_prefix = "new_reference/{assembly}_split"
	shell:
		"python scripts/Chunk_fai.py --fai {input.fai} "
		"--out_prefix {params.out_prefix} --chunks {params.chunks}"

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc fastqc"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	threads:
		2
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref=adapters/{wildcards.sample}.adapters ktrim=r k=21 mink=11 hdist=2 tbo tpe qtrim=rl trimq=10"

rule fastqc_analysis_trimmed:
	input:
		"trimmed_fastqs/{sample}_trimmed_{read}.fastq.gz"
	output:
		"fastqc_trimmed/{sample}_trimmed_{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed {input}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed/{sample}_trimmed_{read}_fastqc.html",
			sample=config["samples"], read=["read1", "read2"])
	output:
		"multiqc_trimmed/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc_trimmed fastqc_trimmed"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		amb = "new_reference/{genome}.fasta.amb",
		ref = "new_reference/{genome}.fasta"
	output:
		"processed_bams/{sample}.{genome}.mkdup.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samblaster = samblaster_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	shell:
		"{params.bwa} mem -t {params.threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samblaster} | {params.samtools} fixmate -O bam - - "
		"| {params.samtools} sort -O bam -o {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.mkdup.sorted.bam"
	output:
		"processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	output:
		"stats/{sample}.{genome}.mkdup.sorted.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule create_interval_file_for_genomicsdbimport:
	input:
		fai = "new_reference/{genome}.fasta.fai"
	output:
		temp("interval_files/{chrom}_{genome}.bed")
	params:
		chromosome = lambda wildcards: config["scaffold_dict"][wildcards.chrom]
	shell:
		"""python scripts/Create_scaffold_bed_from_fai.py --fai {input.fai} """
		"""--outfile {output} --scaffold_name '{params.chromosome}'"""

rule gatk_gvcf_per_chrom:
	input:
		ref = "new_reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
		interval = "interval_files/{chrom}_{genome}.bed"
	output:
		"gvcfs/{sample}.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		chromosome = lambda wildcards: config["scaffold_dict"][wildcards.chrom]
	threads:
		4
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.interval} """
		"""-ERC GVCF -O {output}"""

rule gatk_genomicsdbimport_per_chrom:
	input:
		gvcfs = lambda wildcards: expand(
			"gvcfs/{sample}.{genome}.{chrom}.g.vcf.gz",
			sample=sample_dict[wildcards.comparison],
			genome=[wildcards.genome],
			chrom=[wildcards.chrom]),
		ref = "new_reference/{genome}.fasta",
		interval = "interval_files/{chrom}_{genome}.bed"
	output:
		"gvcf_databases/{comparison}-{genome}-{chrom}"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		chromosome = lambda wildcards: config["scaffold_dict"][wildcards.chrom]
	threads:
		8
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("-V " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx32g -Djava.io.tmpdir={params.temp_dir}" """
			"""GenomicsDBImport -R {input.ref} {variant_files} """
			"""--genomicsdb-workspace-path {output} -L {input.interval}""")


# rule gatk_gvcf_per_chunk:
# 	input:
# 		ref = "new_reference/{genome}.fasta",
# 		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
# 		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
# 		chunkfile = "new_reference/{genome}_split_chunk{chunk}.bed"
# 	output:
# 		"gvcfs/{sample}.{genome}.{chunk}.g.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path
# 	threads:
# 		4
# 	shell:
# 		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
# 		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.chunkfile} """
# 		"""-ERC GVCF -O {output}"""
#
# rule gatk_combinegvcfs_per_chunk:
# 	input:
# 		ref = "new_reference/{assembly}.fasta",
# 		gvcfs = lambda wildcards: expand(
# 			"gvcfs/{sample}.{assembly}.{chunk}.g.vcf.gz",
# 			sample=sample_dict[wildcards.comparison],
# 			assembly=[wildcards.assembly],
# 			chunk=[wildcards.chunk])
# 	output:
# 		"combined_gvcfs/{comparison}.{assembly}.{chunk}.gatk.combined.g.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path
# 	threads:
# 		8
# 	run:
# 		variant_files = []
# 		for i in input.gvcfs:
# 			variant_files.append("--variant " + i)
# 		variant_files = " ".join(variant_files)
# 		shell(
# 			"""{params.gatk} --java-options "-Xmx32g -Djava.io.tmpdir={params.temp_dir}" """
# 			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")
