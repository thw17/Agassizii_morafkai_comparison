configfile: "gopherus_config.json"

# These paths will need to be changed (also check reference path in gopherus_config.json)

gatk3_path = "GenomeAnalysisTK_37.jar"

temp_directory = "temp"
fastq_directory = "fastqs"

# Runtime values
very_short = "6:00:00"
medium = "12:00:00"
day = "24:00:00"
long = "48:00:00"

# These paths assume you've set up a Conda environment; change if not
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
fastqc_path = "fastqc"
gatk_path = "gatk"
mosdepth_path = "mosdepth"
multiqc_path = "multiqc"
picard_path = "picard"
samtools_path = "samtools"

# Organize sample lists
gmor = config["morafkai"]
gaga = config["agassizii"]
gber = config["berlandieri"]

sample_dict = {"gmor": gmor, "gaga": gaga, "gber": gber, "all": config["samples"]}

fastq_prefixes = [
	config[x]["fq1"][:-9] for x in config["samples"]] + [
		config[x]["fq2"][:-9] for x in config["samples"]]


# Set chromosomes to analyze
chroms_to_analyze = config["scaffolds_to_analyze"]

# Pipeline
rule all:
	input:
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html",
		expand(
			"stats/{sample}.{genome}.mkdup.sorted.realigned.bam.stats",
			sample=config["samples"], genome=["gopaga20"]),
		expand(
			"mosdepth_results/{sample}.{genome}.total.total.nodups_mapq30.mosdepth.summary.txt",
			sample=config["samples"], genome=["gopaga20"]),
		expand(
			"vcf_genotyped/{genome}.{chrom}.gatk.called.raw.vcf.gz",
			genome=["gopaga20"], chrom=chroms_to_analyze)
		# "angsd_results/all_angsd_CONCATENATED.covar"
		# expand(
		# 	"gvcf_databases/{comparison}-{genome}-{chrom}",
		# 	comparison=["gmor", "gaga", "all"],
		# 	genome=["gopaga20"],
		# 	chrom=config["scaffolds_no_semi_colon"]),
		# expand(
		# 	"final_vcfs/{comparison}.{genome}.combined.raw.vcf.gz",
		# 	comparison=["gmor", "gaga", "all"],
		# 	genome=["gopaga20"])

rule download_reference:
	output:
		"reference/{genome}.fasta"
	params:
		web_address = lambda wildcards: config["fasta_address"][wildcards.genome],
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"wget {params.web_address} -O {output}"

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path,
		threads = 4,
		mem = 16,
		t = medium
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

# rule chunk_reference:
# 	input:
# 		fai = "reference/{assembly}.fasta.fai"
# 	output:
# 		expand("reference/{{assembly}}_split_chunk{num}.bed", num=chunk_range)
# 	params:
# 		chunks = num_chunks,
# 		out_prefix = "reference/{assembly}_split"
# 	shell:
# 		"python scripts/Chunk_fai.py --fai {input.fai} "
# 		"--out_prefix {params.out_prefix} --chunks {params.chunks}"

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		bbduksh = bbduksh_path,
		threads = 2,
		mem = 8,
		t = very_short
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=10"

rule fastqc_analysis_trimmed:
	input:
		"trimmed_fastqs/{sample}_trimmed_{read}.fastq.gz"
	output:
		"fastqc_trimmed/{sample}_trimmed_{read}_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc_trimmed fastqc_trimmed"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		amb = "reference/{genome}.fasta.amb",
		ref = "reference/{genome}.fasta"
	output:
		"processed_bams/{sample}.{genome}.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.bwa} mem -t {params.threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - "
		"| {params.samtools} sort -O bam -o {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		metrics = "stats/{sample}.{genome}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path,
		threads = 8,
		mem = 24,
		t = long,
		tmp_dir = temp_directory
	shell:
		"{params.picard} -Xmx14g -Djava.io.tmpdir={params.tmp_dir} "
		"MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics} ASSUME_SORT_ORDER=coordinate"

rule index_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule gatk_indel_target_creator:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai",
		ref = "reference/{genome}.fasta"
	output:
		"indel_targets/{sample}.{genome}.intervals"
	params:
		gatk3 = gatk3_path,
		threads = 4,
		mem = 12,
		t = long
	shell:
		"java -jar -Xmx12g {params.gatk3} -T RealignerTargetCreator "
		"-R {input.ref} -I {input.bam} -o {output}"

rule gatk_indel_realignment:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai",
		ref = "reference/{genome}.fasta",
		indels = "indel_targets/{sample}.{genome}.intervals"
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.realigned.bam"
	params:
		gatk3 = gatk3_path,
		threads = 4,
		mem = 12,
		t = long
	shell:
		"java -jar -Xmx12g {params.gatk3} -T IndelRealigner "
		"-R {input.ref} -I {input.bam} -targetIntervals {input.indels} -o {output}"

# rule gatk_leftalignindels:
# 	input:
# 		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
# 		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai",
# 		ref = "reference/{genome}.fasta"
# 	output:
# 		"processed_bams/{sample}.{genome}.sorted.mkdup.realigned.bam"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path,
# 		threads = 4,
# 		mem = 12,
# 		t = long
# 	shell:
# 		"""{params.gatk} --java-options "-Xmx12g -Djava.io.tmpdir={params.temp_dir}" """
# 		"""LeftAlignIndels -R {input.ref} -I {input.bam} -O {output}"""

rule index_realigned_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.realigned.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.realigned.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.realigned.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.realigned.bam.bai"
	output:
		"stats/{sample}.{genome}.mkdup.sorted.realigned.bam.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule mosdepth_total:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai"
	output:
		"mosdepth_results/{sample}.{genome}.total.total.nodups_mapq30.mosdepth.summary.txt"
	params:
		mosdepth = mosdepth_path,
		prefix = "mosdepth_results/{sample}.{genome}.total.total.nodups_mapq30",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.mosdepth} --fast-mode -F 1024 --mapq 30 {params.prefix} {input.bam}"

rule gatk_gvcf_per_chrom:
	input:
		ref = "reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai"
	output:
		"gvcf/{sample}.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		chr = "{chrom}",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chr} """
		"""-ERC GVCF --do-not-run-physical-phasing -O {output}"""

rule combine_gvcfs_per_chrom:
	input:
		ref = "reference/{genome}.fasta",
		gvcfs = lambda wildcards: expand(
			"gvcf/{sample}.{genome}.{chrom}.g.vcf.gz",
			sample=config["samples"],
			genome=[wildcards.genome], chrom=[wildcards.chrom])
	output:
		v = "gvcf_combined/combined.{genome}.{chrom}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = long
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output.v}""")

rule gatk_genotypegvcf_per_chrom:
	input:
		ref = "reference/{genome}.fasta",
		gvcf = "gvcf_combined/combined.{genome}.{chrom}.g.vcf.gz"
	output:
		"vcf_genotyped/{genome}.{chrom}.gatk.called.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""GenotypeGVCFs --include-non-variant-sites -R {input.ref} -V {input.gvcf} -O {output}"""

# rule create_angsd_bam_list:
# 	input:
# 		bam = lambda wildcards: expand(
# 			"processed_bams/{sample}.{genome}.mkdup.sorted.realigned.bam",
# 			samples=config["samples"],
# 			genome=wildcards.scaffold),
# 	output:
# 		"angsd_results/{genome}_angsd_bam_list.txt"
# 	run:
# 		shell("echo -n > {output}")
# 		for i in input.bam:
# 			shell("echo {} >> {{output}}".format(i))
#
# rule angsd_by_chom_all_sites:
# 	input:
# 		bam = lambda wildcards: expand(
# 			"processed_bams/{sample}.{genome}.mkdup.sorted.realigned.bam",
# 			samples=config["samples"],
# 			genome=wildcards.scaffold),
# 		bai = expand(
# 			"processed_bams/{sample}.{genome}.mkdup.sorted.realigned.bam.bai",
# 			samples=config["samples"],
# 			genome=wildcards.scaffold),
# 		bam_list = "angsd_results/{genome}_angsd_bam_list.txt"
# 	output:
# 		geno = "angsd_results/gaga_gmor_angsd_{scaffold}.geno.gz",
# 		maf = "angsd_results/gaga_gmor_angsd_{scaffold}.mafs.gz"
# 	params:
# 		angsd = angsd_path,
# 		full_scaffold_name = lambda wildcards: config["scaffold_dict"][wildcards.scaffold],
# 		ninds = len(config["samples"]),
# 		minind = 10,
# 		mindepth = 40
# 	threads:
# 		8
# 	shell:
# 		"{params.angsd} -bam {input.bam_list} -GL 1 "
# 		"-out angsd_results/gaga_gmor_angsd_{wildcards.scaffold} -P {threads} "
# 		"-doMaf 1 -doMajorMinor 1 -doGeno 32 -doPost 1 -nind {params.nind} "
# 		"-minMapQ 20 -minQ 20 -minInd {params.minind} -skipTriallelic 1"
# 		"-setMinDepth {params.mindepth} -r '{params.full_scaffold_name}':"
#
# rule gunzip_angsd_mafs:
# 	input:
# 		"angsd_results/gaga_gmor_angsd_{scaffold}.mafs.gz"
# 	output:
# 		"angsd_results/gaga_gmor_angsd_{scaffold}.mafs"
# 	shell:
# 		"gunzip {input}"
#
# rule gunzip_angsd_geno:
# 	input:
# 		"angsd_results/gaga_gmor_angsd_{scaffold}.geno.gz"
# 	output:
# 		temp("angsd_results/gaga_gmor_angsd_{scaffold}.geno")
# 	shell:
# 		"gunzip {input}"
#
# rule concat_angsd_geno:
# 	input:
# 		expand(
# 			"angsd_results/gaga_gmor_angsd_{scaffold}.geno",
# 			scaffold=scaffolds_to_analyze)
# 	output:
# 		"angsd_results/all_angsd_CONCATENATED.geno"
# 	shell:
# 		"cat {input} > {output}"
#
# rule combine_angsd_mafs:
# 	input:
# 		expand(
# 			"angsd_results/gaga_gmor_angsd_{scaffold}.mafs",
# 			scaffold=scaffolds_to_analyze)
# 	output:
# 		"angsd_results/all_angsd_CONCATENATED.mafs"
# 	run:
# 		first = input[0]
# 		shell(
# 			"zcat {first} | head -n 1 > {output} && "
# 			"awk 'FNR-1' {input} >> {output}")
#
# rule ngs_covar_by_chrom:
# 	input:
# 		geno = "angsd_results/all_angsd_CONCATENATED.geno",
# 		mafs = "angsd_results/all_angsd_CONCATENATED.mafs"
# 	output:
# 		"angsd_results/all_angsd_CONCATENATED.covar"
# 	params:
# 		ngscovar = ngscovar_path,
# 		ninds = len(config["samples"]),
# 	shell:
# 		"{params.ngscovar} -probfile {input.geno} -outfile {output} "
# 		"-nsites $(tail -n +2 {input.mafs}) -nind {params.ninds} -call 0 -norm 0 "
# 		"-block_size 100000 -minmaf 0.01"

# rule create_interval_file_for_genomicsdbimport:
# 	input:
# 		fai = "new_reference/{genome}.fasta.fai"
# 	output:
# 		temp("interval_files/{chrom}_{genome}.bed")
# 	params:
# 		chromosome = lambda wildcards: config["scaffold_dict"][wildcards.chrom]
# 	shell:
# 		"""python scripts/Create_scaffold_bed_from_fai.py --fai {input.fai} """
# 		"""--outfile {output} --scaffold_name '{params.chromosome}'"""
#
# rule gatk_gvcf_per_chrom:
# 	input:
# 		ref = "new_reference/{genome}.fasta",
# 		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
# 		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
# 		interval = "interval_files/{chrom}_{genome}.bed"
# 	output:
# 		"gvcfs/{sample}.{genome}.{chrom}.g.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path,
# 		chromosome = lambda wildcards: config["scaffold_dict"][wildcards.chrom]
# 	threads:
# 		4
# 	shell:
# 		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
# 		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.interval} """
# 		"""-ERC GVCF -O {output}"""
#
# rule gatk_genomicsdbimport_per_chrom:
# 	input:
# 		gvcfs = lambda wildcards: expand(
# 			"gvcfs/{sample}.{genome}.{chrom}.g.vcf.gz",
# 			sample=sample_dict[wildcards.comparison],
# 			genome=[wildcards.genome],
# 			chrom=[wildcards.chrom]),
# 		ref = "new_reference/{genome}.fasta",
# 		interval = "interval_files/{chrom}_{genome}.bed"
# 	output:
# 		"gvcf_databases/{comparison}-{genome}-{chrom}"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path,
# 		chromosome = lambda wildcards: config["scaffold_dict"][wildcards.chrom]
# 	threads:
# 		8
# 	run:
# 		variant_files = []
# 		for i in input.gvcfs:
# 			variant_files.append("-V " + i)
# 		variant_files = " ".join(variant_files)
# 		shell(
# 			"""{params.gatk} --java-options "-Xmx32g -Djava.io.tmpdir={params.temp_dir}" """
# 			"""GenomicsDBImport -R {input.ref} {variant_files} """
# 			"""--genomicsdb-workspace-path {output} -L {input.interval}""")
#
# rule gatk_genotypegvcf_genomicsdb:
# 	input:
# 		gvcf = "gvcf_databases/{comparison}-{genome}-{chrom}",
# 		ref = "new_reference/{genome}.fasta"
# 	output:
# 		"vcfs/{comparison}.{genome}.{chrom}.raw.vcf.gz"
# 	params:
# 		temp_dir = temp_directory,
# 		gatk = gatk_path
# 	threads:
# 		8
# 	shell:
# 		"""{params.gatk} --java-options "-Xmx32g -Djava.io.tmpdir={params.temp_dir}" """
# 		"""GenotypeGVCFs -R {input.ref} -V gendb://{input.gvcf} -O {output}"""
#
# rule concat_vcfs:
# 	input:
# 		vcfs = lambda wildcards: expand(
# 			"vcfs/{comparison}.{genome}.{chrom}.raw.vcf.gz",
# 			comparison=wildcards.comparison,
# 			genome=wildcards.genome,
# 			chrom=config["scaffolds_no_semi_colon"])
# 	output:
# 		"final_vcfs/{comparison}.{genome}.combined.raw.vcf.gz"
# 	params:
# 		bcftools = bcftools_path
# 	shell:
# 		"{params.bcftools} concat -o {output} -Oz {input.vcfs}"

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
