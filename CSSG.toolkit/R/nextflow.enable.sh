nextflow.enable.dsl=2

params.fasta_dir = pwd()
params.adapters_path = file("${which("IndexPip")}/requirements_file/Adapters.fa")
params.results_dir = file("results")
params.tmp_dir = file("tmp")

process PrepareDirectories {
    output:
    path "results/matrices", emit: results_matrices
    path "tmp", emit: tmp_dir

    script:
    """
    mkdir -p results
    mkdir -p results/matrices
    mkdir -p tmp
    """
}

process ProcessFASTQ {
    tag "$sample"
    input:
    path fastq1
    path fastq2
    path adapters
    val umi_length
    val reads_length
    val species
    val qc_reads
    val annotation_names
    val annotation_side
    val multi
    val CPU
    path genome_dir
    path tmp_dir
    path results_dir

    output:
    path "${sample}_processed.bam", emit: processed_bam

    script:
    """
    fastp -i $fastq1 -I $fastq2 -o $tmp_dir/trimmed_$(basename $fastq1) \
          -O $tmp_dir/trimmed_$(basename $fastq2) --adapter_fasta $adapters \
          --trim_poly_x --length_required $umi_length --html $results_dir/${sample}_quality_report.html

    umi_tools extract -I $tmp_dir/trimmed_$(basename $fastq2) --bc-pattern="$(printf 'N%.0s' $(seq 1 $umi_length))" \
                 --read2-in=$tmp_dir/trimmed_$(basename $fastq1) --stdout=$tmp_dir/extracted_$(basename $fastq2) \
                 --read2-out=$tmp_dir/umi_$(basename $fastq1)

    STAR --outReadsUnmapped Fastx --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 3 --runThreadN $CPU --genomeDir $genome_dir/index/$reads_length \
         --readFilesIn ${tmp_dir}/umi_$(basename $fastq1) --outFileNamePrefix $tmp_dir/${sample}_

    samtools index $tmp_dir/${sample}_Aligned.sortedByCoord.out.bam
    umi_tools dedup -I $tmp_dir/${sample}_Aligned.sortedByCoord.out.bam -S $results_dir/${sample}_deduplicated.bam
    """
}

process GenerateMatrix {
    tag "$sample"
    input:
    path processed_bam
    path genome_dir
    val annotation_names
    val annotation_side
    val multi
    val CPU
    path results_dir

    output:
    path "results/matrices/${sample}_count_matrix.txt", emit: count_matrix

    script:
    """
    for name in ${annotation_names.split(',')}; do
        for side in ${annotation_side.split(',')}; do
            if [[ $multi == "TRUE" ]]; then
                featureCounts -a $genome_dir/correct_annotation.gtf -o $results_dir/matrices/${name}_${side}_${sample}_genes_count_matrix.txt \
                             -f -M -t $side -g $name -T $CPU $processed_bam
            else
                featureCounts -a $genome_dir/correct_annotation.gtf -o $results_dir/matrices/${name}_${side}_${sample}_genes_count_matrix.txt \
                             -f -t $side -g $name -T $CPU $processed_bam
            fi
        done
    done
    """
}

workflow {
    prepareDirs = PrepareDirectories()
    genome_dir = file("${which("IndexPip")}/genome/${params.species}")
    fastq_files = Channel.fromPath("*.fastq.gz").groupBy { it.baseName.replaceAll(/_R[12]/, '') }

    processed_bams = fastq_files.flatMap { sample, files ->
        def fastq1 = files.find { it.name.contains("_R1") }
        def fastq2 = files.find { it.name.contains("_R2") }
        ProcessFASTQ(fastq1, fastq2, params.adapters_path, params.umi_length, params.reads_length, 
                     params.species, params.qc_reads, params.annotation_names, 
                     params.annotation_side, params.multi, params.CPU, genome_dir, 
                     prepareDirs.tmp_dir, prepareDirs.results_matrices)
    }

    processed_bams.view()
    processed_bams | GenerateMatrix(genome_dir, params.annotation_names, params.annotation_side, params.multi, params.CPU, prepareDirs.results_matrices)
}
