#!/home/datnguyen/nextflow nextflow

// Process for Quality Control of raw reads
process FASTQC {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    """
    fastqc -f fastq --threads ${task.cpus} -q ${reads}
    """
}

process TRIMMOMATIC {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_paired*"), emit: paired
    tuple val(sample_id), path("*_unpaired*"), emit: unpaired
    tuple val(sample_id), path("*.summary.txt"), emit: summary
    tuple val(sample_id), path("*.log"), emit:log


    script:
    fq_1_paired = sample_id + "*_paired.1.fq.gz"
    fq_1_unpaired = sample_id + "_unpaired.1.fq.gz"
    fq_2_paired = sample_id + "_paired.2.fq.gz"
    fq_2_unpaired = sample_id + "_unpaired.2.fq.gz"

    """
    trimmomatic PE \
    -threads ${task.cpus} \
    -trimlog ${sample_id}.log \
    -summary ${sample_id}.summary.txt \
    ${reads[0]} ${reads[1]} \
    ${fq_1_paired} ${fq_1_unpaired} \
    ${fq_2_paired} ${fq_2_unpaired} \
    SLIDINGWINDOW:4:25 \
    MINLEN:100
    """
}


process FASTQC_2 {
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("*")

    script:
    """
    fastqc -f fastq -q ${reads} -t ${task.cpus}
    """
}

process KRAKEN_READ {
    input:
    tuple val(sample_id), path(reads)
    // path(k2_database)

    output:
    tuple val(sample_id), path("${sample_id}_classification.txt"), emit: classification
    tuple val(sample_id), path("${sample_id}_report.txt"), emit: report
    
    script:
    def k2_classification = "${sample_id}_classification.txt"
    def k2_report = "${sample_id}_report.txt"

    """
    kraken2 --paired \
    --db ${params.k2_database} \
    --threads ${task.cpus} \
    --output ${k2_classification} \
    --report ${k2_report} \
    ${reads[0]} ${reads[1]}
    """
}

process SPADES {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_contigs.fasta"), emit: contigs
    tuple val(sample_id), path("*_scaffolds.fasta"), emit: scaffolds
    tuple val(sample_id), path("*_graph.fastg"), emit: graph

    script:
    def contigs = "${sample_id}_contigs.fasta"
    def scaffolds = "${sample_id}_scaffolds.fasta"
    def graph = "${sample_id}_graph.fastg"

    """
    spades.py --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} --threads ${task.cpus} -o ./
    mv contigs.fasta ${contigs}
    mv scaffolds.fasta ${scaffolds}
    mv assembly_graph.fastg ${graph}
    """
}

process BANDAGE {
    input:
    tuple val(sample_id), path(graph)

    output:
    tuple val(sample_id), path("*_info.tsv"), emit: info
    tuple val(sample_id), path("*_image.png"), emit: image

    script:
    def info = "${sample_id}_info.tsv"
    def image = "${sample_id}_image.png"

    """
    Bandage info ${graph} --tsv > ${info}
    Bandage image ${graph} ${image}
    """
}

process SEQKIT {
    input:
    tuple val(sample_id), path(contigs)
    tuple val(sample_id), path(scaffolds)

    output:
    tuple val(sample_id), path("*_contigs.filtered.fasta"), emit: filtered_contigs
    tuple val(sample_id), path("*_scaffolds.filtered.fasta"), emit: filtered_scaffolds


    script:
    def filtered_contigs = "${sample_id}_contigs.filtered.fasta"
    def filtered_scaffolds = "${sample_id}_scaffolds.filtered.fasta"
    """
    seqkit seq -m 500 ${contigs} > ${filtered_contigs}
    seqkit seq -m 500 ${scaffolds} > ${filtered_scaffolds}
    """
}

process QUAST {
    input:
    tuple val(sample_id), path(filtered_contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    def contigs_qc = "${sample_id}_contigs"
    def scaffolds_qc = "${sample_id}_scaffolds"

    """
    quast.py -o ${contigs_qc} --threads ${task.cpus} ${filtered_contigs}
    quast.py -o ${scaffolds_qc} --threads ${task.cpus} ${filtered_scaffolds}
    """
}

process BUSCO {
    input:
        tuple val(sample_id), path(filtered_contigs)

    output:
        tuple val(sample_id), path("${sample_id}/*.txt")

    script:
    """
    busco -i ${filtered_contigs} \
    -m genome -c ${task.cpus} \
    -o ${sample_id} \
    --out_path ./ \
    --auto-lineage-prok \
    -f
    """
}

process PROKKA {
    input:
        tuple val(sample_id), path(filtered_contigs)

    output:
        tuple val(sample_id), path("*.faa"), emit: FAA
        tuple val(sample_id), path("*")

    script:
    """
    prokka --kingdom Bacteria \
    --prefix ${sample_id}_annotation \
    --cpus ${task.cpus} \
    --force \
    --outdir ./ \
    ${filtered_contigs}
    """
}

process ABRICATE {
    input:
    tuple val(sample_id), path(filtered_contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    abricate --threads ${task.cpus} --db argannot --csv ${filtered_contigs} > ${sample_id}_argannot.AMR.csv
    abricate --threads ${task.cpus} --db card --csv ${filtered_contigs} > ${sample_id}_card.AMR.csv
    abricate --threads ${task.cpus} --db megares --csv ${filtered_contigs} > ${sample_id}_megares.AMR.csv
    abricate --threads ${task.cpus} --db ncbi --csv ${filtered_contigs} > ${sample_id}_ncbi.AMR.csv
    abricate --threads ${task.cpus} --db resfinder --csv ${filtered_contigs} > ${sample_id}_resfinder.AMR.csv
    abricate --threads ${task.cpus} --db plasmidfinder --csv ${filtered_contigs} > ${sample_id}_plasmidfinder.csv
    abricate --threads ${task.cpus} --db vfdb --csv ${filtered_contigs} > ${sample_id}_vfdb.csv
    """
}

process MLST {
    input:
        tuple val(sample_id), path(filtered_contigs)

    output:
        tuple val(sample_id), path("*")

    script:
    """
    mlst ${filtered_contigs} > ${sample_id}_mlst.csv
    """
}

process EGGNOGMAPPER {
    input:
    tuple val(sample_id), path(faa)

    output:
    tuple val(sample_id), path("*.annotations"), emit:annotations
    tuple val(sample_id), path("*")

    script:
    """
    emapper.py -m diamond --itype proteins -i ${faa} -o ${sample_id}.emapper --cpu ${task.cpus} 
    """
}

process COGCLASSIFIER {
    input:
    tuple val(sample_id), path(eggnog_table)

    output:
    tuple val(sample_id), path("*_table.tsv"), emit: table
    tuple val(sample_id), path("*.png"), emit: image

    script:
    def cog_table = "${sample_id}_table.tsv"
    def cog_barchart = "${sample_id}_barchart.png"
    def cog_piechart = "${sample_id}_piechart.png"
    """
    python3 ./COG_analysis.py -i ${eggnog_table} -o ${cog_table}
    python3 ./COG_chart.py -i ${cog_table} -o ${cog_barchart}
    python3 ./COG_piechart.py -i ${cog_table} -o ${cog_piechart}
    """
}

process GOCLASSIFIER {
    input:
    tuple val(sample_id), path(eggnog_table)

    output:
    tuple val(sample_id), path("*_table.GO.tsv"), emit:table
    tuple val(sample_id), path("*.png"), emit: image

    script:
    def go_table = "${sample_id}_table.GO.tsv"
    def go_chart = "${sample_id}_barchart.GO.png"
    """
    python3 ./GO_analysis.py -i ${eggnog_table} -o ${go_table} -db ${params.go_database}
    python3 ./GO_chart.py -i ${go_table} -o ${go_chart}
    """
}

process BWA_INDEX {
    input:
    path(reference)

    output:
    // path("*.fna"), emit: reference
    path("*.amb"), emit: amb
    path("*.ann"), emit: ann
    path("*.0123"), emit: test
    path("*.bwt.*"), emit: bwt
    path("*.pac"), emit: pac

    script:
    """
    bwa-mem2 index -p ${reference.getBaseName()} ${reference}
    """
}

process BWA_MEM {
    input:
    tuple val(sample_id), path(reads)
    path(index)

    output:
    tuple val(sample_id), path("*.sam"), emit: alignment

    script:
    def output = "${sample_id}_alignment.sam"
    """
    bwa-mem2 mem -o ${output} -t ${task.cpus} ${index} ${reads[0]} ${reads[1]} 
    """
}

process SAMTOOLS {
    input:
    tuple val(sample_id), path(alignment)

    output:
    tuple val(sample_id), path("*.sorted.bam"), emit: sorted_alignment

    script:
    def sorted_output = "${sample_id}_alignment.sorted.bam"
    """
    samtools sort ${alignment} -o ${sorted_output} -@ ${task.cpus}
    """
}

process FREEBAYES {
    input:
    tuple val(sample_id), path(sorted_alignment)

    output:
    tuple val(sample_id), path("*.vcf") ,emit: vcf

    script:
    def output = "${sample_id}.vcf"
    """
    freebayes -b ${sorted_alignment} -f ${params.reference} -v ${output} -p 1 
    """
}

process SNPEFF_BUILD {
    output:
    path("data")
    path("snpEff.config"), emit: config

    script:
    """
    touch snpEff.config
    cat > snpEff.config <<EOF
    #-----
    # Core
    #-----

    #-------------------------------
    # SnpEff configuration file
    #-------------------------------

    data.dir = ./data/

    #-------------------------------
    # Codon Tables
    #-------------------------------

    codon.Standard								: TTT/F , TTC/F , TTA/L , TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Vertebrate_Mitochondrial				: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/M+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/* , AGG/* , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Yeast_Mitochondrial					: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/T , CTC/T , CTA/T , CTG/T , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/M+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Mold_Mitochondrial					: TTT/F , TTC/F , TTA/L+, TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Protozoan_Mitochondrial				: TTT/F , TTC/F , TTA/L+, TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Coelenterate							: TTT/F , TTC/F , TTA/L+, TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Mitochondrial							: TTT/F , TTC/F , TTA/L+, TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Mycoplasma							: TTT/F , TTC/F , TTA/L+, TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Spiroplasma							: TTT/F , TTC/F , TTA/L+, TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Invertebrate_Mitochondrial			: TTT/F , TTC/F , TTA/L , TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/M+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/S , AGG/S , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Ciliate_Nuclear						: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/Q , TAG/Q , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Dasycladacean_Nuclear					: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/Q , TAG/Q , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Hexamita_Nuclear						: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/Q , TAG/Q , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Echinoderm_Mitochondrial				: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/N , AAG/K , AGT/S , AGC/S , AGA/S , AGG/S , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Flatworm_Mitochondrial				: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/N , AAG/K , AGT/S , AGC/S , AGA/S , AGG/S , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Euplotid_Nuclear						: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/C , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Bacterial_and_Plant_Plastid			: TTT/F , TTC/F , TTA/L , TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I+, ATA/I+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Alternative_Yeast_Nuclear				: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/S+, CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Ascidian_Mitochondrial				: TTT/F , TTC/F , TTA/L , TTG/L+, TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/M+, ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/G , AGG/G , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Alternative_Flatworm_Mitochondrial	: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/Y , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/N , AAG/K , AGT/S , AGC/S , AGA/S , AGG/S , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Blepharisma_Macronuclear				: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/Q , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Chlorophycean_Mitochondrial			: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/L , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Trematode_Mitochondrial				: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/W , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/M , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/N , AAG/K , AGT/S , AGC/S , AGA/S , AGG/S , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Scenedesmus_obliquus_Mitochondrial	: TTT/F , TTC/F , TTA/L , TTG/L , TCT/S , TCC/S , TCA/* , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/L , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I , ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V , GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G
    codon.Thraustochytrium_Mitochondrial		: TTT/F , TTC/F , TTA/* , TTG/L , TCT/S , TCC/S , TCA/S , TCG/S , TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/* , TGG/W , CTT/L , CTC/L , CTA/L , CTG/L , CCT/P , CCC/P , CCA/P , CCG/P , CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , ATT/I+, ATC/I , ATA/I , ATG/M+, ACT/T , ACC/T , ACA/T , ACG/T , AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , GTT/V , GTC/V , GTA/V , GTG/V+, GCT/A , GCC/A , GCA/A , GCG/A , GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G

    # Databases
    reference.genome : Bacillus subtilis str168 v250725
    reference.M.codonTable : Standard

    EOF

    mkdir -p data/reference
    wget ${params.ref_path} -O data/reference/sequences.fa.gz
    wget ${params.annot_path} -O data/reference/genes.gtf.gz
    snpEff build -gtf22 -noCheckCds -noCheckProtein -c snpEff.config -v reference

    """
}

process SNPEFF_ANN {
    input:
    path(config)
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("*.csv"), emit: stats_csv
    tuple val(sample_id), path("*.html"), emit: stats_html
    tuple val(sample_id), path("*.ann.vcf"), emit: ann_vcf

    script:
    def stats_csv = "${sample_id}.stats.csv"
    def stats_html = "${sample_id}.stats.html"
    def ann_vcf = "${sample_id}.ann.vcf"
    """
    snpEff -csvStats ${stats_csv} -i vcf -o vcf -htmlStats ${stats_html} -c ${config} -v reference ${vcf} > ${ann_vcf}
    """
}

workflow {
    Channel.fromFilePairs(params.reads, checkIfExists: true).set { raw_reads_channel }

    fastqc_channel = FASTQC(raw_reads_channel)
    trim_channel = TRIMMOMATIC(raw_reads_channel)

    fastqc_channel_2 = FASTQC_2(TRIMMOMATIC.out.paired)
            
    kraken_channel = KRAKEN_READ(TRIMMOMATIC.out.paired)
    spades_channel = SPADES(TRIMMOMATIC.out.paired)

    bandage_channel = BANDAGE(SPADES.out.graph)
    seqkit_channel = SEQKIT(SPADES.out.contigs, SPADES.out.scaffolds)

    quast_channel = QUAST(seqkit_channel)
    busco_channel = BUSCO(SEQKIT.out.filtered_contigs)

    prokka_channel = PROKKA(SEQKIT.out.filtered_contigs)

    abricate_channel = ABRICATE(SEQKIT.out.filtered_contigs)
    mlst_channel = MLST(SEQKIT.out.filtered_contigs)
    eggnog_channel = EGGNOGMAPPER(PROKKA.out.FAA)

    cog_channel = COGCLASSIFIER(EGGNOGMAPPER.out.annotations)
    go_channel = GOCLASSIFIER(EGGNOGMAPPER.out.annotations)
    // go_classifier_channel = GOCLASSIFIER(go_input_channel)

    Channel.value(params.reference).set { ref_channel }

    index_channel = BWA_INDEX(ref_channel)
    alignment_channel = BWA_MEM(TRIMMOMATIC.out.paired, BWA_INDEX.out.ann)
    samtools_channel = SAMTOOLS(alignment_channel)
    freebayes_channel = FREEBAYES(samtools_channel)
    build_channel = SNPEFF_BUILD()
    annot_channel = SNPEFF_ANN(SNPEFF_BUILD.out.config, freebayes_channel)
}



/*
==========================================================================
Workflow Event Handler
*/

workflow.onComplete {

    println(workflow.success ? """
    Pipeline execution summary
    ------------------------------
    Completed at: ${workflow.complete}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    workDir: ${workflow.workDir}
    exit Status: ${workflow.exitStatus}
    """:"""
    Failed: ${workflow.errorReport}
    exit status: ${workflow.exitStatus}
    
    """
    )
}
