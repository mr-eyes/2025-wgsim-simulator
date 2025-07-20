import os
import pandas as pd

# Load configuration
configfile: "config.yaml"

# Load sample information
samples_df = pd.read_csv(config["samples_file"], sep="\t", comment="#")
SAMPLES = samples_df["sample"].tolist()


rule all:
    input:
        # Raw simulated reads
        expand("results/{sample}/{sample}_R1.fq", sample=SAMPLES),
        expand("results/{sample}/{sample}_R2.fq", sample=SAMPLES),
        # Quality control
        expand("results/{sample}/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("results/{sample}/fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        # Summary reports
        "reports/multiqc_report.html",
        "reports/simulation_summary.txt"

rule simulate_reads:
    """
    Simulate paired-end reads using wgsim
    """
    input:
        ref = config["reference_genome"]
    output:
        r1 = "results/{sample}/{sample}_R1.fq",
        r2 = "results/{sample}/{sample}_R2.fq"
    params:
        sample_config = lambda wildcards: samples_df[samples_df["sample"] == wildcards.sample].iloc[0]
    log:
        "logs/{sample}_wgsim.log"
    threads: 16
    resources:
        mem_mb = 80000,
        time_min = 480  # 8 hours
    shell:
        """
        echo "Starting wgsim simulation for {wildcards.sample}" > {log}
        echo "Parameters:" >> {log}
        echo "  Read pairs: {params.sample_config[n_reads]}" >> {log}
        echo "  Read length: {params.sample_config[read_length]}" >> {log}
        echo "  Error rate: {params.sample_config[error_rate]}" >> {log}
        echo "  Mutation rate: {params.sample_config[mutation_rate]}" >> {log}
        echo "  Seed: {params.sample_config[seed]}" >> {log}
        echo "" >> {log}
        
        # Unzip reference if needed
        if [[ {input.ref} == *.gz ]]; then
            gunzip -c {input.ref} > results/{wildcards.sample}/ref_temp.fasta
            REF_FILE="results/{wildcards.sample}/ref_temp.fasta"
        else
            REF_FILE="{input.ref}"
        fi
        
        # Create output directory
        mkdir -p results/{wildcards.sample}
        
        # Run wgsim
        wgsim \
            -N {params.sample_config[n_reads]} \
            -1 {params.sample_config[read_length]} \
            -2 {params.sample_config[read_length]} \
            -d {params.sample_config[fragment_mean]} \
            -s {params.sample_config[fragment_std]} \
            -e {params.sample_config[error_rate]} \
            -r {params.sample_config[mutation_rate]} \
            -R {params.sample_config[indel_fraction]} \
            -X {params.sample_config[indel_extended]} \
            -S {params.sample_config[seed]} \
            $REF_FILE \
            {output.r1} \
            {output.r2} \
            2>> {log}
        
        # Clean up temporary reference
        if [[ {input.ref} == *.gz ]]; then
            rm -f results/{wildcards.sample}/ref_temp.fasta
        fi
        
        echo "Simulation completed successfully" >> {log}
        """

rule fastqc:
    """
    Run FastQC on simulated reads
    """
    input:
        "results/{sample}/{sample}_{read}.fq"
    output:
        html = "results/{sample}/fastqc/{sample}_{read}_fastqc.html",
        zip = "results/{sample}/fastqc/{sample}_{read}_fastqc.zip"
    params:
        outdir = "results/{sample}/fastqc"
    log:
        "logs/{sample}_{read}_fastqc.log"
    threads: 16
    resources:
        mem_mb = 20000,
        time_min = 60
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {threads} {input} > {log} 2>&1
        """

rule multiqc:
    """
    Aggregate FastQC reports with MultiQC
    """
    input:
        expand("results/{sample}/fastqc/{sample}_{read}_fastqc.zip", 
               sample=SAMPLES, read=["R1", "R2"])
    output:
        "reports/multiqc_report.html"
    params:
        outdir = "reports"
    log:
        "logs/multiqc.log"
    threads: 16
    resources:
        mem_mb = 40000,
        time_min = 30
    shell:
        """
        multiqc results/*/fastqc/ -o {params.outdir} -f > {log} 2>&1
        """

rule simulation_summary:
    """
    Generate summary of simulation results
    """
    input:
        expand("results/{sample}/{sample}_R1.fq", sample=SAMPLES),
        expand("results/{sample}/{sample}_R2.fq", sample=SAMPLES)
    output:
        "reports/simulation_summary.txt"
    log:
        "logs/simulation_summary.log"
    shell:
        """
        echo "wgsim Dog Genome Simulation Summary" > {output}
        echo "==================================" >> {output}
        echo "Generated on: $(date)" >> {output}
        echo "" >> {output}
        
        for sample in {SAMPLES}; do
            echo "Sample: $sample" >> {output}
            echo "  Directory: results/$sample/" >> {output}
            
            if [[ -f results/$sample/${{sample}}_R1.fq ]]; then
                r1_size=$(du -h results/$sample/${{sample}}_R1.fq | cut -f1)
                r1_reads=$(cat results/$sample/${{sample}}_R1.fq | wc -l | awk '{{print $1/4}}')
                echo "  R1 file: ${{sample}}_R1.fq ($r1_size, $r1_reads reads)" >> {output}
            fi
            
            if [[ -f results/$sample/${{sample}}_R2.fq ]]; then
                r2_size=$(du -h results/$sample/${{sample}}_R2.fq | cut -f1)
                r2_reads=$(cat results/$sample/${{sample}}_R2.fq | wc -l | awk '{{print $1/4}}')
                echo "  R2 file: ${{sample}}_R2.fq ($r2_size, $r2_reads reads)" >> {output}
                
                # Calculate estimated coverage
                total_bases=$((r1_reads * 2 * 150))  # Assuming 150bp reads
                coverage=$(echo "scale=1; $total_bases / 2410976875" | bc -l)
                echo "  Estimated coverage: ${{coverage}}X" >> {output}
            fi
            
            echo "" >> {output}
        done
        
        echo "Summary complete" > {log}
        """
