#!/usr/bin/env nextflow

nextflow.enable.dsl=2  // Enable DSL 2 syntax

process wTO {
    // Output file
    output:
    file "wto_${task.name}.out"

    // Automatically publish outputs to the directory
    publishDir "Results/Raw_wTO", mode: 'move'

    script:
    """
    # Activate Conda environment
    eval "\$(${params.conda}/bin/conda shell.bash hook)"
    conda activate R

    # Run the R script
    R CMD BATCH --vanilla "--args ${params.bootstrap} input.txt ${task.name}.txt" Calls_wTO.R wto_${task.name}.out
    """
}

workflow {
    wTO()
}
