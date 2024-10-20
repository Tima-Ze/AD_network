#!/usr/bin/env nextflow

nextflow.enable.dsl=2  // Enable DSL 2 syntax

process wTO {
    // Output file
    output:
    file "wto_${task.name}.txt"

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
// Process 2: Check the topology of the wTO networks
process RunCheckTopology1 {
    input:
    path Results/Raw_wTO  // Input file to be processed
    output:
    path Results/CheckTopology/topology_output.txt  // Output of the first R script
    // Automatically publish outputs to the directory
    publishDir "Results/CheckTopology/", mode: 'move'

    script:
    """
    # Activate Conda environment
    eval "\$(${params.conda}/bin/conda shell.bash hook)"
    conda activate R

    # Run the first R script with input and output directories passed as arguments
    R CMD BATCH --vanilla "--args $input_file $output_file" check_topology_wTO.R
    """
}
workflow {
    wTO()
}
