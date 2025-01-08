#!/usr/bin/env nextflow

nextflow.enable.dsl=2  // Enable DSL 2 syntax

 process wTO {
    // Input files
    input:
    path file from Channel.fromPath("Data/*")  // This automatically handles all files in Data directory

    // Output file
    file "${file.getBaseName()}.txt"

    // Automatically publish outputs to the directory
    publishDir "Results/Raw_wTO", mode: 'move'

    script:
    """
    # Activate Conda environment
    R CMD BATCH --vanilla "--args ${params.bootstrap} $file ${file.getBaseName()}.txt" Calls_wTO.R wto_${file.getBaseName()}.out
    """
}


// Process 2: Check the topology of the wTO networks
process RunCheckTopology1 {
    input:
    path Results/Raw_wTO  // Input file to be processed
    output:
    path Results/CheckTopology/topology_wTO.txt  // Output of the first R script
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

process filter_wTO
{
    input: //files from outpot process wTO
    path wTO_output from wTO.out  // Output of the first process as input
    output:
    path Results/Filtered_wTO
    publishDir "Results/Filtered_wTO", mode: 'move'
    script:
    """
    # Activate Conda environment
    eval "\$(${params.conda}/bin/conda shell.bash hook)"
    conda activate R

    # Run the first R script with input and output directories passed as arguments
    R CMD BATCH --vanilla "--args $input_file $output_file" filter_wTO.R
    """
}
workflow {
    input_channel = Channel.fromPath("Data/*")  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process

    RunCheckTopology1()
}
