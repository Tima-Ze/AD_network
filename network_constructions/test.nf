#!/usr/bin/env nextflow

// Process 2: make wTO networks
 process wTO {
    // Input files
    input:
    path file

    // Output file
    output:
    path publishDir "${workflow.projectDir}"/Results/Raw_wTO/, mode: 'copy'

    script:
    """
    Calls_wTO.R ${params.bootstrap} ${file} ${workflow.projectDir}> ${file.baseName}.out

    """
}


// Process 2: Check the topology of the wTO networks
process RunCheckTopology1 {
    input:
    path Results/Raw_wTO
    output:
    path Results/CheckTopology/topology_wTO.txt
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
    def input_channel = Channel.fromPath('Data/*')  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process

    RunCheckTopology1()
}



