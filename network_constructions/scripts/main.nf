#!/usr/bin/env nextflow

 process wTO {
    // Input files
    input:
    path file

    // Output file
    ???

    // Automatically publish outputs to the directory
    publishDir "Results/Raw_wTO", mode: 'move'

    script:
    """
    Rscript ${workflow.projectDir}/scripts/Calls_wTO.R ${params.bootstrap} ${file}> ${file}_wto.out

    """
}

workflow {
    def input_channel = Channel.fromPath('Data/*')  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process
}
