#!/usr/bin/env nextflow

process wTO {
    // Input files
    input:
    path file

    // Output files
    output:
    path "${file.baseName}.txt"  // Specify the output file to be captured
    
    script:
    """
    Calls_wTO.R ${params.bootstrap} ${file} ${workflow.projectDir} > ${file.baseName}.out
    """
}

workflow {
    def input_channel = Channel.fromPath('Data/*')  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process
}
