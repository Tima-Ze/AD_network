#!/usr/bin/env nextflow

// Process 2: make wTO networks
 process wTO {
    // Input files
    input:
    path file

    // Output file
    output:
    path publishDir "${workflow.projectDir}"/results/raw_wTO/, mode: 'copy'

    script:
    """
    Calls_wTO.R ${params.bootstrap} ${file} ${workflow.projectDir}> ${file.baseName}.out

    """
}

workflow {
    def input_channel = Channel.fromPath('data/*')  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process
}
