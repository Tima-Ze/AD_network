#!/usr/bin/env nextflow

// Define parameters
params.output_dir = 'results/raw_wTO/'

process wTO {
    input:
    path file

    output:

     path "${file.baseName}.out"

    script:
    """
    Calls_wTO.R ${params.bootstrap} data/${file} ${params.output_dir} ${workflow.projectDir} > ${file.baseName}.out
    """
}

workflow {
    def input_channel = Channel.fromPath('data/*')  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process
}
