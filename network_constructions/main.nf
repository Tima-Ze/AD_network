#!/usr/bin/env nextflow

// Define parameters
params.output_raw = 'results/raw_wTO/'
params.output_filter = 'results/filter_wTO/'

//calculate wTO networks
process wTO {
    input:
    path file

    output:

     path "${file}" into wto_raw_channel

    script:
    """
    Calls_wTO.R ${params.bootstrap} data/${file} ${params.output_raw} ${workflow.projectDir} > ${file.baseName}.out
    """
}

//subset links with absulute value of wTO >= 0.5
process filter_wTO
{
    input: //files from outpot process wTO
    path file from wto_raw_channel

    output:
    path ${file} into wto_filter_channel

    script:
    """
    filter_wTO.R ${params.output_raw.file} ${params.output_filter} ${workflow.projectDir}
    """
}

workflow {
    def input_channel = Channel.fromPath('data/*')  // Define the input channel
    wTO(input_channel)  // Pass the input channel to the process
    filter_wTO(wto_raw_channel)
}