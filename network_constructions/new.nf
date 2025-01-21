#!/usr/bin/env nextflow

// Define parameters
params.output_raw = 'results/raw_wTO/'
params.output_filter = 'results/filter_wTO/'
params.output_filter = 'results/cns/'

//calculate wTO networks
process wTO {
    input:
    path file

    output:

     path "${file.baseName}.out"

    script:
    """
    Calls_wTO.R ${params.bootstrap} data/${file} ${params.output_raw} ${workflow.projectDir} > ${file.baseName}.out
    """
}

//subset links with absulute value of wTO >= 0.5
process filter_wTO {
    input: //files from outpot process wTO
    path file

    output:
    path "${file}"

    script:
    """
    filter_wto.R ${params.output_raw}${file} ${params.output_filter} ${workflow.projectDir}
    """
}

workflow {
    def input_channel = Channel.fromPath('data/*')  // Define the input channel
    def wto_raw_channel = Channel.fromPath("${params.output_raw}*")
    wTO(input_channel)  // Pass the input channel to the process
    filter_wTO(wto_raw_channel)
}