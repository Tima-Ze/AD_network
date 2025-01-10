#!/usr/bin/env nextflow

// Process 2: make wTO networks
 process wTO {
    // Input files
    input:
    path file

    // Output file
    output:
    path  'results/raw_wTO/${file.basename}.txt'

    script:
    """
    Calls_wTO.R ${params.bootstrap} ${file} ${workflow.projectDir}> ${file.basename}.out

    """
}

workflow{
    main:
    wTO_net = wTO('/data/*')

    publish:
    wTO_net >> 'raw_wTO'
}