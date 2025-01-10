#!/usr/bin/env nextflow

// Process 2: make wTO networks
 process wTO {
    // Input files
    input:
    path file

    // Output file
    output:
    path  'results/raw_wTO/${file.baseName}.txt'

    script:
    """
    Calls_wTO.R ${params.bootstrap} ${file} ${workflow.projectDir}> ${file.baseName}.out

    """
}

workflow {
    main:
    wTO_net = wTO("${workflow.projectDir}"/data/)

    publish:
    wTO_net >> 'raw_wTO'
}