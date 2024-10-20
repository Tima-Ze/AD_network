#!/usr/bin/env nextflow

// nextflow.enable.dsl=2

//This process will execute Calls.wTO R script on SLURM, the script then internally calls two more R script
process Calls_wTO {
    // Output: The file generated by the R script
    output:
    path "Results/Raw_wTO/wto_${task.name}.out"

    script:
    """
    // Create the output directory if it does not exist
    mkdir -p Results/Raw_wTO

    # Activate Conda environment
    eval "\$(${params.conda}/bin/conda shell.bash hook)"
    conda activate R

    # Run the R script
    R CMD BATCH --vanilla "--args ${params.Bootstrap} input.txt ${task.name}.txt" Calls_wTO.R wto_${task.name}.out
    """
}

workflow {
    // Workflow definition, this will trigger the process
    Calls_wTO()
}
