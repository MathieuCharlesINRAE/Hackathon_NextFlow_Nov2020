#!/bin/env nextflow

//Not working
// how to concatenate a parameter value and a string ?

/* 
params.matrixFolder = "$baseDir/matrix_firstAnalysis"

matrix_ch = Channel.fromPath(params.matrixFolder"/*.pfm")
                   .map { file -> [ file.name.replace(".pfm",""), file ] }
                  .view()
*/ 

params.matrixFolder = "$baseDir/matrix_firstAnalysis/*.pfm"

matrix_ch = Channel.fromPath(params.matrixFolder)
                   .map { file -> [ file.name.replace(".pfm",""), file ] }
                  .view()

process prepare_matrix {
    
    publishDir "results/matrix_processed"
    
    input :
    tuple val(matrix_name), file(matrix_file) from matrix_ch
    
    output:
    file "*" into matrices_processed_ch
    
    script:
    """
    touch $matrix_name".pfm" 
    # touch $matrix_name".score_distrib" 
    # touch $matrix_name".ratio_distrib" 
    """
}

// personnal note : to not run again exact same process, we must use the -resume option.
// How to group the 3 files in one element of the channel matrices_processed_ch

// NOT WORKING, No such variant x

variant_ch = Channel.from("variant1")
                    .view()

process do_something {

    echo true
    
    input:
    vat(x) from variant_ch
    file(pfm) from matrices_processed_ch
    
    output:
    val(score) into stdout
    
    script:
    """
    echo compute_impact.pl --variant $x --matrice $pfm
    """
    
}

// Now I want to use and other matrixFolder
params.matrixFolder = "$baseDir/matrix_secondAnalysis"

// but I do not want that matrices M00001 and M00002 are processed again