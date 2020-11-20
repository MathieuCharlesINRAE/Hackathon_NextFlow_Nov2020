#!/bin/env nextflow

//Not working
// how to concatenate a parameter value and a string ?

params.matrixPreprocessedFolder = "$baseDir/results/matrix_processed"
matrixPreprocessedFolder_ch = Channel.fromPath("${params.matrixPreprocessedFolder}/*.pfm")
                                     .map { file -> [ file.name.replace(".pfm",""), file,  "${params.matrixPreprocessedFolder}/" + file.name.replace(".pfm",".score_distrib"), "${params.matrixPreprocessedFolder}/" + file.name.replace(".pfm",".ratio_distrib") ] }
                                     //.view()


//params.matrixFolder = "$baseDir/matrix_firstAnalysis"
params.matrixFolder = "$baseDir/matrix_secondAnalysis"

matrix_ch = Channel.fromPath("${params.matrixFolder}/*.pwm")
                   .map { file -> [ file.name.replace(".pwm",""), file ] }
                 //.view()

process check_if_preprocessed{
    
    input:
    tuple val(matrix_name), file( pwm) from matrix_ch
    
    output:
    tuple val(matrix_name), file( "*.pwm")into matrices_toProcessed_ch
    
    script:
    """
    if [ ! -e ${params.matrixPreprocessedFolder}"/"$matrix_name".pfm" ]
    then
        echo $matrix_name, $pwm
    fi
    """
}

matrices_toProcessed_ch.view()

/*
process prepare_matrix {
    
    publishDir "$baseDir/results/matrix_processed", mode: 'link'
    
    input :
    tuple val(matrix_name), file(matrix_file) from matrix_ch
    
    output:
    tuple val(matrix_name), file( "*.pfm"), file("*.score_distrib"), file("*.ratio_distrib") into matrices_processed_ch

    
    script:
    """
    touch ${matrix_name}.pfm
    touch ${matrix_name}.score_distrib 
    touch ${matrix_name}.ratio_distrib
    """
}

//matrices_processed_ch.view()

// NOT WORKING, No such variant x

variant_list = ["variant1","var2","var3"]

variant_ch = Channel.fromList(variant_list)

compare_ch = variant_ch.combine(matrices_processed_ch)

process do_something {
    
    input:
    tuple val(x), val(matrix_name), file(pfm), file(score), file(ratio) from compare_ch
    
    output:
    env impact into result_ch
    
    script:
    """
    impact=\$( echo compute_impact.pl --variant $x --matrice $pfm --score $score --ratio $ratio)
    """
    
}

result_ch.view()

// Now I want to use and other matrixFolder


// but I do not want that matrices M00001 and M00002 are processed again
// add an other channel for previously process matrices ?
// use mix and collect
*/