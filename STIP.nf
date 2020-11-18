#!/usr/bin/env nextflow

params.genome  = 'genome/genome.chr1.1M.fa'
genome_ch = Channel.fromPath(params.genome)

params.gtf  = 'genome/genome.chr1.1M.gtf'
gtf_ch = Channel.fromPath(params.gtf)

params.snp  = 'snp/bos_taurus.SNP.biallelic.100k.bed'
snp_ch = Channel.fromPath(params.snp)


process generateWindowsFromGtf {
    input:
    file genome from gtf_ch
    
    output:
    file 'windows.bed' into generateWindowsFromGtfOut

    shell:
    '''
    awk -F "\t" 'BEGIN{OFS="\t"}{ if ($3 == "gene"){ if ($7 == "+"){ print $1,$4-2000,$4}else{print $1,$5,$5+2000}}}' !{genome} > windows.bed
    '''
}

process intersectVariantWithWindows {
    conda "bioconda::bedtools"
    
    input:
    file windows from generateWindowsFromGtfOut
    file snp from snp_ch
    
    output:
    file 'snp.filteredByWindows.bed' into intersectVariantWithWindowsOut
    
    script:
    """
    bedtools intersect -a $snp -b $windows > snp.filteredByWindows.bed
    """

}

process extractFastaFromSelectSNP {
    conda "bioconda::bedtools"
    
    input:
    file genome from genome_ch 
    file snp from intersectVariantWithWindowsOut
    
    output:
    file 'snp.filteredByWindows.fasta' into intersectVariantWithWindowsFastaOut
    
    script:
    """
    awk '{print \$1"\t"\$2-14"\t"\$2+14;}' $snp > $snp"_modif"
    bedtools getfasta -fi $genome -bed $snp"_modif" | paste - -  > $snp"_modif2"
    paste $snp $snp"_modif2" > snp.filteredByWindows.fasta
    """
}

//generateWindowsFromGtfOut.view()
//intersectVariantWithWindowsOut.view()
intersectVariantWithWindowsFastaOut.view()

