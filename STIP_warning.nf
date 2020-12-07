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
    file 'warning.txt' into warning 
    
    script:
    """
    awk '{print \$1"\t"\$2-30"\t"\$2+29;}' $snp > $snp"_modif"
    bedtools getfasta -fi $genome -bed $snp"_modif" | paste - -  > $snp"_modif2"
    paste $snp $snp"_modif2" | awk -v warning='warning.txt' 'BEGIN{w=0}{if (toupper(substr(\$8,30,1)) != toupper(\$5)){w++;} if (substr(\$4,1,2) == "rs"){ name = \$4;}else{name = "snp_"\$1"_"\$2;} print \$1"\t"\$2"\t"\$5"\t"\$6"\t"name"\t"\$8;}END{print w > warning}' > snp.filteredByWindows.fasta
    """
}

//generateWindowsFromGtfOut.view()
//intersectVariantWithWindowsOut.view()
intersectVariantWithWindowsFastaOut.view()
warning.view()
