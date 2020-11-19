# Hackathon_NextFlow_Nov2020

**Goal**

The goal is to identified SNP that have impact on transcription factor binding site.

* inputs
    - Genome reference file: Fasta and annotation in GTF file
    - Variant file : in VCF format
    - TFBS pattern matrices : directory of PWM files (one file per TFBS)

* outputs
    - TSV file with regulatory variant, and metadata such as TFBS name, score of affinity, impacting score

* How it works

step1:
select variant in regulatory regions. We could either provide known region or by default STIP will extract 2kb upstream region from genes start.

step2:
extract flanking sequences 14 bp from both side.

step3: 
conversion of TFBS PWM matrices into PFM matrices.

step4: 
each sequence (and its mutated version) in each strand will be compared to each PFM matrice. An afinitiy signal is computed and if the signal is significant a ratio is compute to evaluate the impact of the variant on this TFBS

