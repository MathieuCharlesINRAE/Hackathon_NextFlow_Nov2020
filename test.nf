#!/usr/bin/env nextflow

params.seq='snp.formatted.fasta'
params.matrix_batch='pfm_0000'
params.maxtrix_dir='Hackathon_NextFlow_Nov2020/matrix_stat'
params.script='test.pl'

seq_ch = Channel.fromPath(params.seq)
matrix_ch = Channel.fromPath(params.matrix_batch)
matrixDir_ch = Channel.fromPath(params.maxtrix_dir)
script_ch = Channel.fromPath(params.script)

process toto {
    conda 'STIPenv.yaml'
    echo true
    //conda "eumetsat::pdl"
    //conda "bioconda::perl-bioperl"
    //conda 'STIPenv.yaml'
    
    input:
    file seq from seq_ch
    file matrix from matrix_ch
    file matrixDir from matrixDir_ch
    file script from script_ch
    
    output:
    file 'imp_0000' into imp
    
    script:
    """
    export PERL5LIB=/home/ec2-user/environment/Hackathon_NextFlow_Nov2020/TFBS-0.7.1/:$PERL5LIB
    echo "test"
    perl $script > imp_0000
    """
    
}

imp.view()
//perl STIP_hackathon.pl -matrix_list $matrix -seq_file $seq -matrix_dir $matrixDir -output_imp 'imp_0000'
// /home/ec2-user/environment/work/conda/stip2-a3525cbe5a0dd3c600f5bfb61da26d99/lib/site_perl/5.26.0/x86_64-linux-thread-multi/:$PERL5LIB
// /home/ec2-user/environment/work/conda/stip2-a3525cbe5a0dd3c600f5bfb61da26d99/lib/site_perl/5.26.2/:$PERL5LIB
    