#!/usr/bin/perl
use strict;
use Getopt::Long;
#use lib "/work/project/sigenae/Mathieu/STIP/STIPlib/";
#use lib "/work/project/sigenae/Mathieu/STIP/STIPlib/PDL";
use TFBS::Matrix::PFM;


my $help;
my $matrix_file_list;
my $seq_file;
my $matrix_dir;
my $output_imp;
my $output_tfbs;

GetOptions(
	"help"	=> \$help,
	"matrix_list=s"	=> \$matrix_file_list,
	"seq_file=s"	=> \$seq_file,
	"matrix_dir=s" => \$matrix_dir,
	"output_imp=s" => \$output_imp,
	"output_tfbs=s" => \$output_tfbs
);

my $threshold = 0.001;
my $ratio_threshold = 0.05;
my $output_buffer= "";
my $extra_output_buffer = "";

#unless (-d "$out_dir"){
#	system("mkdir $out_dir");
#}

my @matrix_id_list;
my @matrix_db_list;

open (IN, $matrix_file_list) or die ("Can't open $matrix_file_list");
while (my $line =<IN>){
	chomp($line);
	next if $line =~ /^\s*$/;
	my @fields = split(/\t/,$line);
	my $matrix_id = $fields[0];
	my $matrix_db = ($fields[1]) ? $fields[1] : "Unknown";
	push(@matrix_id_list,$matrix_id);
	push(@matrix_db_list,$matrix_db);
}

close IN;


open(OUT, ">$output_imp") or die("Can't write $output_imp\n");					
open(SUP, ">$output_tfbs") or die("Can't write $output_tfbs\n");


for (my $i=0;$i<=$#matrix_id_list;$i++){
#for (my $i=0;$i<=0;$i++){

	my $id = $matrix_id_list[$i];
	my $db = $matrix_db_list[$i];
	
	
	#####################

	my $pwm = &getPWMfromPFMFile("$matrix_dir/$id.pfm");
	print "$id\t",$pwm->ID,"\t",$pwm->length,"\n";
	my $max_profile_width = $pwm->length;

	unless (-f "$matrix_dir/$id.pwm"){
		open(MAT,">$matrix_dir/$id.pwm") or die("Can't write to $matrix_dir/$id.pwm\n");
		print MAT $pwm->rawprint();
		close(MAT);
	}

	#####################

	unless (-f "$matrix_dir/$id".".score_distrib"){
		&create_distrib_score_file("$matrix_dir/$id",$pwm);
	}

	unless (-f "$matrix_dir/$id".".ratio_distrib"){
		&create_distrib_ratio_file("$matrix_dir/$id",$pwm);
	}


	my $score_DTBs = dtb_hashing("score", $pwm, "$matrix_dir/$id".".score_distrib");
	my $ratio_DTBs = dtb_hashing("ratio", $pwm, "$matrix_dir/$id".".ratio_distrib");
    
    #print "Score\t".$score_DTBs."\n"; ###
    #print "Ratio\t".$score_DTBs."\n"; ###

	######################
	open(SEQ,"$seq_file") or die("Can't open $seq_file\n");
	while (my $seq =<SEQ>){
		chomp $seq;	
		# skip blank lines and comment lines
		next unless $seq;
		next if $seq =~ /^#/;

		# splitting line elements

		my @elem = split /\t/, $seq;
		my $chrom = $elem[0];
		my $position = $elem[1];
		my $rs_id = $elem[4];
		my $ref_nlcs = uc $elem[2];
		my $alt_nlcs = uc $elem[3];
		my $ref_seq = uc $elem[5];

		# Format checking 
		if ($ref_nlcs =~ /,/) {print STDERR "WARNING $seq: multiple reference values! Check your VCF file.\n";}

		next if ($ref_seq =~/N/i);
        next if ($ref_nlcs =~/[\.\*\,]/i);
        next if ($alt_nlcs =~/[\.\*\,]/i);

		my $ref = $ref_nlcs;
		my $ref_length = length $ref;

		if ($chrom =~ /^chr(\S+)/) {$chrom = $1;}

		my @alts = split ',', $alt_nlcs;# Possible multiple comma separated alternate alleles.

		# Check for indels
		my $is_indel = 0;
		foreach my $alt (@alts) {
			my $alt_length = length $alt;

			if ($ref_length != $alt_length) {
				$is_indel = 1;
			}
		}
        #print "OK";print "$is_indel\t$ref_nlcs\t$alt_nlcs\n";exit;
		if ($is_indel == 0){
            
			#my $start = $position - $max_profile_width + 1;
			#my $end   = $position + $ref_length + $max_profile_width - 2;
			my $matrix = $pwm;
			my $rel_position = (length($ref_seq)+1)/2;
			my $ref_site_start = $rel_position - $matrix->length + 1;
			my $ref_length = length $ref_nlcs;
			
			my %ref_best_site_plus;
			my %ref_best_site_minus;
			
			my ($ref_best_site_plus,$ref_best_site_minus) = &best_stranded_site_quick($matrix,$ref_seq,$ref_site_start,$score_DTBs);
            
            #print join("\t",$matrix,$ref_seq,$ref_site_start,$score_DTBs,"\n");
            #print join("\t",$ref_best_site_plus,$ref_best_site_minus,"\n");

			foreach my $alt (@alts) {
				
				#my $id = $matrix->ID;



				my $alt_seq = $ref_seq;
				substr($alt_seq, $rel_position - 1, $ref_length, $alt);
				my $alt_site_start = $ref_site_start;

				if ($ref_seq!~/^[ATGC]+$/i){next;}


				my ($alt_best_site_plus,$alt_best_site_minus) = &best_stranded_site_quick($matrix,$alt_seq,$alt_site_start,$score_DTBs);
				if ($alt_seq!~/^[ATGC]+$/i){next;}
				
				#print $id,"\t",$ref_nlcs,"\t",$alt,"\t",$ref_seq,"\n";
				#print "\n";
				
				#print join("\t",$ref_best_site_plus->{"seq"},$ref_best_site_plus->{"start"},$ref_best_site_plus->{"end"},$ref_best_site_plus->{"score"},$ref_best_site_plus->{"pvalue"},"\n");
				#print join("\t",$ref_best_site_minus->{"seq"},$ref_best_site_minus->{"start"},$ref_best_site_minus->{"end"},$ref_best_site_minus->{"score"},$ref_best_site_minus->{"pvalue"},"\n");
				
				#print join("\t",$alt_best_site_plus->{"seq"},$alt_best_site_plus->{"start"},$alt_best_site_plus->{"end"},$alt_best_site_plus->{"score"},$alt_best_site_plus->{"pvalue"},"\n");
				#print join("\t",$alt_best_site_minus->{"seq"},$alt_best_site_minus->{"start"},$alt_best_site_minus->{"end"},$alt_best_site_minus->{"score"},$alt_best_site_minus->{"pvalue"},"\n");
					
				predict_impact_quick($ref_best_site_plus, $alt_best_site_plus, $chrom, $position, $rs_id, $ref_nlcs, $rel_position, $alt, $matrix, "+",$id,$db,$score_DTBs,$ratio_DTBs);
				predict_impact_quick($ref_best_site_minus, $alt_best_site_minus, $chrom, $position, $rs_id, $ref_nlcs, $rel_position, $alt, $matrix, "-",$id,$db,$score_DTBs,$ratio_DTBs);
					
			}
			
			
			print OUT $output_buffer;
			$output_buffer = "";

			print SUP $extra_output_buffer;
			$extra_output_buffer = "";


		}
	}
	close SEQ;
}

close OUT;
close SUP;



sub predict_impact_quick{
	my ($ref_best_site, $alt_best_site, $chrom, $position, $rs_id, $ref_nlcs, $rel_position, $alt, $matrix, $strand,$id,$db,$score_DTBs,$ratio_DTBs) = @_;
	# Return if both score p-values are above the score p-value threshold
	return 1 if ($ref_best_site->{"pvalue"} > $threshold && $alt_best_site->{"pvalue"} > $threshold);
	
	my $ref_best_site_start = $position - $rel_position + $ref_best_site->{"start"};
	my $ref_best_site_end = $position - $rel_position + $ref_best_site->{"end"};
	my $alt_best_site_start = $position - $rel_position + $alt_best_site->{"start"};
	my $alt_best_site_end = $position - $rel_position + $alt_best_site->{"end"};
	# Report detected TFBS if extra option is specified
	$extra_output_buffer .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\n",
	$chrom,
	$position,
	$rs_id,
	$ref_nlcs,
	$alt,
	$matrix->name,$db,
	$strand,
	$ref_best_site_start."-".$ref_best_site_end.":".$ref_best_site->{"seq"},
	$ref_best_site->{"pvalue"},
	0,#100*$ref_best_site->rel_score
	$alt_best_site_start."-".$alt_best_site_end.":".$alt_best_site->{"seq"},
	$alt_best_site->{"pvalue"},
	0,#100*$alt_best_site->rel_score,
	$id;
	
	# Compute the ratio and convert it to p-value		
	my $ratio = ratio($ref_best_site->{"pvalue"}, $alt_best_site->{"pvalue"});
	my $ratio_pvalue = $ratio_DTBs->{$matrix->ID}{$ratio};
	my $status;
	
	# Return if ratio p-value is above the ratio p-value threshold.
	# Otherwise, depending on the ratio, specify the loss or gain state.
	if ($ref_best_site->{"pvalue"} < $alt_best_site->{"pvalue"} && (1-$ratio_pvalue) < ($ratio_threshold/2)){
		$status = "Loss";
		$ratio_pvalue = 1-$ratio_pvalue;
	} elsif ($ref_best_site->{"pvalue"} > $alt_best_site->{"pvalue"} && $ratio_pvalue < ($ratio_threshold/2)){
		$status = "Gain";
	}
	else {
		return 1;
	}
	
	# Report detected TFBS affected by the variant in the output file
	$output_buffer .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\n",
	$chrom,
	$position,
	$rs_id,
	$ref_nlcs,
	$alt,
	$matrix->name,$db,
	$status,
	$ratio_pvalue,
	$strand,
	$ref_best_site_start."-".$ref_best_site_end.":".$ref_best_site->{"seq"},
	$ref_best_site->{"pvalue"},
	0,#100*$ref_best_site->rel_score
	$alt_best_site_start."-".$alt_best_site_end.":".$alt_best_site->{"seq"},
	$alt_best_site->{"pvalue"},
	0,#100*$alt_best_site->rel_score,
	$id;
	
	#print join("\t",$chrom,$position,$rs_id,$ref_nlcs,$alt,$matrix->name,$db,$status,$ratio_pvalue,$strand,$ref_best_site_start."-".$ref_best_site_end.":".$ref_best_site->{"seq"},$ref_best_site->{"pvalue"},0,$alt_best_site_start."-".$alt_best_site_end.":".$alt_best_site->{"seq"},$alt_best_site->{"pvalue"},0,$id,"\n");
	#exit;
}

# Function to find closest numeric hash key to a score
sub findClosestInHash{
	my $DTB_ref = shift;
	my $id = shift;
	my $score = shift;
	my $pvalue;
	
	my $upValue;
	my $downValue = -2000;
	my $closest;
	
	foreach my $DTB_score (sort {$a <=> $b} (keys(%{$DTB_ref->{$id}}))){
		$upValue = $DTB_score;
		if ($score <= $upValue && $score >= $downValue){
			$closest = findClosest($upValue, $downValue, $score);
		} 
		else {
			$downValue = $upValue;
		}
	}
	unless ($closest){
		$closest = $downValue;
	}
	return $DTB_ref->{$id}{$closest};
}

# Function to determine which one of two values is the closest to 
# a third value
sub findClosest{
	my($up, $down, $value) = @_;
	
	if (abs($up-$value) <= abs($down-$value)){
		return $up;
	}
	else {
		return $down;
	}	
}

sub ratio {
	my $ref = shift;
	my $alt = shift;
	
	my $ratio = round(.1,(log2($ref/$alt)));
	return $ratio;
}

sub get_Matrix_Score_Hash{
	#my $pwm = shift; 
	my @matrix = split(/\n/,shift);
	my %matrix_score;
	my $line_score_A = $matrix[0];
	$line_score_A =~ s/^\s+|\s+$//g; ##PAs oublié le trim avant!!! dans le cas d'un score positif, il y a un espace devant la valeur dans la pwm
	my @scoreA = split(/\s+/,$line_score_A);
	$matrix_score{"A"} = \@scoreA;
	
	my $line_score_C = $matrix[1];
	$line_score_C=~ s/^\s+|\s+$//g;
	my @scoreC = split(/\s+/,$line_score_C);
	$matrix_score{"C"} = \@scoreC;
	
	my $line_score_G = $matrix[2];
	$line_score_G=~ s/^\s+|\s+$//g;
	my @scoreG = split(/\s+/,$line_score_G);
	$matrix_score{"G"} = \@scoreG;
	
	my $line_score_T = $matrix[3];
	$line_score_T=~ s/^\s+|\s+$//g;
	my @scoreT = split(/\s+/,$line_score_T);
	$matrix_score{"T"} = \@scoreT;
	
	return \%matrix_score;
}

sub best_stranded_site_quick{
	my %best_site_plus;
	my %best_site_minus;
	my $matrix = shift;
	my $ref_seq = shift;
	my $ref_site_start = shift;
	my $score_DTBs = shift;
	
	#print "Best ", $matrix->ID,"\n";
	
	my $matrix_score = &get_Matrix_Score_Hash($matrix->rawprint);
	
	for (my $i=0;$i<=$matrix->length-1;$i++){
		my $current_start = $ref_site_start+$i; 
		my $current_end = $ref_site_start+ $i + $matrix->length -1;
		my $current_seq = substr($ref_seq,$current_start-1,$matrix->length);
		my $current_seq_rev="";
		
		my $score_plus = 0;
		my $score_minus = 0;

		my @seq = split(//,$current_seq);

		for (my $j=0;$j<=$#seq;$j++){
			if ($seq[$j] eq "A"){
				$score_plus += $matrix_score->{"A"}->[$j];
				$current_seq_rev = "T".$current_seq_rev;
				$score_minus += $matrix_score->{"T"}->[$#seq-$j];
			}
			elsif ($seq[$j] eq "C"){
				$score_plus += $matrix_score->{"C"}->[$j];
				$current_seq_rev = "G".$current_seq_rev;
				$score_minus += $matrix_score->{"G"}->[$#seq-$j];
			}
			elsif ($seq[$j] eq "G"){
				$score_plus += $matrix_score->{"G"}->[$j];
				$current_seq_rev = "C".$current_seq_rev;
				$score_minus += $matrix_score->{"C"}->[$#seq-$j];
			}
			elsif ($seq[$j] eq "T"){
				$score_plus += $matrix_score->{"T"}->[$j];
				$current_seq_rev = "A".$current_seq_rev;
				$score_minus += $matrix_score->{"A"}->[$#seq-$j];
			}
			else{
				die("Erreur 1\n".$seq[$j]."\n$j\n".$current_seq."\n".$ref_seq."\t".$current_start."\t".$matrix->length."\n");
			}
		}
		
		my $pvalue_plus = sc2pv($score_DTBs, $matrix->ID, round(.1, $score_plus), $matrix);
		my $pvalue_minus = sc2pv($score_DTBs, $matrix->ID, round(.1, $score_minus), $matrix);
		if ((!%best_site_plus)||($score_plus > $best_site_plus{"score"})){
			$best_site_plus{"seq"} = $current_seq;
			$best_site_plus{"score"} = round(.1, $score_plus);
			$best_site_plus{"start"} = $current_start; 
			$best_site_plus{"end"} = $current_end;
			$best_site_plus{"pvalue"} = sc2pv($score_DTBs, $matrix->ID, round(.1, $score_plus), $matrix);

		}
		if ((!%best_site_minus)||($score_minus > $best_site_minus{"score"})){
			$best_site_minus{"seq"} = $current_seq_rev;
			$best_site_minus{"score"} = round(.1, $score_minus);
			$best_site_minus{"start"} = $current_start; 
			$best_site_minus{"end"} = $current_end;
			$best_site_minus{"pvalue"} = sc2pv($score_DTBs, $matrix->ID, round(.1, $score_minus), $matrix);
		}
		

		#print $current_seq,"\t",$score_plus,"\t",$pvalue_plus,"\t",$current_start,"\t",$current_end,"\n";
		#print $current_seq_rev,"\t",$score_minus,"\t",$pvalue_minus,"\t",$current_start,"\t",$current_end,"\n";
        #exit;

	}
	
	return 	\%best_site_plus,\%best_site_minus;
}

sub sc2pv{
	my $DTB_ref = shift;
	my $id = shift;
	my $score = shift;
	my $matrix = shift;
	my $pvalue;
	my $pvalueTFM;
	
	if($DTB_ref->{$id}{$score}){
		$pvalue = $DTB_ref->{$id}{$score};
	} else {
		$pvalue  = findClosestInHash($DTB_ref, $id, $score);
	}
	return $pvalue;
}

sub dtb_hashing
{
	my $type = shift;
	my $matrix = shift;
	my $DTB_file_name = shift;
	my %distribH;
	
	my $id = $matrix->ID;
	#my $DTB_file_name = "$id"."_"."$type"."_distrib";
	
	open(DTB,"$DTB_file_name") or die ("Couln't open $DTB_file_name\n");
						
	while (my $DTB_line = <DTB>){
		chomp $DTB_line;
		my @DTB_elem = split(/ /,$DTB_line);
		my $DTB_score = $DTB_elem[0];
		my $DTB_pvalue = $DTB_elem[1];
		$distribH{$id}{$DTB_score} = $DTB_pvalue;
	}	
	close(DTB);
	
	return \%distribH;		
}

sub create_distrib_score_file{
	my $file = shift;
	my $out_string = `TFMpvalue-distrib -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m $file.pwm -s -100 -S 100 -G 0.0001 -w`;
	my @out_lines = split(/\n/,$out_string);
	my $lines_nb = @out_lines;
	my %score_pv_tab;
	my %score_freq_tab;
	for (my $i; $i < $lines_nb -1 ; $i++){
		my @elem = split(/\ /,$out_lines[$i]);
		my $score = round(.1, $elem[0]);
		if ($score_pv_tab{$score}){
			if ($score_pv_tab{$score} < $elem[1]){
				$score_pv_tab{$score} = $elem[1];
			}
		}else{
			$score_pv_tab{$score} = $elem[1];
		}
		$score_freq_tab{$score} += $elem[2];
	}
	open(DTB,">$file.score_distrib");
	foreach my $score (sort {$a <=> $b}(keys(%score_pv_tab))){
		printf DTB "%s %s %s\n",$score, $score_pv_tab{$score}, $score_freq_tab{$score};
	}
	close(DTB);
}

sub create_distrib_ratio_file{
	my $file = shift;
	my $pwm = shift;
	my %freq_tab;
	my %pval_tab;
	my %ratio_tab;
	my $line_number = 0;
	open(DTB,"$file.score_distrib");
	# reading the score distribution and storing informations
	while (my $line = <DTB>){
		chomp $line;
		$line_number++;
		my @elem = split(' ',$line);
		$pval_tab{$line_number} = $elem[1];
		$freq_tab{$line_number} = $elem[2]*(4**$pwm->length);
	}
	# The double foreach loop
	foreach my $line1 (keys(%pval_tab)){
		my $pval1 = $pval_tab{$line1};
		my $freq1 = $freq_tab{$line1};
		foreach my $line2 (keys(%pval_tab)){
			my $freq2 = $freq_tab{$line2};
			my $pval2 = $pval_tab{$line2};
			
			my $ratio = (log2($pval1/$pval2));
			# like the score distrib, the ratio distrib is rounded
			$ratio = round(.1,$ratio);
			# The number of times the division happens
			$ratio_tab{$ratio} += $freq1 * $freq2;
		}
	}
	close(DTB);
	#
	# We have the ratio distribution and the frequency of each ratio. We
	# can compute the pvalue of each ratio using the following formula 
	# pvalue = (total_number_of_ratios - 
	# number_of_ratios_under_the_observed_one) / total_number_of_ratios
	#
	open(DTB,">$file.ratio_distrib");
	my $total_ratio = 16**($pwm->length);
	my $somme = 0;
	foreach my $ratio (sort {$a <=> $b}(keys(%ratio_tab))){
		my $pval = ($total_ratio - $somme) / $total_ratio;
		$somme += $ratio_tab{$ratio};
		print DTB $ratio," ",$pval," ", ($ratio_tab{$ratio}/$total_ratio),"\n";
	}
	close(DTB);
}

sub log2 {
	my $n = shift;
	return log($n)/log(2);
}
#
# For a given matrix, round every value within it to the given granularity
#
sub round_matrix{
	my $granularity = shift;
	my $matrix = shift;
	my @mat_array = $matrix->matrix;
	my @new;
	
	for (my $i = 0; $i < 4; $i++){
		for (my $j = 0; $j < $matrix->length; $j++){
			$new[0][$i][$j] = round($granularity, $mat_array[0][$i][$j]);
		}
	} 
	$matrix->set_matrix(@new);
	return $matrix;
}

#
# For a given value, round it to the given granularity
#
sub round{
	my $granularity = shift;
	my $value = shift;
	my $rounded_value = $granularity * int($value/$granularity)
}

sub getPWMfromPFMFile{
	my $PFM_file = shift;
	my $pwm;
	
	open(FH, $PFM_file) or die ("Can't open $PFM_file\n");

	my $id 			  = '';
	my $name		  = '';
	my $matrix_string = '';
	my $line_count	= 0;

	while (my $line = <FH>) {
		chomp $line;
		next if !$line;
		if ($line =~ /^>(\S+)\s(\S+)/g) {
			$id = $1;
			$name = $2;
		}
		else {
			if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
				# line of the form: A [ # # # ... # ]
				$matrix_string .= "$1\n";
			} elsif ($line =~ /^\s*\d+/) {
				# line of the form: # # # ... #
				$matrix_string .= "$line\n";
			} else {
				next;
			}
			$line_count++;

			if ($line_count == 4) {
				my $pfm = TFBS::Matrix::PFM->new(
					-matrixstring   => $matrix_string,
					-name		=> $name,
					-ID	        => $id
				);
				
				$pwm = round_matrix(.0001, $pfm->to_PWM);
				return $pwm;	
			}
		}
	}
	close(FH);
	return $pwm;
}
