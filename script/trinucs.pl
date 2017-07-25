#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use feature qw/ say /;
use FindBin '$Script';
use Bio::SeqIO;
use File::Basename;

use Getopt::Long qw/ GetOptions /;

my $genome_file;
my $vcf_file;
my $chrom_out_file = 'chroms.trinucs.txt';
my $genome_out_file = 'GW.trinucs.txt';
my $snv_dist_file = 'GW.snv.dist.txt';
my $debug;
my $quiet;
my $help;
my $in_file;

# Should add score threshold option
GetOptions( 'genome=s'				=>			\$genome_file,
            'infile=s'        =>      \$in_file,
            'vcf=s'						=>			\$vcf_file,
     			  'all-snvs=s'    	=>	 		\$genome_out_file,
     			  'chrom-snvs=s'    =>	 		\$chrom_out_file,
     			  'help'         		=>   		\$help,
     			  'quiet'        		=>   		\$quiet,
           	'debug'        		=>    	\$debug
) or die usage();

if ($help)  { exit usage() }
if ($quiet) { say "Running in quiet mode" }
if ($debug) { say "Running in debug mode" }

unless ($in_file or $vcf_file ) { exit usage() }

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my $genome_ref = get_genome($genome_file);

my %genome = %{$genome_ref};

my $data_ref;

$data_ref = parse_vcf($vcf_file) if $vcf_file;
$data_ref = parse_varscan($in_file) if $in_file;

my ($filtered_data_ref) = get_context($data_ref);
my ($sample, $all_snvs_count, $genome_wide_snvs_ref, $snvs_by_chrom_ref, $snp_count_ref, $snp_freq_ref, $tri_count_ref, $snv_dist_ref) = count($filtered_data_ref);

# write_snvs_per_chrom($chrom_out_file, $sample, $snvs_by_chrom_ref, $snp_count_ref);
# write_snvs_genome_wide($all_snvs_count, $sample, $genome_wide_snvs_ref);
write_dataframe($all_snvs_count, $sample, $snv_dist_ref);


sub get_genome {
  my $genome_file = shift;
  my $seqio = Bio::SeqIO->new('-file' => "$genome_file", '-format' => 'Fasta');

  my $genome_length;

  for (keys %chroms){
    ($genome_length) += $chroms{$_};
  }

  say "Reading in genome: $genome_file";

  while(my $seq = $seqio->next_seq) {
    my $nucs = $seq->seq;
    my $chr = $seq->id;
    next unless $chroms{$chr};
    $genome{$chr} = $nucs;
  }
  return(\%genome);
}


sub parse_varscan {
  my $in_file = shift;
  say "Reading in varscan file: $in_file";
  open my $VAR_in, '<', $in_file or die $!;

  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($in_file), 2);
  my ($sample) = split(/_/, $name, 0);

  say "Parsing varscan native file...";
  my @vars;
  while(<$VAR_in>){
    chomp;
    my ($chr, $pos, $ref, $alt, $n_freq, $t_freq, $type) = (split)[0,1,2,3,6,10,12];
    push @vars, [$sample, $chr, $pos, $ref, $alt, $n_freq, $t_freq, $type];
  }
  return(\@vars);
}

my %data;

sub parse_vcf {
  my $vcf_file = shift;

  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($vcf_file), 2);
  my $VCF_in;
  my ($sample) = split(/_/, $name, 0);

  if ( $extention eq 'gz'){
    say "Reading in compressed VCF file: $vcf_file";
    open $VCF_in, "zcat $vcf_file |" or die $!;
  }
  else{
    say "Reading in VCF file: $vcf_file";
    open $VCF_in, '<', $vcf_file or die $!;
  }

  say "Parsing VCF file...";
  my @vars;

  while(<$VCF_in>){
    chomp;
    next if /^#/;
    my ($chr, $pos, $ref, $alt) = (split)[0,1,3,4];
    next if $ref eq 'N';
    my $type = 'somatic';
    if ($_ =~ /VT=(.*?);/){
      $type = $1;
    }
    push @vars, [$sample, $chr, $pos, $ref, $alt, $type];
  }
  return(\@vars);
}


sub get_context {
  my $var_ref = shift;
  my (@filtered_vars, %snvs);

  foreach my $var ( @$var_ref ) {
    my ($sample, $chr, $pos, $ref, $alt, $type) = @$var;

   	if ( length $ref == 1 and length $alt == 1 and $chroms{$chr} ) {

      $snvs{$chr}{$pos} = [$ref, $alt];

      my ($trinuc) = substr( $genome{$chr}, $pos - 2, 3 );

      if ($trinuc =~ /N/){
        say "excluding $trinuc";
        next;
      }

      my ($trans_trinuc, $grouped_ref, $grouped_alt) = group_muts($trinuc, $ref, $alt);

      push @filtered_vars, [$sample, $chr, $pos, $ref, $alt, $trinuc, $trans_trinuc, $grouped_ref, $grouped_alt, $type];
    }
  }
  return(\@filtered_vars);
}


sub group_muts {
  my ($trinuc, $ref, $alt) = @_;
  my ($new_ref, $new_alt) = ($ref, $alt);

  if ($ref eq 'G'){
    $new_ref = 'C';
    $new_alt = 'A' if $alt eq 'T';
    $new_alt = 'G' if $alt eq 'C';
    $new_alt = 'T' if $alt eq 'A';
    $trinuc = rev_comp($trinuc);
  }
  elsif ($ref eq 'A'){
    $new_ref = 'T';
    $new_alt = 'A' if $alt eq 'T';
    $new_alt = 'C' if $alt eq 'G';
    $new_alt = 'G' if $alt eq 'C';
    $trinuc = rev_comp($trinuc);
  }

  return($trinuc, $new_ref, $new_alt);
}

sub rev_comp {
  my $trinuc = shift;
  my $rev_tri = reverse($trinuc);
  # complement the reversed DNA sequence
  $rev_tri =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
  return $rev_tri;
}


sub count {
  my ($var_ref) = shift;

  my $all_snvs_count = 0;

  my (%genome_wide_snvs, %snvs_by_chrom, %snp_count, %snp_freq, %tri_count);
  my @snv_dist;
  my $sample;

  foreach my $var ( @$var_ref ) {
    my ($samp, $chr, $pos, $ref, $alt, $trinuc, $trans_trinuc, $grouped_ref, $grouped_alt, $type) = @$var;
    $sample = $samp;

    $genome_wide_snvs{$trinuc}{"$ref>$alt"}++;		# count genome-wide trinuc transformations: "A>C"

    $snvs_by_chrom{$chr}{$trinuc}{"$ref>$alt"}++;	# count chromosome-wide trinuc transformations: 2L "A>C"

    $snp_freq{$chr}{$ref}{$alt}++; 					      # count total transformations by type (e.g. `A` -> `G`) per chromosome
    $snp_count{$chr}++; 							            # count total transformations per chromosome
    $all_snvs_count++; 								            # count total transformations
    $tri_count{$chr}++;								            # count total trinucs per chromosome

    # Record location of each snv per sample
    push @snv_dist, [ $sample, $chr, $pos, $ref, $alt, $trinuc, "$ref>$alt", $trans_trinuc, "$grouped_ref>$grouped_alt", $type ];

    my $snp_count = $snp_freq{$chr}{$ref}{$alt};
    # my ($mut_cont) = eval sprintf('%.1f', $snp_count/$snp_count{$chr} * 100);

    debug($chr, $pos, $ref, $alt, $trinuc) if $debug;

  }
  return($sample, $all_snvs_count, \%genome_wide_snvs, \%snvs_by_chrom, \%snp_count, \%snp_freq, \%tri_count, \@snv_dist);
}


sub write_snvs_per_chrom {
  my ($chrom_out_file, $sample, $snvs_by_chrom_ref, $snp_count_ref) = @_;

  my %snvs_by_chrom = %{ $snvs_by_chrom_ref};
  my %snp_count = %{ $snp_count_ref};

  open my $chrom_out, '>>',  $chrom_out_file or die $!;
  say "Printing out snvs per chromosome to '$chrom_out_file' for $sample...";

  for my $chr (sort keys %snvs_by_chrom){
    for my $tri (sort keys %{$snvs_by_chrom{$chr}} ) {
      for my $ref_alt (sort keys %{$snvs_by_chrom{$chr}{$tri}} ){
        my ($count) = $snvs_by_chrom{$chr}{$tri}{$ref_alt};
        my ($freq) = eval sprintf('%.2f', ( $count/$snp_count{$chr} ) * 100);
        print $chrom_out join("\t", $chr, $tri, $ref_alt, $freq, $sample) . "\n";
        # print join("\t", $chr, $tri, $ref_alt, $freq, $sample) . "\n";
      }
    }
  }
  say "...done";
}


sub write_snvs_genome_wide {
  my ($all_snvs_count, $sample, $genome_wide_snvs_ref) = @_;
  my %genome_wide_snvs = %{ $genome_wide_snvs_ref };

  open my $genome_out, '>>',  $genome_out_file or die $!;

  say "Printing out genome-wide snvs '$genome_out_file' for $sample...";

  for my $tr (sort keys %genome_wide_snvs){
    for my $ref_alt (sort keys %{ $genome_wide_snvs{$tr} } ) {
      my ($count) = $genome_wide_snvs{$tr}{$ref_alt};
      my ($freq) = eval sprintf('%.2f', ( $count/$all_snvs_count ) * 100);
      print $genome_out join("\t", $tr, $ref_alt, $freq, $sample) . "\n";
      # print join("\t", $tr, $ref_alt, $freq, $sample) . "\n";
    }
  }
  say "...done";
}

sub write_dataframe {
  my ($all_snvs_count, $sample, $snv_dist_ref) = @_;
  open my $snv_dist, '>>',  $snv_dist_file or die $!;

  say "Printing out genome-wide snv distribution '$snv_dist_file' for $sample...";

  foreach my $var ( @$snv_dist_ref ) {
    my ($samp, $chr, $pos, $ref, $alt, $trinuc, $trans, $decomp_trinuc, $grouped_trans, $type) = @$var;
    $sample = $samp;
    print $snv_dist join("\t", $sample, $chr, $pos, $ref, $alt, $trinuc, $trans, $decomp_trinuc, $grouped_trans, $type) . "\n";
  }
}

sub debug {
  my ($chr, $pos, $ref, $alt, $trinuc, $mut_cont) = @_;
  printf "%-30s %-s\n", "SNV:", $ref;
  printf "%-30s %-s\n", "Position:", "$chr\:$pos";
  printf "%-30s %-s\n", "Transition:", "$ref>$alt";
  printf "%-30s %-s\n", "Trinucleotide:", $trinuc;
  say "***********";
}

sub usage {
  print
"
usage: $Script [-h] [-v VCF_IN] [-g GENOME] [-a ALLSNVS] [-c CHROMSNVS] [-q QUIET] [-d DEBUG]

trinucs
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v0.1
description: Get trinucleotide context from VCF file

arguments:
  -h, --help            show this help message and exit
  -v VCF_IN, --vcf      vcf input file
  -g GENOME, --genome
                        genome fasta file
  -a ALLSNVS, --all-snvs
                        specify name of output file for all snvs
  -c CHROMSNVS, --chrom-snvs
                        specify name of output file for snvs per chromosome
  -q QUIET, --quiet     run in quet mode
  -d DEBUG, --debug     run in debug mode
"
}
