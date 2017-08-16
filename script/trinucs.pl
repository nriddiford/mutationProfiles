#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use feature qw/ say /;
use FindBin '$Script';
# use FindBin qw($Bin);
# use File::Spec;
# use lib File::Spec->catdir($FindBin::Bin, '..', 'lib/');

use vcfParse;

use Bio::SeqIO;
use File::Basename;

use Getopt::Long qw/ GetOptions /;

my $genome_file;
my $vcf_file;
my $snv_dist_file = 'combined_snvs.txt';
my $out_dir = './';
my $help;
my $in_file; # Varscan native

# Should add score threshold option
GetOptions( 'genome=s'				  =>			\$genome_file,
            'infile=s'          =>      \$in_file,
            'vcf=s'						  =>			\$vcf_file,
     			  'out-file=s'        =>	 		\$snv_dist_file,
            'dir=s'             =>      \$out_dir,
     			  'help'         		  =>   		\$help
) or die usage();

if ($help)  { exit usage() }

unless ($in_file or $vcf_file ) { exit usage() }

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my $genome_ref = get_genome($genome_file);
my %genome = %{$genome_ref};
my $data_ref;
$data_ref = parse_vcf($vcf_file) if $vcf_file;
# $data_ref = parse_varscan($in_file) if $in_file;
my ($filtered_data_ref) = get_context($data_ref);
my ($sample, $snv_dist_ref) = count($filtered_data_ref);

# Write to R-friendly dataframe
write_dataframe($sample, $snv_dist_ref);

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

sub parse_vcf {
  my $vcf_file = shift;

  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($vcf_file), 2);
  my ($sample) = split(/_/, $name, 0);

  say "Parsing VCF file...";
  my (@vars, $caller);

  my ($snpData, $info, $filtered_vars ) = vcfParse::parse($vcf_file);

    for ( sort { @{ $snpData->{$a}}[0] cmp @{ $snpData->{$b}}[0] or
          @{ $snpData->{$a}}[1] <=> @{ $snpData->{$b}}[1]
        }  keys %{ $snpData } ){
      my ( $chrom, $pos, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $filters, $samples ) = @{ $snpData->{$_} };

      my (%sample_info)  = @{ $info->{$_}->[6] };

      my $af;

      if ($sample_info{$_}{'TUMOR'}{'AF'}){
        $caller = 'mutect2';
        $af = $sample_info{$_}{'TUMOR'}{'AF'};
      }
      elsif ($sample_info{$_}{'TUMOR'}{'FREQ'}){
        $caller = 'varscan2';
        $af = $sample_info{$_}{'TUMOR'}{'FREQ'};
        ($af) =~ s/%//;
        $af = $af/100;
      }
      push @vars, [$sample, $chrom, $pos, $ref, $alt, $af, $caller];
    }

  return(\@vars);
}


sub get_context {
  my $var_ref = shift;
  my (@filtered_vars, %snvs);

  foreach my $var ( @$var_ref ) {
    my ($sample, $chrom, $pos, $ref, $alt, $af, $caller) = @$var;

   	if ( length $ref == 1 and length $alt == 1 and $chroms{$chrom} ) {

      $snvs{$chrom}{$pos} = [$ref, $alt];

      my ($trinuc) = substr( $genome{$chrom}, $pos - 2, 3 );

      if ($trinuc =~ /N/){
        say "excluding $trinuc";
        next;
      }

      my ($trans_trinuc, $grouped_ref, $grouped_alt) = group_muts($trinuc, $ref, $alt);

      push @filtered_vars, [ $sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $trans_trinuc, $grouped_ref, $grouped_alt ];
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
    my ($samp, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $trans_trinuc, $grouped_ref, $grouped_alt) = @$var;
    $sample = $samp;

    push @snv_dist, [ $sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, "$ref>$alt", $trans_trinuc, "$grouped_ref>$grouped_alt" ];

    # debug($chrom, $pos, $ref, $alt, $trinuc) if $debug;
  }
  return($sample, \@snv_dist);
}

sub write_dataframe {
  my ($sample, $snv_dist_ref) = @_;

  my $outlocation = $out_dir . $snv_dist_file;
  open my $snv_dist, '>>',  $outlocation or die $!;

  say "Printing out genome-wide snv distribution '$outlocation' for $sample...";

  foreach my $var ( @$snv_dist_ref ) {

    my ($sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $trans, $decomp_trinuc, $grouped_trans ) = @$var;
    $sample = $sample;
    print $snv_dist join("\t", $sample, $chrom, $pos, $ref, $alt, $trinuc, $trans, $decomp_trinuc, $grouped_trans, $af, $caller ) . "\n";
  }
}

# sub debug {
#   my ($chr, $pos, $ref, $alt, $trinuc, $mut_cont) = @_;
#   printf "%-30s %-s\n", "SNV:", $ref;
#   printf "%-30s %-s\n", "Position:", "$chr\:$pos";
#   printf "%-30s %-s\n", "Transition:", "$ref>$alt";
#   printf "%-30s %-s\n", "Trinucleotide:", $trinuc;
#   say "***********";
# }

sub usage {
  print
"
usage: $Script [-h] [-v VCF_IN] [-g GENOME]

trinucs
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v0.1
description: Get trinucleotide context from VCF file

arguments:
  -h, --help              show this help message and exit
  -v vcf_file, --vcf      vcf input file
  -g genome, --genome
                          genome fasta file
  -o out-file, --out
                          name of file to write [Default 'combined_snvs.txt']
  -d out-directory, --dir
                          directory to write to [Default cwd]
"
}
