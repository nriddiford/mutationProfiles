#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Data::Printer;
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
my $out_dir = '.';
my $help;
my $in_file; # Varscan native
my $snpeff = 0;
my $type = 'snv';
my $caller;

# Should add score threshold option
GetOptions( 'genome=s'				  =>			\$genome_file,
            'infile=s'          =>      \$in_file,
            'caller=s'          =>      \$caller,
            'vcf=s'						  =>			\$vcf_file,
            'snpEff'            =>      \$snpeff,
            'type=s'            =>      \$type,
     			  'out-file=s'        =>	 		\$snv_dist_file,
            'dir=s'             =>      \$out_dir,
     			  'help'         		  =>   		\$help
) or die usage();

if ($help)  { exit usage() }

unless ( ($in_file or $vcf_file) and $genome_file) {
  say "Both input file and a genome file required";
  exit usage()
}

my %chroms = qw / 2L 23513712 2R 25286936 3L 28110227 3R 32079331 4 1348131 X 23542271 Y 3667352 /;

my $genome_ref = get_genome($genome_file);
my %genome = %{$genome_ref};
my $data_ref;
my $statements;

($data_ref, $statements) = parse_vcf($vcf_file, $caller) if $vcf_file;

if (scalar %{$statements}){
  p(%{$statements});
}

$data_ref = parse_varscan($in_file) if $in_file;
my ($filtered_data_ref) = get_context($data_ref, $type);
my ($sample, $snv_dist_ref) = count($filtered_data_ref, $type);
# Write to R-friendly dataframe

write_dataframe($sample, $snv_dist_ref, $type);

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
  say $name;
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
  my ($vcf_file, $caller) = @_;
  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($vcf_file), 2);
  my ($sample) = split(/_/, $name, 0);

  say "Parsing VCF file... for sample $sample";
  my @vars;
  my %statements;

  my ($snpData, $info, $filtered_vars, $heads, $sams) = vcfParse::parse($vcf_file);

  for(@{$heads}){
    if (m/mutect/i){
      say "This is a mutect file";
      $caller = 'mutect';
      last;
    }
    elsif (m/varscan/i){
      say "This is a varscan file";
      $caller = 'varscan';
      last;
    }
  }

    for ( sort { @{ $snpData->{$a}}[0] cmp @{ $snpData->{$b}}[0] or
          @{ $snpData->{$a}}[1] <=> @{ $snpData->{$b}}[1]
        }  keys %{ $snpData } ){
      my ( $chrom, $pos, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $filters, $samples ) = @{ $snpData->{$_} };

      my (%info)  = @{ $info->{$_}->[5] };
      my (%sample_info)  = @{ $info->{$_}->[6] };
      my (@samples) = @{ $sams };

      my ($tumour, $normal) = @samples;
      if ($caller eq 'varscan'){
        ($tumour, $normal) = reverse @samples;
      }

      my ($variant_type, $status, $hit_gene, $other) = ("", "", "", "");

      if ($info{$_}{ANN}){
        $statements{'caller'} = "This file has been annotated by SnpEff";
        my @hits = split(/,/, $info{$_}{ANN});
        my @annotated_parts = split(/\|/, $hits[0]);
        ($variant_type, $status, $hit_gene) = @annotated_parts[1..3];
        next unless length $hit_gene;
      }
      my $af;

      if ($caller eq 'varscan' and $sample_info{$_}{$tumour}{FREQ}){
        $af = $sample_info{$_}{$tumour}{FREQ};
        ($af) =~ s/%//;
        $af = $af/100;
      }
      elsif ($caller eq 'mutect' and $sample_info{$_}{$tumour}{AF}){
        $af = $sample_info{$_}{$tumour}{AF};
      }

      push @vars, [$sample, $chrom, $pos, $ref, $alt, $af, $caller, $variant_type, $status, $hit_gene];
    }

  return(\@vars, \%statements);
}


sub get_context {
  my ($var_ref, $type) = @_;

  my (@filtered_vars, %snvs);

  foreach my $var ( @$var_ref ) {
    my ($sample, $chrom, $pos, $ref, $alt, $af, $caller, $variant_type, $status, $hit_gene) = @$var;

   	if ( (length $ref == 1 and length $alt == 1 and $chroms{$chrom}) or $type eq 'indel' ) {

      $snvs{$chrom}{$pos} = [$ref, $alt];

      my ($trinuc) = substr( $genome{$chrom}, $pos - 2, 3 );

      if ($trinuc =~ /N/){
        say "excluding $trinuc";
        next;
      }

      my ($trans_trinuc, $grouped_ref, $grouped_alt) = group_muts($trinuc, $ref, $alt);

      push @filtered_vars, [ $sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $trans_trinuc, $grouped_ref, $grouped_alt, $variant_type, $status, $hit_gene ];
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
  my ($var_ref, $type) = @_;

  my $all_snvs_count = 0;

  my (%genome_wide_snvs, %snvs_by_chrom, %snp_count, %snp_freq, %tri_count);
  my @snv_dist;
  my $sample;

  foreach my $var ( @$var_ref ) {
    my ($samp, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $trans_trinuc, $grouped_ref, $grouped_alt, $variant_type, $status, $hit_gene) = @$var;
    $sample = $samp;

    if ($type eq 'indel'){
      my $mut_type = 'INS';
      if ($alt =~ /^\-/){
        $mut_type = 'DEL';
      }

      push @snv_dist, [ $sample, $chrom, $pos, $ref, "\'$alt\'", $af, $caller, $trinuc, $mut_type, $trans_trinuc, $variant_type, $status, $hit_gene ];
    }
    else{
      push @snv_dist, [ $sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, "$ref>$alt", $trans_trinuc, "$grouped_ref>$grouped_alt", $variant_type, $status, $hit_gene ];
    }


    # debug($chrom, $pos, $ref, $alt, $trinuc) if $debug;
  }
  return($sample, \@snv_dist);
}

sub write_dataframe {
  my ($sample, $snv_dist_ref, $type) = @_;

  my $outlocation = $out_dir . "/" . $snv_dist_file;
  open my $snv_dist, '>>',  $outlocation or die $!;

  # if ( -z $snv_dist ) {
  #   say "Adding header to file";
  # }

  say "Printing out genome-wide snv distribution '$outlocation' for $sample...";

  foreach my $var ( @$snv_dist_ref ) {

    if ($type eq 'indel'){
      my ($sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $mut_type, $trans_trinuc, $variant_type, $status, $hit_gene) = @$var;
      print $snv_dist join("\t", $sample, $chrom, $pos, $ref, $alt, $trinuc, $mut_type, $trans_trinuc, $af, $caller, $variant_type, $status, $hit_gene ) . "\n";
    }
    else {
      my ($sample, $chrom, $pos, $ref, $alt, $af, $caller, $trinuc, $trans, $decomp_trinuc, $grouped_trans, $variant_type, $status, $hit_gene) = @$var;
      print $snv_dist join("\t", $sample, $chrom, $pos, $ref, $alt, $trinuc, $trans, $decomp_trinuc, $grouped_trans, $af, $caller, $variant_type, $status, $hit_gene ) . "\n";
    }


    # snv2gene = $_ + $hit_feature, $hit_gene, $hit_id
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
  -v, --vcf               vcf input file
  -i, --in_file           varscan native file
  -g, --genome            genome fasta file
  -c, --caller            specify caller used ['varscan|mutect']

  -o, --out
                          name of file to write [Default 'combined_snvs.txt']
  -d, --out_dir
                          directory to write to [Default cwd]
  -h, --help              show this help message and exit
"
}
#
#
# GetOptions( 'genome=s'				  =>			\$genome_file,
#             'infile=s'          =>      \$in_file,
#             'caller=s'          =>      \$caller,
#             'vcf=s'						  =>			\$vcf_file,
#             'snpEff'            =>      \$snpeff,
#             'type=s'            =>      \$type,
#      			  'out-file=s'        =>	 		\$snv_dist_file,
#             'dir=s'             =>      \$out_dir,
#      			  'help'         		  =>   		\$help
