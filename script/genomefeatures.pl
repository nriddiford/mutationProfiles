#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Data::Printer;
use Data::Dumper;
use feature qw/ say /;

use List::Util 'sum';

use Getopt::Long qw/ GetOptions /;

my $chroms;
my $annontation_file = '/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf';


# Should add score threshold option
GetOptions( 'chroms=s'				  =>			\$genome_file,
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

my $gtf_file = $annontation_file;


my @chroms = qw/ 2L 2R 3L 3R X Y 4 /;
my @lengths = qw/ 23513712 25286936 28110227 32079331 23542271 3667352 1348131 /;

my %genome;
@genome{@chroms} = @lengths;

say "Genome:";
say "  $_: $genome{$_}" for sort { $genome{$b} <=> $genome{$a} } keys %genome;
print "\n";

my $gtf_file = '/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf';

my $gene_to_inspect = shift;

open my $gtf_in, '<', $gtf_file;
open my $gene_lengths, '>', 'data/gene_lengths.txt';
open my $tss_pos, '>', 'data/tss_positions.txt';

print $gene_lengths join("\t", 'gene', 'length', 'chrom', 'start', 'end', 'tss', 'scaling_factor', 'id') . "\n";
print $tss_pos join("\t", 'gene', 'chrom', 'tss') . "\n";

my %exons;

my (%features, %transcript_length, %genes, %feature_lengths, %longest_feature_per_gene, %gene_info, %tss_seen);

my (%exons_per_gene);

while(<$gtf_in>){
  chomp;
  my ($chrom, $feature, $start, $stop, $id, $gene) = (split)[0,2,3,4,9,11];

  next unless grep /$chrom/, @chroms;

  ($gene) = $gene =~ /\"(.*)\";/;
  ($id) = $id =~ /\"(.*)\";/;
  my ($gene_length, $transcript,$feature_length);

  next if $feature eq 'CDS' or $feature eq 'mRNA';

  if ($feature ne 'gene'){
    $transcript = (split)[13];
    ($transcript) = $transcript =~ /\"(.*)\";/;
  }

  $feature_length = ($stop - $start);

  if ($feature eq 'gene'){
    $gene_info{$gene}{'gene'} = [ $chrom, $start, $stop, $feature_length, '-', $id ];
    $gene_length = ($stop - $start);
  }

  if ($feature eq 'start_codon'){
    my $tss = (($start+$stop)/2);
    # if there's a tss for this gene, replace '-' at element 4 with 1 element $tss
    splice(@{$gene_info{$gene}{'gene'}}, 4, 1, $tss);
    if (not $tss_seen{$gene}{$tss}++){
      print $tss_pos join("\t", $gene, $chrom, $tss) . "\n";
    }
  }

  if ($feature eq 'exon'){
    $exons{$transcript}{$feature} += $feature_length;
    $feature_length = $exons{$transcript}{$feature};
    $exons_per_gene{$gene}++;
  }

  if ( not exists $longest_feature_per_gene{$gene}{$feature} or ($longest_feature_per_gene{$gene}{$feature} < $feature_length) ){
    $longest_feature_per_gene{$gene}{$feature} = $feature_length;

  }
}

# print Dumper \%tss_locs;
my $total_features_length;

my (%genome_features, %intron_length);
for my $gene ( sort keys %longest_feature_per_gene ){
  my ($gene_length, $per_gene_feature_length);

  $per_gene_feature_length = 0;
  $longest_feature_per_gene{$gene}{'exon'} -= $longest_feature_per_gene{$gene}{'3UTR'} if $longest_feature_per_gene{$gene}{'3UTR'};
  $longest_feature_per_gene{$gene}{'exon'} -= $longest_feature_per_gene{$gene}{'5UTR'} if $longest_feature_per_gene{$gene}{'5UTR'};

  for my $feature ( sort keys $longest_feature_per_gene{$gene} ) {

    if ( $feature ne 'gene' ){

      $per_gene_feature_length += $longest_feature_per_gene{$gene}{$feature};

      $total_features_length += $longest_feature_per_gene{$gene}{$feature};

      $genome_features{$feature} += $longest_feature_per_gene{$gene}{$feature};
    }

    else {
      $gene_length = $longest_feature_per_gene{$gene}{'gene'};
      my $scaling_factor = 1/$gene_length;

      my ( $chrom, $start, $stop, $len, $tss, $id ) = @{$gene_info{$gene}{'gene'}};

      print $gene_lengths join("\t", $gene, $gene_length, $chrom, $start, $stop, $tss, $scaling_factor, $id) . "\n";
    }
  }
  my $intron_length = ($gene_length - $per_gene_feature_length);

  $longest_feature_per_gene{$gene}{'intron'} = $intron_length;

  $genome_features{'intron'} += $intron_length;
  $total_features_length += $intron_length;
}

my $gen_length = sum @lengths;

say "Genome length = " . $gen_length;
say "Total features length = $total_features_length";
print "\n";

open my $genomic_features, ">", 'data/genomic_features.txt';

my $intergenic = ($gen_length - $total_features_length);
$genome_features{'intergenic'} = $intergenic;

print $genomic_features "feature\tlength\tpercentage\n";

say "Feature\tlength\tpercentage";
for (sort { $genome_features{$b} <=> $genome_features{$a} } keys %genome_features){
  print $genomic_features "$_\t$genome_features{$_}\t";
  my $percent = sprintf "%.2f", (($genome_features{$_}/$gen_length)*100) ;
  print $genomic_features "$percent\n";

  print "$_\t$genome_features{$_}\t$percent\n";
}

if ($gene_to_inspect){
  say "$gene_to_inspect";
  say "  $_: $longest_feature_per_gene{$gene_to_inspect}{$_}" for sort { $longest_feature_per_gene{$gene_to_inspect}{$b} <=> $longest_feature_per_gene{$gene_to_inspect}{$a} } keys $longest_feature_per_gene{$gene_to_inspect};
}
