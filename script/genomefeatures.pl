#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Data::Printer;
use Data::Dumper;
use feature qw/ say /;

use List::Util 'sum';

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

print $gene_lengths join("\t", 'gene', 'length', 'chrom', 'start', 'end', 'scaling_factor') . "\n";

my %exons;

my (%features, %transcript_length, %genes, %feature_lengths, %longest_feature_per_gene, %gene_info);

my (%exons_per_gene);

while(<$gtf_in>){
  chomp;
  my ($chrom, $feature, $start, $stop, $gene) = (split)[0,2,3,4,11];

  next unless grep /$chrom/, @chroms;

  ($gene) = $gene =~ /\"(.*)\";/;
  my ($gene_length, $transcript,$feature_length);

  next if $feature eq 'CDS' or $feature eq 'mRNA';

  if ($feature ne 'gene'){
    $transcript = (split)[13];
    ($transcript) = $transcript =~ /\"(.*)\";/;
  }

  $feature_length = ($stop - $start);

  if ($feature eq 'gene'){
    $gene_info{$gene}{'gene'} = [$chrom, $start, $stop, $feature_length ];
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

      my ( $chrom, $start, $stop ) = @{$gene_info{$gene}{'gene'}};

      print $gene_lengths join("\t", $gene, $gene_length, $chrom, $start, $stop, $scaling_factor) . "\n";
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
