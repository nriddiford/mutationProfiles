#!/usr/bin/perl

use strict;
use warnings;

use feature qw/ say /;
use Data::Dumper;
use Data::Printer;
use File::Basename;
use autodie;

my $vcf_file = shift;

my ($snvs, $info, $filtered_data) = parse($vcf_file);

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

my ($id) = (split(/_/, $name))[0];

open my $out, '>', $id . "_filtered.snv.vcf";

my %filtered_snvs = %{ $filtered_data };
my $sv_count = 0;

for (sort {$a <=> $b} keys %filtered_snvs){
  my $line = $filtered_snvs{$_};

  if ($line =~ /^#/){
    print $out $line . "\n";
  }
  else {
    $sv_count++;
    my @cols = split("\t", $line);
    print $out join("\t", @cols[0..5], "PASS", @cols[7..$#cols]) . "\n";
  }
}

say "$sv_count snvs in $id";

sub parse {
  my $file = shift;
  open my $in, '<', $file or die $!;

  my @headers;

  my (%snvs, %info, %filtered_snvs);
  my ($tumour_name, $control_name);
  my %format_long;
  my %info_long;
  my $filter_count;
  my @samples;
  my $replacement_id = 1;

  while(<$in>){
    chomp;

    if (/^#{2}/){
       push @headers, $_;
      $filtered_snvs{$.} = $_;

      if (/##FORMAT/){

        my ($format_long) = $_ =~ /\"(.*?)\"/;
        my ($available_format_info) = $_ =~ /ID=(.*?),/;
        $format_long{$available_format_info} = $format_long;
      }

      if (/##INFO/) {
        my ($info_long) = $_ =~ /\"(.*?)\"/;
        my ($available_info) = $_ =~ /ID=(.*?),/;
        $info_long{$available_info} = $info_long;
      }
      next;
    }

    if (/^#{1}/){
       push @headers, $_;
      $filtered_snvs{$.} = $_;
      my @split = split;
      push @samples, $_ foreach @split[9..$#split];

      $control_name = $samples[0];
      $tumour_name = $samples[1];
      next;
    }

    my @fields = split;

    my ($chr, $pos, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, @sample_info) = @fields;

    if ( $id eq '.' ){
      $id = $replacement_id++;
    }

    my %sample_parts;

    push @{$sample_parts{$samples[$_]}}, split(/:/, $sample_info[$_]) for 0..$#samples;

    my @normal_parts   = split(/:/, $sample_info[0]);
    my @tumour_parts   = split(/:/, $sample_info[1]);

    my @format        = split(/:/, $format_block);
    my @info_parts    = split(/;/, $info_block);

    my %sample_info;

    for my $sample (@samples){
      for my $info (0..$#format){
        $sample_info{$id}{$sample}{$format[$info]} = $sample_parts{$sample}[$info];
      }
    }

    my @filter_reasons;

    #######################
    # Allele depth filter #
    #######################

    # Filter if alt AD < 2
    if ($sample_info{$id}{'TUMOR'}{'AD'}){
      my ($ref_ad, $alt_ad) = split(/,/, $sample_info{$id}{'TUMOR'}{'AD'});
      if ($alt_ad < 2){
        # say "Alt allele depth filt: $alt_ad";
        push @filter_reasons, "alt\_AD=" . $sample_info{$id}{'TUMOR'}{'AD'};
      }

      # Filter if alt AD + ref AD < 5
      if ($alt_ad + $ref_ad < 5){
        my $sample_depth = $alt_ad + $ref_ad;
        # say "Sample depth filt: $sample_depth";
        push @filter_reasons, "AD=" . $sample_depth;
      }
      # Filter if normal alt AD != 0
      my ($n_ref_ad, $n_alt_ad ) = split(/,/, $sample_info{$id}{'NORMAL'}{'AD'});
      if ($n_alt_ad != 0){
        # say "Normal alt allele depth filt: $n_alt_ad";
        push @filter_reasons, "NORM_alt_AD=" . $n_alt_ad;
      }
    }

    if ($sample_info{$id}{'TUMOR'}{'AF'}){
      my $af = $sample_info{$id}{'TUMOR'}{'AF'};
      if ($af < 0.075){
        # say "Allele freq filt: $af";
        push @filter_reasons, "TUM\_AF=" . $sample_info{$id}{'TUMOR'}{'AF'};
      }
    }

    my %information;

    foreach(@info_parts){

      my ($info_key, $info_value);

      if (/=/){
        ($info_key, $info_value) = $_ =~ /(.*)=(.*)/;
      }

      else {
        ($info_key) = $_ =~ /(.*)/;
        $info_value = "TRUE";
      }
      $information{$id}{$info_key} = $info_value;
    }

    $snvs{$id} = [ @fields[0..10], \@filter_reasons, \@samples ];

    $info{$id} = [ [@format], [%format_long], [%info_long], [@tumour_parts], [@normal_parts], [%information], [%sample_info] ];

    if ( scalar @filter_reasons == 0 ){
      $filtered_snvs{$.} = $_;
    }

  }

  return (\%snvs, \%info, \%filtered_snvs);
}
