#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my $tally_file;

GetOptions(
    'tally_file|f=s'     =>\$tally_file,
);


my %peptide_seq=();

#--------MAIN PROGRAM----------------------------------------
&check_incoming_data();
&tally_to_fasta($tally_file);
#-----------------------------------------------------------
 
 sub check_incoming_data{
  if(!$tally_file){
    print STDERR "Give me a tally file\n";
  &usage();
  exit(1)
 }
}

sub tally_to_fasta{
  my($f)=@_;
  my $corrected = "";
  my @line;
  my $cont =1;
  my $head = "";
  my $len = 0;
  
    open IDS, "$f" or die "Cannot open file  $f $!\n";
    
    #my $header=<IDS>;#delete the header
    foreach $_(<IDS>){
        chomp $_;
        (@line) = split(/\t/,$_);
        $corrected = &fasta_correction($line[0]);
        $len = scalar(@line)-1;
        $head = join(',', @line[1..$len]);
        print ">s$cont,$head\n$corrected\n";
        $cont++;
    }
    close IDS;
}

sub fasta_correction{
  my ($original) = @_;
  my $corrected = "";
  my $aux=0;
  
  my $len = length($original);
  for (my $i = 0; $i < $len; $i++){
    $corrected.=substr($original,$i,1);
    $aux++;
    if($aux==60){
      $corrected .= "\n";
      $aux = 0;
    }
  }
  
  return($corrected);
}



sub usage{
    print STDERR "\n This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR <<EOF;
NAME
    tally_to_fasta.pl   

USAGE
    tally_to_fasta.pl -f file.tally > file.fa

OPTIONS
    
Version 1.0, Copyright (C) 26/feb/18 EO

EOF
}

sub license{
    print STDERR <<EOF;

Copyright (C) 2018 EO
e-mail: obed\_eo\@hotmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
EOF
exit;
}