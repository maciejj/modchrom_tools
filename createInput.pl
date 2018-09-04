#!/usr/bin/env perl

use strict;
use URI::Escape;
use IO::Handle;
use Email::Valid;
#use DBI;

#require "/blue/web/services/feigservices.pl"

#my $args= shift @ARGV;
#my @farg=split(/ /,$args);

#my $service=$farg[0];
#my $remotehost=$farg[1];
#my $encdata=$farg[2];
#my $data=uri_unescape($encdata);

#my $parameter={};

#my @fdata=split(/:::/,$data);
#foreach my $f (@fdata) {
 # my @t=split(/=:=/,$f);
 # $parameter->{$t[0]}=$t[1];
#}
my $service = "chromosome";
if ($service eq "chromosome"){
 my $copies=1;
 my $branchlenmin=60;
 my $branchlenmax=60;
 my $domainmin=60;
 my $domainmax=60;
 my $restraintsmin=600;
 my $restraintsmax=600;
 my $branchingpointsmin=300;
 my $branchingpointsmax=300;

my $FH;
 for(my $i=1; $i<=$copies;$i++){
  for (my $j=$branchlenmin; $j<=$branchlenmax; $j+=10){
   for (my $k=$domainmin; $k<=$domainmax; $k+=20){
    for (my $l=$restraintsmin; $l<=$restraintsmax; $l+=100){
     for (my $m=$branchingpointsmin; $m<=$branchingpointsmax; $m+=50){
        my $filename= "input/$i-$j-$k-$l-$m.inp";
	open $FH,'>',$filename;
	
	 print $FH "branchingpoints: ";
	 print $FH $m;
         print $FH "\nrestraints: ";
         print $FH $l;
         print $FH "\ndomain: ";
         print $FH $k;
         print $FH "\nbranchlength: ";
         print $FH $j;
         print $FH "\ncopy: ";
         print $FH $i; 	 
     }
    }
   } 
  }
 }
}
#close $FH;


