#!/usr/bin/env perl

use strict;
use warnings;

my ($pdb, $ppdb, $harmlist, $out, $out2) = @ARGV;
select((select(OUT), $|=1) [0]);

open (PDB, "$ARGV[0]") || die("unable to open $pdb.");
my @raw_pdb=<PDB>;
close (PDB);

open (HARM, "$ARGV[1]") || die("unable to open $harmlist.");
my @raw_harm=<HARM>;
close (HARM);

my $xcoor={};
my $ycoor={};
my $zcoor={};
my $npdb = scalar @raw_pdb;

my @split_harm;
my $hinx1={};
my $hinx2={};
my $hdist0={};

my @split_lb;
my $lbinx1={};
my $lbinx2={};

for (my $i=0; $i < $npdb; $i++) {
  $xcoor->{$i+1}=substr($raw_pdb[$i],30,8)*10;
  $ycoor->{$i+1}=substr($raw_pdb[$i],38,8)*10;
  $zcoor->{$i+1}=substr($raw_pdb[$i],46,8)*10;
}

for (my $i=0; $i < scalar @raw_harm; $i++) {
  @split_harm=split('\s',$raw_harm[$i]);
  $hinx1->{$i}=int($split_harm[0]/15+0.94);
  $hinx2->{$i}=int($split_harm[1]/15+0.94);
  $hdist0->{$i}=$split_harm[3]*1;
}

my $minoffset;
my $minenergy=1000000;
my $atom;
my $xdist;
my $ydist;
my $zdist;
my $hdist;
my $harminx1;
my $harminx2;

my $oxcoor={};
my $oycoor={};
my $ozcoor={};

open(OUT,">>$ARGV[2]");
my $step=1;
#MJ
#now every one bead but only in very narrow region
#8840 8860

#MJ was
#for (my $k=0; $k < $npdb; $k=$k+10) {

for (my $k=8840; $k < 8860; $k=$k+1) {
  for (my $i=0; $i < $npdb; $i++) {
    $atom=($i+1+$k)%$npdb;
    if ($atom==0) { $atom=$npdb; }
    $oxcoor->{$atom}=$xcoor->{$i+1};
    $oycoor->{$atom}=$ycoor->{$i+1};
    $ozcoor->{$atom}=$zcoor->{$i+1};
  }
  my $harmener=0;
  for (my $i=0; $i < scalar @raw_harm; $i++) {
    $harminx1=$hinx1->{$i};
    $harminx2=$hinx2->{$i};
    if ($harminx1 < $npdb && $harminx2 < $npdb) {
      $xdist=($oxcoor->{$harminx1}-$oxcoor->{$harminx2})**2;
      $ydist=($oycoor->{$harminx1}-$oycoor->{$harminx2})**2;
      $zdist=($ozcoor->{$harminx1}-$ozcoor->{$harminx2})**2;
      $hdist=sqrt($xdist+$ydist+$zdist);
      $harmener+=5.0*($hdist-$hdist0->{$i})**2;
    }
  }
  my $energy=($harmener)/1000000;
  printf  OUT "Step: %d Shift: %d HarmE: %.3f Total: %.3f\n",$step,$k,$harmener,$energy;
  if ($energy < $minenergy) {
    $minenergy=$energy;
    $minoffset=$k;
  }
  $step++;
}

printf OUT "Minimum energy %.1f is found at offset %d\n",$minenergy,$minoffset;
close(OUT);

my $neworigin=$npdb-$minoffset+1;

my $inx=0;
open(OUT2,">>$ARGV[3]");
for (my $i=$neworigin; $i<=$npdb; $i++) {
  $inx++;
  if ($inx<10000) {
    printf OUT2 "ATOM %6s  CG  BD  A%4s    %8.3f%8.3f%8.3f  0.00  0.00%10s\n",$inx,$inx,$xcoor->{$i}/10,$ycoor->{$i}/10,$zcoor->{$i}/10,"";
  } elsif ($inx<100000) {
    printf OUT2 "ATOM %6s  CG  BD  A%5s   %8.3f%8.3f%8.3f  0.00  0.00%10s\n",$inx,$inx,$xcoor->{$i}/10,$ycoor->{$i}/10,$zcoor->{$i}/10,"";
  } else {
    printf OUT2 "ATOM %6s  CG  BD  A%6s  %8.3f%8.3f%8.3f  0.00  0.00%10s\n",$inx,$inx,$xcoor->{$i}/10,$ycoor->{$i}/10,$zcoor->{$i}/10,"";
  }
}

for (my $i=1; $i<$neworigin; $i++) {
  $minoffset++;
  if ($minoffset<10000) {
    printf OUT2 "ATOM %6s  CG  BD  A%4s    %8.3f%8.3f%8.3f  0.00  0.00%10s\n",$minoffset,$minoffset,$xcoor->{$i}/10,$ycoor->{$i}/10,$zcoor->{$i}/10,"";
  } elsif ($minoffset<100000) {
    printf OUT2 "ATOM %6s  CG  BD  A%5s   %8.3f%8.3f%8.3f  0.00  0.00%10s\n",$minoffset,$minoffset,$xcoor->{$i}/10,$ycoor->{$i}/10,$zcoor->{$i}/10,"";
  } else {
    printf OUT2 "ATOM %6s  CG  BD  A%6s  %8.3f%8.3f%8.3f  0.00  0.00%10s\n",$minoffset,$minoffset,$xcoor->{$i}/10,$ycoor->{$i}/10,$zcoor->{$i}/10,"";
  }
}
close(OUT2);

