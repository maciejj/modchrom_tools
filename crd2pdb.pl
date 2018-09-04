#!/usr/bin/env perl

while (<>) {
  chop;
  if ( /^ *[0-9]/ ) {
    if ($natom==0) { 
      $natom=$_;
    } else {
      s/ +/ /g;
      @f=split();
      if ($f[1]<10000) {
        printf "ATOM %6d %-4s %-4s %4d    %8.3f%8.3f%8.3f %5.2f %5.2f      %-4s\n",
        $f[0],$f[3],$f[2],$f[1],$f[4],$f[5],$f[6],0.0,0.0,$f[7];
      } else {
        printf "ATOM %6d %-4s %-4s%5d    %8.3f%8.3f%8.3f %5.2f %5.2f      %-4s\n",
        $f[0],$f[3],$f[2],$f[1],$f[4],$f[5],$f[6],0.0,0.0,$f[7];
      }

    }
  }
}
