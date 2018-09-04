#!/usr/bin/env perl

use strict;

use FileHandle;
use IO::Handle;
use IO::File;
#use DBI;
use IPC::Open2;
use Sys::Hostname;
#use Net::SSH::Perl;
use POSIX qw(setsid);
use POSIX 'WNOHANG';

require "/blue/web/services/feigservices.pl";

my $load=&getLoad();
my $cores=&getCores();
my $available= $cores-$load;
my $eightypercent=$cores*0.8;
if ($available>=$eightypercent){
 #CALL SCRIPT 
}

print "$available\n";

