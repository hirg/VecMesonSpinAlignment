#!/usr/bin/perl -w 
use DBI;
use strict;
use Time::Local;

my $prod = (defined($ARGV[0])) ? $ARGV[0] : "not set!";
my $year = (defined($ARGV[1])) ? $ARGV[1] : "not set!";
#print "production: $prod\n";
#print "year: $year\n";
if ( ($prod =~ "not set!") || ($year =~ "not set!") ) {
  print "production flag or year not set!\n";
  die "usage: perl ...pl AuAu200_production 2011\n";
}

my $database;
if ( $year =~ 2010 ) { 
  $database = "dbbak.starp.bnl.gov:3409";
}
elsif ( $year =~ 2011 ) { 
  $database = "dbbak.starp.bnl.gov:3410";
}
else {
  die "year $year not implemented!"
}

#print "$database\n";

my $dbUser = '';
my $dbh_rts = DBI->connect("DBI:mysql:Conditions_rts:${database}",$dbUser,"")
  || die "cannot connect to server $DBI::errstr\n";

my $dbh_RunLog = DBI->connect("DBI:mysql:RunLog:${database}",$dbUser,"")
  || die "cannot connect to server $DBI::errstr\n";

my $dbh_RDO = DBI->connect("DBI:mysql:RunLog_onl:dbx.star.bnl.gov:3316",$dbUser,"");

my $query_base = "select run_number,base_trg_rate from run where TRG_SETUP_name='${prod}'";
my $sth_base = $dbh_rts->prepare($query_base);
my %base = ();
$sth_base->execute();
while (my $stuff = $sth_base->fetchrow_arrayref) {
  my $run = $stuff->[0];
  my $base = $stuff->[1];
  $base{$run} = $base;
}
my $query_ev = "select runNumber,numberOfEvents from daqSummary";
my $sth_ev = $dbh_RunLog->prepare($query_ev);
my %nev = ();
$sth_ev->execute();
while (my $stuff = $sth_ev->fetchrow_arrayref) {
  my $run = $stuff->[0];
  my $nev = $stuff->[1];
  $nev{$run} = $nev;
}
my @trignames = ();
if ($prod =~ /200/) {
  if ($year =~ 2010) {
    @trignames = qw(Central vpd-mb NPE_18 NPE_15 NPE_11 NPEHT_25);
  }
  if ($year =~ 2011) {
    @trignames = qw(vpd-zdc-mb-protected NPE_18 NPE_15 NPE_11 NPE_18_ftp);
  }
}
if ($prod =~ /62/) {
  @trignames = qw(min-bias min-bias-slow ht_11_mb central);
}
if ($prod =~ /39/) {
  @trignames = qw(mb mbslow ht-11);
}
if ($year =~ 2010) {
  if ($prod =~ /11/) {
    @trignames = qw(mb);
  }
  if ($prod =~ /7/) {
    @trignames = qw(mb);
  }
}
if ($year =~ 2011) {
  if ($prod =~ /27/) {
    @trignames = qw(mb1-fast vpd-mon-tac bbc-small-mon-narrow zdc-mon-tac);
  }
  if ($prod =~ /19/) {
    @trignames = qw(mb1-fast HLT-tracks vpd-mon bbc-small-mon-narrow zdc-mon-tac);
  }
}
#print @trignames;
#print "\n";
#print "---------------------\n";

my %trigs = ();
for my $trigname (@trignames) {
  my $query = "select a.runNumber,a.numberOfEvents from l0TriggerSet a,runDescriptor b,runStatus c, detectorSet d where a.runNumber = b.runNumber and b.runNumber = c.runNumber and c.runNumber = d.runNumber and a.name='$trigname' and b.glbSetupName='$prod' and c.rtsStatus=0 and d.detectorId=20";
  #  c.shiftLeaderStatus<1 and  and a.offlineTriggerID>1000 
  my $sth = $dbh_RunLog->prepare($query);
  $sth->execute();
  while (my $stuff = $sth->fetchrow_arrayref) {
    my $run = $stuff->[0];
    my $nev = $stuff->[1];
    $trigs{$run}{$trigname} = $nev;
  }
}

#print "Run\t";
#for my $name (@trignames) {
#  print "$name\t";
#}
#print "Total\t";
#print "Base\t";
#print "Nbadrdo\t";
#print "\n";
for my $run (sort {$a<=>$b} keys %trigs) {
  print $run,"\t";
  for my $trigname (@trignames) {
    if (defined($trigs{$run}{$trigname})) {
      print $trigs{$run}{$trigname},"\t"
    }
    else {
      print "0\t";
    }
  }
  if (defined($nev{$run})) {
    print $nev{$run},"\t";
  }
  else {
    print "0\t";
  }
  print "\n";
}
