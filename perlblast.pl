#!/usr/bin/env perl
use 5.010;
use warnings;
use strict;

# $/= "";

# read in from a file a query string Q
print "Input filename: ";
open(IN, <>);
my $Q = "temp";
while(my $line = <IN>)
{
  $Q = $line;
  chomp $Q;
  last; # read only the first string
}
my $Qcopy = $Q;
my @Qchars = split("", $Q);
close IN;


# for k=4, make a hash table to find all locations of each different k-mer in Q
my $k = 3;
my %Qhash = ();
my $i = 1;
while (length($Q) >= $k)
{
  $Q =~ m/(.{$k})/;
  if (!defined $Qhash{$1})
  {
    $Qhash{$1} = [$i-1];
  }
  else
  {
    push (@{$Qhash{$1}}, $i-1)
  }
  $i++;
  $Q = substr($Q, 1, length($Q)-1);
}

print "Input threshold t: ";
my $t = <>;

# read in one string at a time from perlblastdata.txt 
# when a string S is read in, scan through its 4-mers, making a hash table for it
open(IN, "perlblastdata.txt");
while(my $S = <IN>)
{
  chomp $S;
  my $Scopy = $S;
  my @Schars = split("", $S);
  my %Shash = ();
  my $i = 1;
  while(length($S) >= $k)
  {
    $S =~ m/(.{$k})/;
    if (!defined $Shash{$1})
    {
      $Shash{$1} = [$i-1];
    }
    else
    {
      push (@{$Shash{$1}}, $i-1)
    }
    $i++;
    $S = substr($S, 1, length($S)-1);
  }



  # when a kmer is in both Q and S, extract the value from Qhash
  # extend the kmer in Q and S as long as there are matching characters
  # if L>10, print a message
  my %stringhash = ();
  foreach my $Qhashkey (keys(%Qhash)) # for every unique kmer in Q
  {
    if (defined $Shash{$Qhashkey}) # if S also has that kmer
    {
      foreach my $Qkmer (@{$Qhash{$Qhashkey}}) # for every time the kmer appears in Q
      {
        foreach my $Skmer (@{$Shash{$Qhashkey}}) # for every time the kmer appears in S
        {
          my $Qkmerpos = $Qkmer;
          my $Skmerpos = $Skmer;
          my $kmerstartpos = $Skmerpos;
          my $kmerendpos = $Skmerpos + $k - 1;
          my $L = $k; 

          my $counter = 1;
          while($Qchars[$Qkmerpos-$counter] eq $Schars[$Skmerpos-$counter])
          {
            $L++;
            $counter++;
            $kmerstartpos--;
          }
          $counter = 0; 
          while($Qchars[$Qkmerpos + $k + $counter] eq $Schars[$Skmerpos + $k + $counter])
          {
            $L++;
            $counter++;
            $kmerendpos++;
          }
          if($L > $t)
          {
            if (!defined $stringhash{$kmerstartpos})
            {
              my $substring = "";
              for (my $i = $kmerstartpos; $i <= $kmerendpos; $i++) 
              {
                $substring = $substring . $Schars[$i];
              }
              $stringhash{$kmerstartpos} = $substring;
              print "$substring\n";
            }
          }
        }
      }
    }
  }
}
close IN;