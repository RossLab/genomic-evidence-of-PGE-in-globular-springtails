#! /usr/bin/perl

# there are million ways to do this, I took this perl solution from
# https://www.biostars.org/p/172179/#172295

my $window = shift @ARGV
  or die "Specify a window size!";

while (<>) {
  my $pos = (split)[1];

  $cov += (split)[2];

  unless ($pos % $window) {
    printf "%s\t%s-%s\t%s\n", (split)[0], ($pos - $window), $pos, ($cov / $window);
    $cov = 0;
  }
}
