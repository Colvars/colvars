#!/usr/bin/env perl
# Write a page that redirects links to the former single-page HTML manual
# (colvars-refman-<engine>.html#<LaTeX label>) to the multi-page manual
# built by LaTeXML (html-<engine>/<page>.html#<id>).
#
# Usage: make_redirect.pl <LaTeXML dbfile> <engine> <base URL> > page.html
#   dbfile:   written by latexmlc --dbfile; maps each LaTeX label to the page
#             and fragment where it ended up
#   engine:   gromacs, lammps, namd, tinkerhp or vmd
#   base URL: location of the multi-page manual relative to the redirect
#             page, e.g. "html-namd/" or "../master/html-namd/"
#
# Links with an unknown fragment (or none) go to the front page.

use strict;
use warnings;
use File::Basename qw(dirname);
use JSON::PP;
use LaTeXML::Util::ObjectDB;

my ($dbfile, $engine, $base) = @ARGV;
die "Usage: $0 <dbfile> <engine> <base URL>\n" unless defined $base;
die "$dbfile: not found\n" unless -f $dbfile;

my %names = (gromacs => 'GROMACS', lammps => 'LAMMPS', namd => 'NAMD',
  tinkerhp => 'Tinker-HP', vmd => 'VMD');
my $name = $names{$engine} // $engine;

my $db = LaTeXML::Util::ObjectDB->new(dbfile => $dbfile);
# Page locations are relative to LaTeXML's site directory (the current
# directory of latexmlc); make them relative to the manual's folder, which
# also holds the dbfile
my $dbdir = dirname($dbfile);
my %map;
foreach my $key ($db->getKeys) {
  next unless $key =~ /^LABEL:(.+)$/;
  my $label = $1;
  my $entry = $db->lookup($key)         or next;
  my $id    = $entry->getValue('id')    or next;
  my $obj   = $db->lookup("ID:$id")     or next;
  my $page  = $obj->getValue('location') or next;
  $page =~ s{^\Q$dbdir\E/}{} unless $dbdir eq '.';
  my $frag  = $obj->getValue('fragid');
  $map{$label} = defined $frag ? [$page, $frag] : [$page]; }
die "$dbfile: no labels found\n" unless %map;

my $json = JSON::PP->new->canonical->ascii->encode(\%map);
# Don't let the data close the <script> element
$json =~ s{</}{<\\/}g;
my $base_json = JSON::PP->new->allow_nonref->encode($base);
(my $base_html = $base) =~ s/&/&amp;/g;
$base_html =~ s/"/&quot;/g;

print <<"EOF";
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Collective Variables Module - Reference manual for $name</title>
<script>
// Former single-page manual: forward #<LaTeX label> to its new location
(function () {
  var base = $base_json;
  var map = $json;
  var label = "";
  try { label = decodeURIComponent(window.location.hash.slice(1)); } catch (e) {}
  var target = map[label];
  var url = base;
  if (target) {
    url += target[0];
    if (target.length > 1) url += "#" + encodeURI(target[1]);
  }
  window.location.replace(url);
})();
</script>
<noscript><meta http-equiv="refresh" content="0; url=$base_html"></noscript>
</head>
<body>
<p>The Colvars reference manual for $name has moved to
<a href="$base_html">$base_html</a>.</p>
</body>
</html>
EOF
