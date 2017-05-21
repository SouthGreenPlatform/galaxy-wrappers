#!/usr/bin/perl
use Getopt::Long;
use strict;
use CGI;
use Cwd;
my $dir = getcwd;

my $chromosome;
my $output;
my @tracks;
my @select;

GetOptions ("chromosome=s" => \$chromosome,    # numeric
              "tracks=s"   => \@tracks,      # string
              "select=s"  => \@select,
              "output=s" => \$output)   # flag
or die("Error in command line arguments\n");

open(O,">$output");
my $q = CGI->new;        

print O $q->header('text/html'); # create the HTTP header
print O $q->start_html(''); # start the HTML
system("cp $chromosome \$HOME/galaxy/static/style/blue/ideogram/");
my @part = split("/",$chromosome);

my $base_url = "http://galaxy.southgreen.fr/galaxy/";
if ($dir =~/galaxy_dev/){
	$base_url = "http://cc2-web1.cirad.fr/galaxydev/";
}


my $chromosomelength="$base_url/static/style/ideogram/".$part[$#part];
my $iframe = qq~
<form target="myIframe" action="http://salanque.cirad.fr/visu_genomeharvest/circosJS/demo/index.php" method="post" enctype="multipart/form-data">
<input type="hidden" name="chromosome" value="$chromosomelength" />
~;
foreach my $i (0 .. $#tracks){
	#print "$tracks[$i]";
	system("cp $tracks[$i] \$HOME/galaxy/static/style/blue/ideogram/");
	my @part = split("/",$tracks[$i]);
	my $track_file = "$base_url/static/style/ideogram/".$part[$#part];
	$iframe .= qq~
	<input type="hidden" name="select[]" id="annot" value="$select[$i]" />
	<input type="hidden" name="data[]" id="annot" value="$track_file" />
	~;
}    
$iframe .= qq~
<input type="submit" style="display: none;" id="load" value="reload">
</form>

<iframe src="" id="myIframe" name="myIframe" width=\"100%\" height=\"1700\" style='border:solid 0px black;'></iframe>
<script>
document.getElementById("load").click();
</script>
~;

print O $iframe;

close(O);
