#!/usr/bin/perl

use strict;
use CGI;
use Cwd;
my $dir = getcwd;

my $chromosomelength = $ARGV[0];
my $annotation = $ARGV[1];
my $ploidy = $ARGV[2];
my $out = $ARGV[3];

system("cp $chromosomelength \$HOME/galaxy/static/style/blue/ideogram/");
system("cp $annotation \$HOME/galaxy/static/style/blue/ideogram/");
my $base_url = "http://galaxy.southgreen.fr/galaxy";
if ($dir =~/galaxy_dev/){
        $base_url = "http://cc2-web1.cirad.fr/galaxydev";
}
my @part = split("/",$chromosomelength);
$chromosomelength="$base_url/static/style/ideogram/".$part[$#part];
@part = split("/",$annotation);
$annotation="$base_url/static/style/ideogram/".$part[$#part];


open(O,">$out");
my $q = CGI->new;                        # create new CGI object
print O $q->header('text/html');                    # create the HTTP header
print O $q->start_html(''); # start the HTML

my $iframe = qq~

                <form target="myIframe" action="http://cassava-genome.southgreen.fr/visu_genomeharvest/ideogram/newindex.php" method="post">
                                <input type="hidden" name="data" id="data" value="$chromosomelength" />
				<input type="hidden" name="annot" id="annot" value="$annotation" />
				<input type="hidden" name="select" id="selectorpost" value="$ploidy" />
                                <input type="submit" style="display: none;" id="load" value="reload">
                </form>
           
                <iframe src="" id="myIframe" name="myIframe" width=\"100%\" height=\"1700\" style='border:solid 0px black;'></iframe>
                     <script>
                document.getElementById("load").click();
  				</script>
                ~;
print O $iframe;

close(O);
