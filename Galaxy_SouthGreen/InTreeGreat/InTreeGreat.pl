#!/usr/bin/perl

use strict;
use CGI;
my $input = $ARGV[0];
my $out = $ARGV[1];

my $newick = `cat $input`;

open(O,">$out");
#print O "Content-Type: text/html\n\n";
#print O "<html> <head>\n";
#print O "<title>Hello, world!</title>";
#print O "</head>\n";
#print O "<body>\n";
my $q = CGI->new;                        # create new CGI object
print O $q->header('text/html');                    # create the HTTP header
print O $q->start_html(''); # start the HTML

my $iframe = qq~
                <form target="myIframe" action="http://phylogeny.southgreen.fr/treedisplay_beta/index.php?" method="post">
                                <input type="hidden" name="hiddenfield" value="$newick" />
                                <input type="submit" value="Display tree">
                </form>

                <iframe src="" name="myIframe" width=\"1000\" height=\"700\" style='border:solid 1px black;'></iframe>
                ~;
print O $iframe;

#print O $q->end_html; 

close(O);
