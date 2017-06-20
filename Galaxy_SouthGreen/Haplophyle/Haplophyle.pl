#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

my $HAPLOPHYLE_EXE = "java -Xmx2048m -jar /usr/local/bioinfo/sniplay/Haplophyle/NetworkCreator_fat.jar";
my $NEATO_EXE = "neato";
my $CONVERT_EXE = "convert";
my $RSCRIPT_EXE = "/usr/local/bioinfo/R/default/bin/Rscript";

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <input>
    -o, --output        <output>
    -d, --dotfile       <dotfile>
    -h, --html          <html_output>
~;
$usage .= "\n";

my ($infile,$output,$outfile,$htmlout);


GetOptions(
	"input=s"    => \$infile,
	"output=s"   => \$output,
	"dot=s"      => \$outfile,
	"html=s"     => \$htmlout
);


die $usage
  if ( !$infile);


  
	
my $out_png = "network.png";

my $command = "$HAPLOPHYLE_EXE -in $infile -out $outfile";
system($command);

	
	
open(OUTFILE,"$outfile");
open(OUTFILE2,">$output");
print OUTFILE2 "var nodes = [\n";
open(INFILE,"$infile");
while(<INFILE>)
{
	if (/>haplo(\d+)\|(\d+)/)
	{
		print OUTFILE2 "{id: $1, label: 'Haplo$1',shape:'image',radius:'$2'}\n";
	}
}
close(INFILE);
print OUTFILE2 "];\n";
print OUTFILE2 "var edges = [\n";
my $n = 0;
while(<OUTFILE>)
{
	my $line = $_;
	if (/^haplo(\d+) -- haplo(\d+)/)
	{
		my $from = $1;
		my $to = $2;
		
		print OUTFILE2 "{from: $from, to: $to,length:'7',style: 'line',color:'black',},\n";
	}	
}
close(OUTFILE);
print OUTFILE2 "];\n";
#close(OUTFILE2);


my @colors = ("#ed292a","#ed292a","#82ABA0","#2255a6","#6ebe43","#e76599","#662e91","#c180ff","#ea8b2f","#fff100","#666666","#01ffff","#bfbfbf","#2ac966","#666666");
my $pie_block = "";
my $nb_groups = 3;
for (my $i = 1; $i <= $nb_groups; $i++){
	$pie_block .= "'pie-$i-background-color': '$colors[$i]',\n";
	$pie_block .= "'pie-$i-background-size': 'mapData(group$i, 0, 10, 0, 100)',\n";
}
my $session = int(rand(100000));
my $home = `echo \$HOME`;
$home=~s/\n//g;$home=~s/\n//g;$home=~s/\r//g;
my $out_html = "$home/galaxy/static/style/blue/cytoscape/$session.cytoscape.htm";
print $out_html;
open(HTML_CYTOSCAPE,">$out_html");
                        my $html = qq~<!DOCTYPE html>
<html><head>
<meta http-equiv="content-type" content="text/html; charset=UTF-8">
<link href="Pie_style/style.css" rel="stylesheet">
<meta charset="utf-8">
<meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">
<title>Pie style</title>
<script src="Pie_style/jquery.js"></script>
<script src="Pie_style/cytoscape.js"></script>
<script type="text/javascript">
\$(function(){ // on dom ready

\$('#cy').cytoscape({

  style: cytoscape.stylesheet()
    .selector(':selected')
      .css({
        'background-color': 'black',
        'line-color': 'black',
        'opacity': 1
      })
    .selector('.faded')
      .css({
        'opacity': 0.25,
        'text-opacity': 0
      })
    .selector('edge')
      .css({
                'width': 1,
                'line-color': 'black',
      })
    .selector('node')
      .css({
        'width': 'mapData(width, 0, 10, 0, 100)',
        'height': 'mapData(width, 0, 10, 0, 100)',
        'content': 'data(id)',
        'pie-size': '98%',
        $pie_block
      }),
elements: {
    nodes: [
~;
my $done = 0;
                open(OUTFILE,"$outfile");
                while(<OUTFILE>){
                        if (/(^\w+)\s\[.*width=([\d\.]+),/){
                                my $node = $1;
                                my $size = $2 * 10;
                                #my $ref_hash = $hash{$gene}{$node};
                                my $ref_hash;
                                if ($ref_hash){
                                        my %hash2 = %$ref_hash;
                                        my $s = scalar keys(%hash2);
                                        $html.= "{ data: { id: '$node', width: $size";
                                        for (my $i = 1; $i <= $nb_groups; $i++){
                                                #my $ratio = $hash{$gene}{$node}{$i};
                                                #$html .= ", group$i: $ratio";
                                        }
                                        $html.= " } },\n";
                                }
                                else{
                                        $html.= "{ data: { id: '$node', width: $size} },\n";
                                }
                        }
                        if (/(\w+) -- (\w+)/){
                                if ($done == 0){
                                        $done = 1;
                                        $html .= "],\n";
                                        $html .= "edges: [\n";
                                }
                                $done = 1;
                                $html.= "{ data: { id: '$1$2', weight: 1, source: '$1', target: '$2'} },\n";
                        }
                }
                $html.= qq~
]
  },
layout: {
        name: 'cose',
    padding: 10
  },

  ready: function(){
    window.cy = this;
  }
});

});

</script>
</head>
<body>
<div id="cy">
</div>
~;
print HTML_CYTOSCAPE $html;
close(HTML_CYTOSCAPE);


use Cwd;
my $dir = getcwd;
my $base_url = "http://galaxy.southgreen.fr/galaxy/";
if ($dir =~/galaxy_dev/){
        $base_url = "http://cc2-web1.cirad.fr/galaxydev/";
}
#system("cp -rf $out_html $htmlout");
open(HTML,">$htmlout");
my $iframe = qq~
<a href="$base_url/static/style/cytoscape/$session.cytoscape.htm" target=_blank>Access to the Cytoscape visualisation of haplotype network</a>
~;
print HTML $iframe;
close(HTML);
