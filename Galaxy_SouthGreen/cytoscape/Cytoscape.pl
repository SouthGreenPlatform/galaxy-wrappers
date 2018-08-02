#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

use Cwd;
my $dir = getcwd;


my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <input>
    -h, --html          <html_output>
~;
$usage .= "\n";

my ($infile,$htmlout);


GetOptions(
	"input=s"    => \$infile,
	"html=s"     => \$htmlout,
);


die $usage
  if ( !$infile);


my $datain = "";
open(I,$infile);
while(<I>){
	my $line = $_;
	$datain.=$line;
}
close(I);
# remove brackets at the beginning and end of JSON
my $new_datain = substr($datain,1,length($datain)-2);
$datain = $new_datain;
	
my @colors = ("#ed292a","#ed292a","#82ABA0","#2255a6","#6ebe43","#e76599","#662e91","#c180ff","#ea8b2f","#fff100","#666666","#01ffff","#bfbfbf","#2ac966","#666666");


my $pie_block = "";
for (my $i = 0; $i <= scalar @colors; $i++){
        $pie_block .= "'pie-$i-background-color': '$colors[$i]',\n";
        $pie_block .= "'pie-$i-background-size': 'mapData(group$i, 0, 10, 0, 100)',\n";
}

open(HTML_CYTOSCAPE,">$htmlout");
                        my $html = qq~<!DOCTYPE html>
<html><head>
<meta http-equiv="content-type" content="text/html; charset=UTF-8">
<link href="http://sniplay.southgreen.fr/cytoscape/Pie_style/style.css" rel="stylesheet">
<meta charset="utf-8">
<meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">
<title>Pie style</title>
<script src="http://sniplay.southgreen.fr/cytoscape/Pie_style/jquery.js"></script>
<script src="http://sniplay.southgreen.fr/cytoscape/Pie_style/cytoscape.js"></script>
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
$datain
  ,
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

