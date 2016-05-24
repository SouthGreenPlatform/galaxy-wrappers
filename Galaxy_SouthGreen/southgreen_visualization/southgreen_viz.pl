#!/usr/local/bin/perl
use strict;

my $viz = shift;
my $tree = shift;

my $align = "";
my $ideven = "";
my $expression = "";


if ($viz eq "genfamdashboard"){
    $align = shift;
    $ideven = shift;
    $expression = shift;
}


 my $out = shift;

#Fichier de sortie
open (OUT, ">$out") or die "cannot open out file : $! ";


my $html =" 

<?xml version=\"1.0\" encoding=\"utf-8\" ?>
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">
<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">
 <HTML>
<code>
 <HEAD>
  <TITLE></TITLE>
 </HEAD>
 <BODY>
 ";
 
#$html .="<FORM id=\"repartitor\" METHOD=\"POST\" ACTION=\"http://gohelle.cirad.fr/galaxy/static/repartitor/index.php?inputfile=";

if ($viz eq "intreegreat"){
    printf("test toto");
    $html .="<a href=\"http://cc2-web1.cirad.fr/repartitor/index.php?tree=";
    $html .="$tree";
    $html .="&visu=$viz\">";
}
else {
    $html .="<a href=\"http://cc2-web1.cirad.fr/repartitor/index.php?tree=";
    $html .="$tree";
    $html .="&visu=$viz";
    $html .="&exp=$expression";
    $html .="&ideven=$ideven";
    $html .="&align=$align\">";
}


#$html .="<input type=\"submit\" name =\"Submit\" value=\"submit\" >";

#</FORM> 
$html .="
Click here to display your data with the selected viewer.
</a> 
    </BODY>
</code>
</HTML>";


printf $html; 
print OUT $html;
#print $html;

close OUT;