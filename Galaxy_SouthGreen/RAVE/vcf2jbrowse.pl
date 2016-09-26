#!/usr/bin/perl
 
use Getopt::Long;
use lib;
use Pod::Usage; 
use Net::SCP;
use JSON;
binmode STDOUT, ":utf8";
use utf8;
#  
use File::Temp qw/ tempfile tempdir /;
 
my $vcffile; 
my $tool_directory; 
my $help; 
my $tbi;
my $bgzip; 
my $output_tbi;
my $host = "salanque.cirad.fr";
my $user = "galaxydev";
my $json = "template.json";
my $scp = Net::SCP->new($host);
$scp->login($user); 
GetOptions(
    'vcffile=s'       => \$vcffile, 
    'tool_directory=s'=> \$tool_directory, 
    'bgzip=s'         => \$bgzip,
    'tbi=s'           => \$tbi, 
    'help|h|?'        => \$help
) ; 
my $tabix_cmd;
if (-e $vcffile .".gz") {
    $tabix_cmd = "source " . $tool_directory. "/module_vcf2jbrowse.sh;   tabix -f -p vcf ". $vcffile .".gz";
}
else {
    $tabix_cmd = "source " . $tool_directory. "/module_vcf2jbrowse.sh; bgzip ". $vcffile ." ; tabix -p vcf ". $vcffile .".gz";
}
system("scp ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/trackList.json trackList.json");

my $json;
{
    local $/; #Enable 'slurp' mode
    open my $fh, "<", "trackList.json";
    $json = <$fh>;
    close $fh;
}
my $data = decode_json($json);

system("sed  s:FILE:tmp/galaxy$$.vcf.gz:g ".$tool_directory."/template.json > ".$tool_directory."/galaxy$$.json"); 

push @{$data->{'include'}} , "galaxy$$.json";
open my $fh, ">trackList.json";
print $fh to_json($data,{ utf8 => 1, pretty => 1 });
close $fh;
system($tabix_cmd);
my $file_gz = "galaxy$$.vcf.gz";
my $file_tbi = "galaxy$$.vcf.gz.tbi";

system("scp ". $vcffile.".gz.tbi ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp/" .$file_tbi); 
system("scp ". $vcffile.".gz ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp/". $file_gz);
system("scp ".$tool_directory."/galaxy$$.json ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/");
system("scp ".$tool_directory."/trackList.json ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/");

system("mv ". $vcffile.".gz.tbi " .$tbi);
system("mv ". $vcffile.".gz ". $bgzip);

 


