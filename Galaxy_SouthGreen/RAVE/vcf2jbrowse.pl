#!/usr/bin/perl
 
use Getopt::Long;
use Pod::Usage;
use Net::SSH::Perl;
use Net::SCP;
use File::Temp qw/ tempfile tempdir /;
 
my $vcffile; 
my $tool_directory; 
my $help; 
my $tbi;
my $bgzip; 
my $output_tbi;
my $host = "salanque.cirad.fr";
my $user = "galaxy";
my $json = "template.json";
my $scp = Net::SCP->new($host);
$scp->login($user);
#my $ssh = Net::SSH::Perl->new($host);
#$ssh->login($user) or die;
GetOptions(
    'vcffile=s'       => \$vcffile, 
    'tool_directory=s'=> \$tool_directory, 
    'bgzip=s'         => \$bgzip,
    'tbi=s'           => \$tbi, 
    'help|h|?'        => \$help
) ;
# --threads $threads
my $tabix_cmd;
if (-e $vcffile .".gz") {
    $tabix_cmd = "source " . $tool_directory. "/module_vcf2jbrowse.sh;   tabix -f -p vcf ". $vcffile .".gz";
}
else {
    $tabix_cmd = "source " . $tool_directory. "/module_vcf2jbrowse.sh; bgzip ". $vcffile ." ; tabix -p vcf ". $vcffile .".gz";
}

system("sed  s:FILE:tmp/galaxy$$.vcf.gz: > variants.json")
system($tabix_cmd);
my $file_gz = "galaxy$$.vcf.gz";
my $file_tbi = "galaxy$$.vcf.gz.tbi";

system("scp ". $vcffile.".gz.tbi ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp/" .$file_tbi); 
system("scp ". $vcffile.".gz ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp/". $file_gz);
system("scp variants.json ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/");

system("mv ". $vcffile.".gz.tbi " .$tbi);
system("mv ". $vcffile.".gz ". $bgzip);

 


