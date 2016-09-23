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
my $output;
my $output_bgzip;
my $output_tbi;
my $organism;
my $host = "salanque.cirad.fr";
my $user = "galaxy";
my $scp = Net::SCP->new($host);
$scp->login($user);
my $ssh = Net::SSH::Perl->new($host);
$ssh->login($user) or die;
GetOptions(
    'vcffile=s'       => \$vcffile, 
    'tool_directory=s'=> \$tool_directory,
    'organism=s'      => \$organism, 
    'bgzip=s'         => \$bgzip,
    'tbi=s'           => \$tbi, 
    'output_bgzip=s'  => \$output_bgzip,
    'output_tbi=s'    => \$output_tbi, 
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


system($tabix_cmd);
my $file_gz = "/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp/galaxy$$.vcf.gz";
my $file_tbi = "/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp/galaxy$$.vcf.gz.tbi";


$scp->scp($vcffile.".gz ", $file_gz);
$scp->scp($vcffile.".gz.tbi ", $file_tbi);
system("mv ". $vcffile.".gz.tbi " .$output_tbi);
system("mv ". $vcffile.".gz ". $output_bgzip);



#my $cmd = "/share/apps/bin/fgenesh  /share/apps/share/fgenesh/$models -scip_prom -scip_term $name"; 
#$ssh->cmd($cmd_clean);

