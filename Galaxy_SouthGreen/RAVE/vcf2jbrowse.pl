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


system($tabix_cmd);
my $file_gz = "galaxy$$.vcf.gz";
my $file_tbi = "galaxy$$.vcf.gz.tbi";

system("scp ". $vcffile.".gz.tbi ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp" .$file_tbi);
system("scp ". $vcffile.".gz ". $user ."@". $host.":/opt/projects/jbrowse.southgreen.fr/prod/jbrowse/oryza_sativa_japonica_v7/tmp". $file_gz);
 
system("mv ". $vcffile.".gz.tbi " .$tbi);
system("mv ". $vcffile.".gz ". $bgzip);



#my $cmd = "/share/apps/bin/fgenesh  /share/apps/share/fgenesh/$models -scip_prom -scip_term $name"; 
#$ssh->cmd($cmd_clean);

