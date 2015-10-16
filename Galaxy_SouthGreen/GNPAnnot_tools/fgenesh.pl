#!/usr/bin/perl

use lib '/usr/local/bioinfo/galaxy/galaxy_dist/tools/prediction/lib';
use lib '/usr/local/bioinfo/galaxy/galaxy_dist/tools/prediction/lib/lib/perl5/site_perl/5.8.8/';
use lib '/usr/local/bioinfo/galaxy/galaxy_dist/tools/prediction/lib/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi';
use HTTP::Request::Common qw(POST);
use HTTP::Request::Common qw(GET);
use Bio::Tools::Fgenesh;
use Bio::Tools::Prediction::Gene;
use Sys::Hostname;
use LWP::UserAgent; 
use Bio::SeqIO;
use File::Basename;
use Net::SSH::Perl;
use Net::SCP;

my $file   = shift;
my $models = shift;
# Monocots
# Botrytis         
# MEDICAGO 
# Dicots           
# Magnaporthe      
# Sclerotinia
# Leptosphaeria 
my  $file_fgenesh_gff = &get_result_via_ssh($file,$models);
print $file_fgenesh_gff , "\n";


sub get_result_via_ssh {
    my ($file,$models) = @_;
    my $file_raw_fg = $file .".raw.fg";
    my $file_fgenesh_gff = $file .".fgenesh.gff";
    my $host = "gpi.versailles.inra.fr";
    my $user = "gdroc";
    my $hostname = hostname;
    my $id_dsa = "/usr/local/bioinfo/galaxy/.ssh/id_dsa";
    my $id_dsa_pub = "/usr/local/bioinfo/galaxy/.ssh/id_dsa.pub";
    my $current_id_dsa = join(".",$id_dsa,$hostname);
    my $current_id_dsa_pub = join(".",$id_dsa,$hostname,"pub");
    system("cp $current_id_dsa $id_dsa");
    system("cp $current_id_dsa_pub $id_dsa_pub");
    my $remote_file = basename($file);
    my $scp = Net::SCP->new($host);
    $scp->login($user);
    my $ssh = Net::SSH::Perl->new($host);
    $ssh->login($user) or die;
    my $cmd = "fgenesh2.4  /usr/local/share/fgenesh/$models -scip_prom -scip_term $remote_file";
    print $cmd ,"\n";
    my $cmd_clean = "rm $remote_file"; 
    $scp->put($file)  or die $scp->{errstr};
    my($stdout, $stderr, $exit) = $ssh->cmd($cmd);
    $ssh->cmd($cmd_clean);
    open(OUT,">",$file_raw_fg);
    print OUT $stdout ."\n";
    close OUT;
    #$file_fgenesh_gff = &reformat_fgenesh_gff($file_raw_fg,$file_fgenesh_gff);
    return $file_raw_fg;#fgenesh_gff;
}



sub reformat_fgenesh_gff {
    my ($file_raw_fg,$file_fgenesh_gff) = @_;
    open(IN,$file_raw_fg);
    open(OUTPUT,">".$file_fgenesh_gff);
    while (my $line = <IN>) {
        my @a_line = split (/\s+/, $line);
        if ($#a_line == 12  &&  $a_line[4]=~ /^CDS/) {
            my $seqname = $a_line[1];
            my $strand  = $a_line[2];
            my $feature = $a_line[4];
            my $score   = $a_line[8];
            my $start   = $a_line[5];
            my $end     = $a_line[7];
            $feature = "E.Init" if $feature eq "CDSf";
            $feature = "E.Intr" if $feature eq "CDSi";
            $feature = "E.Term" if $feature eq "CDSl";
            $feature = "E.Sngl" if $feature eq "CDSo";
            print OUTPUT join("\t",$seqname,"fgenesh",$feature,$start,$end,$score,$strand,"."),"\n";
        }
    }
    close OUTPUT;
    close INPUT;
    return $file_fgenesh_gff;
}


1;
