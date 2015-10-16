#!/usr/bin/perl

use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
my $fasta;
my $use_imm;
my $mst;
my $mad;
my $fgenesh;
my $output;

my $blast0;
my $blast1;
my $file_gth0 ;
my $file_gth1 ;
my $file_gth2;
my $file_gth3 ;
my $file_gth4;
my $file_gth5;
my $output_html;
my $tool_directory;
GetOptions(
	'fasta=s'     => \$fasta,
	'use_imm=s'   => \$use_imm,
	'mst=s'       => \$mst,
	'mad=s'       => \$mad,
	'fgenesh=s'   => \$fgenesh,
	'output=s'	  => \$output,
#	'html=s'	  => \$output_html,
	'blast0=s'	  => \$blast0,
	'blast1=s'	  => \$blast1,
	'tool_directory=s' => \$tool_directory,
	'file_gth0=s' => \$file_gth0,    
	'file_gth1=s' => \$file_gth1, 
	'file_gth2=s' => \$file_gth2,
	'file_gth3=s' => \$file_gth3,
	'file_gth4=s' => \$file_gth4,
	'file_gth5=s' => \$file_gth5
);  
my  ($name,$path,$suffix) = fileparse($fasta);
my $seqobj = new Bio::SeqIO(
	-file => $fasta,
	-format => 'fasta'
)->next_seq;
my $file_par    = $use_imm eq "imm" ? $tool_directory."/models/eugene3.2_IMM.par" : $tool_directory."/models/eugene3.2.par";
my $file_matrice = $tool_directory."/models/matrice3.2.mat";
my $file_fasta  = $path.join(".",$seqobj->display_id,"fna");
system("cp $fasta $file_fasta");
my $file_spliceMAD = $path.join(".",$seqobj->display_id , "fna","spliceMAD");
my $file_spliceMSt = $path.join(".",$seqobj->display_id , "fna","spliceMSt");
system("cp $mst $file_spliceMSt");
system("cp $mad $file_spliceMAD");

if ($use_imm eq "imm") {
    my $file_fgenesh = $path.join(".",$seqobj->display_id , "fna","fgenesh.gff");
    system("cp $fgenesh $file_fgenesh"); 
    system($tool_directory."/eugene-3.2/bin/eugene-3.2 -s -A $file_par -m $file_matrice $file_fasta > $output 2>/dev/null");
  	print $tool_directory."/eugene-3.2/bin/eugene-3.2 -s -A $file_par -m $file_matrice $file_fasta > $output  \n";
}
else {

    my $file_fgenesh = $path.join(".",$seqobj->display_id , "fna","fgenesh.gff");
    my $file_blast0 = $path.join(".",$seqobj->display_id , "fna","blast0");
    my $file_blast1 = $path.join(".",$seqobj->display_id , "fna","blast1");
    my $file_est = $path.join(".",$seqobj->display_id , "fna","est");
    system("cp $blast0 $file_blast0");
    system("cp $blast1 $file_blast1");
    system("cat $file_gth0 $file_gth1 $file_gth2 $file_gth3 $file_gth4 $file_gth5 >> $file_est");
    system("cp $fgenesh $file_fgenesh");
    system($tool_directory."/eugene-3.2/bin/eugene-3.2 -s -A $file_par -m $file_matrice $file_fasta > $output 2>/dev/null");
  #  system($tool_directory."/eugene-3.2/bin/eugene-3.2 -ph -s -A $file_par -m $file_matrice $file_fasta > $output_html 2>/dev/null");
    
}
