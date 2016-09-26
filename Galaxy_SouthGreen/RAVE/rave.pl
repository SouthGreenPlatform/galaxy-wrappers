#!/usr/bin/perl
 
use Getopt::Long;
use Pod::Usage;
use File::Temp qw/ tempfile tempdir /;
 
my $bfile;
my $maf;
my $threads;
my $recode;
my $extract;
my $keep;
my $output_format;
my $tool_directory;
my $out;
my $missing;
my $locus;
my $help;
my $variety_list;
my $location_list;
my $locus_list;

GetOptions(
    'bfile=s'             => \$bfile,
    'maf=f'               => \$maf,
    'threads=i'           => \$threads,
    'recode=s'            => \$recode,
    'extract=s{1,}'       => \$extract,
    'keep=s{1,}'          => \$keep,
    'variety_list=s'      => \$variety_list,
    'location_list=s'     => \$location_list,
    'locus_list=s'        => \$locus_list,
    'missing'             => \$missing,
    'locus=s'             => \$locus,
    'out=s'               => \$out,
    'output_format=s{1,}' => \$output_format,
    'tool_directory=s'    => \$tool_directory,
    'help|h|?'            => \$help
) ;
# --threads $threads
my $plink_cmd = "source $tool_directory/module_plink.sh;plink --bfile $bfile --maf $maf --silent";
if ($keep ne "" && $keep ne "None") {
    $plink_cmd .= " --keep $keep";
}
if ($extract ne "" && $extract ne "None") {
    $plink_cmd .= " --extract range $extract";
}
if ($variety_list ne "" && $variety_list ne "None") {
    my ($fh, $filename) = tempfile();
    my @variety = (split(/\_\_cn\_\_/,$variety_list));
    foreach my $variety (@variety) {
        print $fh join(" ",$variety,$variety) ,"\n";
    }
    #print $fh join("\n",@variety),"\n";
    $plink_cmd .= " --keep $filename";
}
if ($location_list ne ""  && $location_list ne "None") {
    my ($fh, $filename) = tempfile();
    my @location = (split(/\_\_cn\_\_/,$location_list));
    print $fh join("\n",@location),"\n";
    $plink_cmd .= " --extract range $filename";
}
if ($locus_list ne "" && $locus_list ne "None") {
    my $locus_file = $tool_directory ."/locus.txt";
    my ($fh, $filename) = tempfile();
    my @locus = (split(/\_\_cn\_\_/,$locus_list));
    my @id;
    foreach my $locus (@locus) {
        next if $locus eq "";
        my $command = `grep $locus $locus_file`;
        chomp($command); 
        $command =~s/\r//g;
        $command =~s/\n//g;
        push @id ,  $command;
    }
    close IN; 
    print $fh join("\n",@id),"\n"; 
    $plink_cmd .= " --extract range $filename";
}

if ($locus ne "" && $locus ne "None") {
    my $locus_file = $tool_directory ."/locus.txt";
    open(IN,$locus);
    my ($fh, $filename) = tempfile();
    my @id;
    while (<IN>) {
        chomp;
        next if $_ eq "";
        my $command = `grep $_ $locus_file`;
        chomp($command); 
        $command =~s/\r//g;
        $command =~s/\n//g;
        push @id ,  $command;
    }
    close IN; 
    print $fh join("\n",@id),"\n"; 
    $plink_cmd .= " --extract range $filename";
    
}

if ($output_format eq "vcf") {
    if ($bfile =~ /HDRA/) {
        $plink_cmd .= " --recode vcf-iid --out $out; sed -i 's/B\tA/C\tA/g' $out.vcf; sed -i 's/A\tB/A\tG/g' $out.vcf;sed -i 's/B\t\./C\t\./g' $out.vcf; sed -i 's/\.\tB/\.\tG/g' $out.vcf; mv $out.vcf $out";
    }
    else {
        $plink_cmd .= " --recode vcf-iid --out $out; mv $out.vcf $out";
    }
}
elsif ($output_format eq "bgz") {
    $plink_cmd .= " --recode vcf-iid bgz --out $out;mv $out.vcf.gz $out";   
}
elsif ($output_format eq "beagle") {
    $plink_cmd .= " --recode beagle --out $out;tar -Pczvf $out.tar.gz $out.chr-*; mv $out.tar.gz $out"; 
}
elsif ($output_format eq "fastphase") {
    $plink_cmd .= " --recode fastphase --out $out; tar -Pczvf $out.tar.gz $out.chr*recode.phase.inp; mv $out.tar.gz $out"; 
}
elsif ($output_format eq "structure") {
    $plink_cmd .= " --recode structure --out $out;mv $out.recode.strct_in $out";  
}

system($plink_cmd);
