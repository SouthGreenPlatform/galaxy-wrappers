#!/usr/bin/perl

#this program works like long_matrix.pl
#PLEASE specify the file from which the data is supposed to be read

$pseudo_count=0.002;

#transforms amino acids and positions into consecutive numbers for an array
sub transform_from_ascii
{
    if($_[0] eq "L")
    {
	return 0;
    }
    if($_[0] eq "I")
    {
	return 1;
    }
    if($_[0] eq "V")
    {
	return 2;
    }
    if($_[0] eq "M")
    {
	return 3;
    }
    if($_[0] eq "F")
    {
	return 4;
    }
    if($_[0] eq "Y")
    {
	return 5;
    }
    if($_[0] eq "G")
    {
	return 6;
    }
    if($_[0] eq "A")
    {
	return 7;
    }
    if($_[0] eq "K")
    {
	return 8;
    }
    if($_[0] eq "R")
    {
	return 9;
    }
    if($_[0] eq "H")
    {
	return 10;
    }
    if($_[0] eq "E")
    {
	return 11;
    }
    if($_[0] eq "D")
    {
	return 12;
    }
    if($_[0] eq "Q")
    {
	return 13;
    }
    if($_[0] eq "N")
    {
	return 14;
    }
    if($_[0] eq "S")
    {
	return 15;
    }
    if($_[0] eq "T")
    {
	return 16;
    }
    if($_[0] eq "C")
    {
	return 17;
    }
    if($_[0] eq "W")
    {
	return 18;
    }
    if($_[0] eq "P")
    {
	return 19;
    }
    if($_[0] eq "a")
    {
	return 0;
    }
    if($_[0] eq "b")
    {
	return 1;
    }
    if($_[0] eq "c")
    {
	return 2;
    }
    if($_[0] eq "d")
    {
	return 3;
    }
    if($_[0] eq "e")
    {
	return 4;
    }
    if($_[0] eq "f")
    {
	return 5;
    }
    if($_[0] eq "g")
    {
	return 6;
    }

}

#background frequencies calculated by "new_coils_profile1.pl"
#(calculated by new_coils_profile1.pl)
$frequency[0]=0.0970546802847317;
$frequency[1]=0.0580518405900064;
$frequency[2]=0.065116630369819;
$frequency[3]=0.0237161865636243;
$frequency[4]=0.0404713617792143;
$frequency[5]=0.0311864754477209;
$frequency[6]=0.0694457968900427;
$frequency[7]=0.0778357174549213;
$frequency[8]=0.0558291289319545;
$frequency[9]=0.0538969848951278;
$frequency[10]=0.0229431029575281;
$frequency[11]=0.0617118018338635;
$frequency[12]=0.051038538214323;
$frequency[13]=0.0397168458723436;
$frequency[14]=0.0436056156039428;
$frequency[15]=0.0719676474135461;
$frequency[16]=0.0569088459698085;
$frequency[17]=0.015900941173675;
$frequency[18]=0.0134111356930123;
$frequency[19]=0.050190722060794;

open(LONG,"temp_files/2stranded");

$residue_counter=0;

while($line=<LONG>)
{
    #get the register
    if($line=~/register:(\w+)/)
    {
	#put the register into an array
	@register=split //, $1;
    }
    
    #get the sequence
    if($line=~/sequence:(\w+)/)
    {
	$residue_counter+=length $1;
	#put the sequence into an array
	@sequence=split //, $1;

	#counts amino acids for each specific position (as in new_coils_profile3.pl)
	$counter=0;
	while(defined($sequence[$counter]))
	{
	    $AA=transform_from_ascii($sequence[$counter]);
	    $position=transform_from_ascii($register[$counter]);
	    if(defined($count[$AA][$position]))
	    {
		$count[$AA][$position]++;
	    }
	    else
	    {	
		$count[$AA][$position]=1;
	    }
	    $counter++;
	}
    }
}

#this printout gives the number of residues that are used for the calculation of the matrix
#print $residue_counter, "\n";


###calculation of the relative occurrences

#this for loop calculates $countpos
#it stands for the number of any amino acids at a particular position
#it is only a help for the real calculation in the next for loop
for($pos=0;$pos<7;$pos++)
{
    $count_pos[$pos]=0;
    for($AA=0;$AA<20;$AA++)
    {
	$count_pos[$pos]+=$count[$AA][$pos];
    }
}

#this is the actual calculation of the relative occurrences
#the relative occurrence of a particular amino acid at a particular position is calculated as
#(number of a particular amino acid at a particular position / number any amino acid at a particular position) / background frequency
#the results are given in a format the corresponds to "new.mat"

print "%\n";
print "% New matrix\n";
print "%\n";
print "% weighted\n";
print "w  14 1.89 0.30 1.04 0.27 20\n";
print "w  21 1.79 0.24 0.92 0.22 25\n";
print "w  28 1.74 0.20 0.86 0.18 30\n";
print "uw 14 1.82 0.28 0.95 0.26 20\n";
print "uw 21 1.74 0.23 0.86 0.21 25\n";
print "uw 28 1.69 0.18 0.80 0.18 30\n";
print "%\n";
print "%   a     b     c     d     e     f     g\n";

$AA=0;
print "L";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=1;
print "I";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=2;
print "V";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=3;
print "M";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=4;
print "F";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=5;
print "Y";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=6;
print "G";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=7;
print "A";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=8;
print "K";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=9;
print "R";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=10;
print "H";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=11;
print "E";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=12;
print "D";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=13;
print "Q";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=14;
print "N";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=15;
print "S";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=16;
print "T";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=17;
print "C";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=18;
print "W";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
$AA=19;
print "P";
for($pos=0;$pos<7;$pos++)
{
    $position_frequency[$AA][$pos]=($count[$AA][$pos]/$count_pos[$pos])/$frequency[$AA];
    printf (" %2.3f", $position_frequency[$AA][$pos]+$pseudo_count);
}
print "\n";
		
print "%\n";

close(LONG);
