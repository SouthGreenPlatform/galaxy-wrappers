#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use lib ".";

my $usage = qq~Usage:$0 <args>
where <args> are:
    -c, --count           <SNP count SNiPloid output>
    -a, --annotation      <annotation file in GFF3>
    -o, --output_png      <output PNG file>
    -s, --scale           <scale. Default:100000>
    -t, --type_analysis   <type of analysis: polyploid_diploid or polyploid_polyploid. Default:polyploid_diploid>
    -m, --max_nb_chrom    <maximum number of chromomsome to display. Default:20>
    -n, --nb_min_snp      <minimal number of SNP to calculate ratio. Default:10>
    -d, --display_cat     <display ratio for each category instead of intra-polyploid and inter-diploid (yes/no). Default:no>
~;
$usage .= "\n";

my ($snp_count,$annotation,$output_png,$global_scale,$max_nb_chrom,$type_analysis);

my $global_scale = 10000;
my $max_nb_chrom = 20;
my $nb_min_snp   = 10;
my $display_cat = "no";


GetOptions(
	"count=s"           => \$snp_count,
	"annotation=s"      => \$annotation,
	"output_png=s"      => \$output_png,
	"scale=s"           => \$global_scale,
	"max_nb_chrom=s"    => \$max_nb_chrom,
	"type_analysis=s"   => \$type_analysis,
	"nb_min_snp=s"      => \$nb_min_snp,
	"display_cat=s"     => \$display_cat
);


die $usage
  if ( !$snp_count || !$annotation || !$output_png || !$type_analysis);
  
my %proportions_categories;
my %ratios; 
my %ratios_poly_diploid;
my %nb_snps;
open(my $COUNT,$snp_count);
<$COUNT>;
while(<$COUNT>)
{
	my $line =$_;
	chomp($line);
	my @infos = split(/\t/,$line);
	
	
	if ($type_analysis eq "polyploid_diploid")
	{
		my $gene = $infos[0];
		my $nb_snp = $infos[2];
		my $nb_1 = $infos[3];
		my $nb_2 = $infos[4];
		my $nb_3or4 = $infos[5];
		my $nb_3 = $infos[6];
		my $nb_4 = $infos[7];
		my $nb_5 = $infos[8];
		my $nb_other = $infos[9];
		my $nb_heterozygot_diploid = $infos[10];
		my $nb_snp_diploid = $infos[11];
		my $nb_snp_polyploid = $infos[12];
		
		$nb_snps{$gene} = $nb_snp;	
		my $sum = $nb_1 + $nb_2 + $nb_3or4 + $nb_5 + $nb_3 + $nb_4 + $nb_other + $nb_heterozygot_diploid;
		
		if ($nb_snp >= $nb_min_snp)
		{
			if ($nb_1)
			{
				$proportions_categories{$gene}{"1"} = $nb_1/$nb_snp;
			}
			if ($nb_2)
			{
				$proportions_categories{$gene}{"2"} = $nb_2/$nb_snp;
			}
			if ($nb_5)
			{
				$proportions_categories{$gene}{"5"} = $nb_5/$nb_snp;
			}
			if ($nb_3or4)
			{
				$proportions_categories{$gene}{"3or4"} = $nb_3or4/$nb_snp;
			}
		}

		my $ratio_g1 = $infos[13];
		my $ratio_g2 = $infos[14];
	
		if ($ratio_g1)
		{
			$ratios{$gene} = $ratio_g1;
		}
	}
	
	if ($type_analysis eq "polyploid_polyploid")
	{
		my $gene = $infos[0];
		my $nb_snp = $infos[2];
		my $nb_equal = $infos[3] + $infos[4];
		my $nb_diff = $infos[7];
		$nb_snps{$gene} = $nb_snp;
		if ($nb_snp >= $nb_min_snp)
		{
			if ($nb_equal)
			{
				$proportions_categories{$gene}{"equal"} = $nb_equal/$nb_snp;
			}
			if ($nb_diff)
			{
				$proportions_categories{$gene}{"difference"} = $nb_diff/$nb_snp;
			}
		}
	}
	
}
close($COUNT);


  

my $max_pos = 0;
my %chrom_sizes;
my $chrom_particule;
my %genes;
my %gene_positions;

open(my $ANNOT,$annotation);
while(<$ANNOT>)
{
	my $line =$_;
	chomp($line);
	if (!/^#/ && /gene/)
	{
		my @infos = split(/\t/,$line);
		my $chrom = $infos[0];
		if ($chrom =~/^(\w+_)(\d+)$/)
		{
			$chrom_particule = $1;
			$chrom = $2;
		}
		
		my $attributes = $infos[8];
		my $gene_name;
		if ($attributes =~/Name=([^;]+);/)
                {
                        $gene_name = $1;
                }
		if (!$gene_name && $attributes =~/ID=([^;]+);/)
                {
                        $gene_name = $1;
                }
		if ($gene_name =~/(.*)_G1/)
		{
			$gene_name = $1;
		}
		else
		{
			next;
		}
		if (not defined $nb_snps{$gene_name})
		{
			next;
		}
		
		my $start = $infos[3];
		my $end = $infos[4];
		my $pos = sprintf("%.0f", ($start + (($end - $start) / 2)));
		
		$end = $end / $global_scale;
		if ($chrom_sizes{$chrom})
		{
			if ($end > $chrom_sizes{$chrom})
			{
				$chrom_sizes{$chrom} = $end;
				if ($end > $max_pos)
				{
					$max_pos = $end;
				}
			}
		}
		else
		{
			$chrom_sizes{$chrom} = $end;
			if ($end > $max_pos)
			{
				$max_pos = $end;
			}
		}
		$genes{$gene_name} = "$chrom:$pos";
		$gene_positions{$chrom}{$pos}= $gene_name;
	}
}
close($ANNOT);



use GD;
use GD::Simple;
use GD::Image;




####################
# drawing
####################

my $scale = 800 / $max_pos;

my $margin_left = 80;
my $margin_right = 50;
my $margin_top = 50;
my $margin_legend = 100;
my $margin_bottom = 10;
my $margin_between_chromosomes = 25;
my $margin_between_section = 50;
my $chrom_width = 10;
my $gene_width = 1;

my $nb_group = 1;

my $width_of_picture = scalar keys(%gene_positions);
if (scalar keys(%gene_positions) > $max_nb_chrom)
{
	$width_of_picture = $max_nb_chrom;
}
                                        
my $diagram_img = GD::Simple->new(($margin_left + $margin_right + ($max_pos*$scale)),
                                         ($margin_top + ((($chrom_width * $nb_group) + ($margin_between_chromosomes * ($nb_group-1))) * $width_of_picture) + ($margin_between_section * $width_of_picture) + $margin_bottom + $margin_legend)
                                        );                                        

my $yellow =  $diagram_img->colorAllocate(247,254,46);
my $orange_light = $diagram_img->colorAllocate(250,204,46);
my $red_light = $diagram_img->colorAllocate(254,100,46);
my $red = $diagram_img->colorAllocate(254,46,46);
my $orange = $diagram_img->colorAllocate(254,154,46);


# draw chromosomes
my $num_chrom = 0;
my @sorted_chrom = sort {$a <=> $b} keys(%gene_positions);

my $nombre_genes = 0;
my $y_end;
foreach my $chrom(@sorted_chrom)
{
		if (!$chrom)
		{
			next;
		}
	        
	    if ($num_chrom > ($max_nb_chrom - 1))
	    {
	    	last;
	    }
        my $ref_hash = $gene_positions{$chrom};
        my %hash = %$ref_hash;

        my $section_size = $chrom_width + (($margin_between_chromosomes + $chrom_width) * ($nb_group - 1));

        # draw chromosome (X number of groups)

        $diagram_img->fgcolor('black');
        $diagram_img->bgcolor('white');
        $diagram_img->setThickness(1); 
        my $chrom_chain = $chrom_particule . $chrom;

        $diagram_img->rectangle( $margin_left,
                                        $margin_top + (($section_size + $margin_between_section) * $num_chrom),
                                        $margin_left + ($chrom_sizes{$chrom}*$scale),
                                        $margin_top + $chrom_width + (($section_size + $margin_between_section) * $num_chrom)
                                        );

        $diagram_img->fgcolor('black');
        $diagram_img->moveTo(5,$margin_top + $chrom_width + (($section_size + $margin_between_section) * $num_chrom) - 1);
        $y_end = $margin_top + $chrom_width + (($section_size + $margin_between_section) * ($num_chrom+1)) - 1;
        $diagram_img->fontsize(12);
        $diagram_img->font('Times');
        $diagram_img->string($chrom_particule . $chrom);
        
        
        my $previous_x_5;
        my $previous_x_1;
        my $previous_x_2;
        my $previous_x_3or4;
		my $previous_y_5;
		my $previous_y_1;
		my $previous_y_2;
		my $previous_y_3or4;
		
		my $previous_x_equal;
        my $previous_x_diff;
		my $previous_y_equal;
		my $previous_y_diff;
		
		my $previous_x_snp_diplo;
        my $previous_x_snp_poly;
		my $previous_y_snp_diplo;
		my $previous_y_snp_poly;

		my $previous_x_ratio_diplo_poly;
		my $previous_y_ratio_diplo_poly;
        
        # draw genes
        foreach my $pos(sort{$a <=> $b}keys(%hash))
        {
        	my $gene = $gene_positions{$chrom}{$pos};
        	if (not defined $nb_snps{$gene} or $nb_snps{$gene} < $nb_min_snp)
        	{
        		next;
        	}
        	
        	if ($type_analysis eq "polyploid_diploid")
			{
	        	#####################################################
	        	# draw ratio (subgenomic contribution)
	        	#####################################################
	        	my $color = "gray";
	        	if ($ratios{$gene})
	        	{
	        		my $ratio_g1 = $ratios{$gene};
	        		if ($ratio_g1 <= 30)
	        		{
	        			$color = $red;
	        		}
	        		elsif ($ratio_g1 > 30 && $ratio_g1 <= 40)
	        		{
	        			$color = $red_light;
	        		}
	        		elsif ($ratio_g1 > 40 && $ratio_g1 <= 60)
	        		{
	        			$color = $orange;
	        		}
	        		elsif ($ratio_g1 > 60 && $ratio_g1 <= 70)
	        		{
	        			$color = $orange_light;
	        		}
	        		elsif ($ratio_g1 > 70)
	        		{
	        			$color = $yellow;
	        		}
	        	}
	
	            $pos = $pos / $global_scale;
	            
	            $diagram_img->fgcolor($color);
	            $diagram_img->bgcolor($color);
	            $diagram_img->rectangle( $margin_left + ($pos*$scale) - ($gene_width / 2),
	                                                $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1,
	                                                $margin_left + ($pos*$scale) + ($gene_width / 2),
	                                                $margin_top + $chrom_width + (($section_size + $margin_between_section) * $num_chrom) - 1
	                                                );
	                                                
	                 
	                                           
	            #####################################################
	        	# draw SNP categories
	        	#####################################################
	        	
				my $proportion_5 = $proportions_categories{$gene}{"5"};
				my $proportion_1 = $proportions_categories{$gene}{"1"};
				my $proportion_2 = $proportions_categories{$gene}{"2"};
				my $proportion_3or4 = $proportions_categories{$gene}{"3or4"};
				my $ratio_poly_diplo = $ratios_poly_diploid{$gene};
		
				my $draw = 0;
				if (defined $previous_x_5)
				{
					$draw = 1;
				}
					
					
				#######################
				# SNP category 5
				#######################
				if ($draw)
				{
					$diagram_img->moveTo($previous_x_5,$previous_y_5);
				}    
				$previous_x_5 = $margin_left + ($pos*$scale) - 1;    
				$diagram_img->setThickness(2);                          
				$diagram_img->fgcolor("red");
				$diagram_img->bgcolor("red");
				$previous_y_5 = $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1 - ($proportion_5 * 20) - 7;
				if ($draw)
				{
			        $diagram_img->lineTo($previous_x_5,$previous_y_5);  
				}
					
				if ($display_cat eq "yes")
				{	
					#######################
					# SNP category 1
					#######################
					if ($draw)
					{
						$diagram_img->moveTo($previous_x_1,$previous_y_1);
					}    
					$previous_x_1 = $margin_left + ($pos*$scale) - 1;                                
			        $diagram_img->fgcolor("orange");
		            $diagram_img->bgcolor("orange");
			        $previous_y_1 = $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1 - ($proportion_1 * 20) - 7;
			        if ($draw)
					{
			        	$diagram_img->lineTo($previous_x_1,$previous_y_1); 
					}
					
					
					#######################
					# SNP category 2
					#######################
					if ($draw)
					{
						$diagram_img->moveTo($previous_x_2,$previous_y_2);
					}    
					$previous_x_2 = $margin_left + ($pos*$scale) - 1;                                
			        $diagram_img->fgcolor("purple");
		            $diagram_img->bgcolor("purple");
			        $previous_y_2 = $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1 - ($proportion_2 * 20) - 7;
			        if ($draw)
					{
			        	$diagram_img->lineTo($previous_x_2,$previous_y_2); 
					}
					
					
					#######################
					# SNP category 3 or 4
					#######################
					if ($draw)
					{
						$diagram_img->moveTo($previous_x_3or4,$previous_y_3or4);
					}    
					$previous_x_3or4 = $margin_left + ($pos*$scale) - 1;                                
			        $diagram_img->fgcolor("green");
		            $diagram_img->bgcolor("green");
			        $previous_y_3or4 = $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1 - ($proportion_3or4 * 20) - 7;
			        if ($draw)
					{
			        	$diagram_img->lineTo($previous_x_3or4,$previous_y_3or4); 
					}
				}
			
			}
			
			if ($type_analysis eq "polyploid_polyploid")
			{
				my $color = "gray";
				$pos = $pos / $global_scale;
	            
	            $diagram_img->fgcolor($color);
	            $diagram_img->bgcolor($color);
	            $diagram_img->rectangle( $margin_left + ($pos*$scale) - ($gene_width / 2),
	                                                $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1,
	                                                $margin_left + ($pos*$scale) + ($gene_width / 2),
	                                                $margin_top + $chrom_width + (($section_size + $margin_between_section) * $num_chrom) - 1
	                                                );
	                                                
	                                                
	            
				my $proportion_equal = $proportions_categories{$gene}{"equal"};
				my $proportion_diff = $proportions_categories{$gene}{"difference"};
	
				my $draw = 0;
				if (defined $previous_x_equal)
				{
					$draw = 1;
				}
				
				
				##################################################
				# SNP category : equality between 2 polyploids
				##################################################
				if ($draw)
				{
					$diagram_img->moveTo($previous_x_equal,$previous_y_equal);
				}    
				$previous_x_equal = $margin_left + ($pos*$scale) - 1;    
				$diagram_img->setThickness(2);                          
		        $diagram_img->fgcolor("red");
	            $diagram_img->bgcolor("red");
		        $previous_y_equal = $margin_top + (($section_size + $margin_between_section) * $num_chrom) + 1 - ($proportion_equal * 20) - 7;
		        if ($draw)
				{
		        	$diagram_img->lineTo($previous_x_equal,$previous_y_equal);  
				}
			}

			$nombre_genes++;
        }

		$num_chrom++;
}     

if ($type_analysis eq "polyploid_polyploid")
{
	$diagram_img->moveTo(5,$y_end);
	$diagram_img->setThickness(2); 
	$diagram_img->fgcolor("red");
	$diagram_img->bgcolor("red");
	$diagram_img->lineTo(25,$y_end);  
	$diagram_img->fgcolor("black");
    $diagram_img->moveTo(30,$y_end + 5);
    $diagram_img->fontsize(12);
    $diagram_img->font('Times');
    $diagram_img->string("% SNP where P1 = P2");
}
elsif ($type_analysis eq "polyploid_diploid")
{
	if ($display_cat eq "yes")
	{
		$diagram_img->moveTo(5,$y_end);
		$diagram_img->setThickness(2); 
		$diagram_img->fgcolor("orange");
		$diagram_img->bgcolor("orange");
		$diagram_img->lineTo(25,$y_end);  
		$diagram_img->fgcolor("black");
	    $diagram_img->moveTo(30,$y_end + 5);
	    $diagram_img->fontsize(12);
	    $diagram_img->font('Times');
	    $diagram_img->string("% SNP type 1");
	    
	    $diagram_img->moveTo(5,$y_end + 20);
		$diagram_img->setThickness(2); 
		$diagram_img->fgcolor("purple");
		$diagram_img->bgcolor("purple");
		$diagram_img->lineTo(25,$y_end + 20);  
		$diagram_img->fgcolor("black");
	    $diagram_img->moveTo(30,$y_end + 25);
	    $diagram_img->fontsize(12);
	    $diagram_img->font('Times');
	    $diagram_img->string("% SNP type 2");
	    
	    $diagram_img->moveTo(5,$y_end + 40);
		$diagram_img->setThickness(2); 
		$diagram_img->fgcolor("green");
		$diagram_img->bgcolor("green");
		$diagram_img->lineTo(25,$y_end + 40);  
		$diagram_img->fgcolor("black");
	    $diagram_img->moveTo(30,$y_end + 45);
	    $diagram_img->fontsize(12);
	    $diagram_img->font('Times');
	    $diagram_img->string("% SNP type 3 or 4");
	    
	    $diagram_img->moveTo(5,$y_end + 60);
		$diagram_img->setThickness(2); 
		$diagram_img->fgcolor("red");
		$diagram_img->bgcolor("red");
		$diagram_img->lineTo(25,$y_end + 60);  
		$diagram_img->fgcolor("black");
	    $diagram_img->moveTo(30,$y_end + 65);
	    $diagram_img->fontsize(12);
	    $diagram_img->font('Times');
	    $diagram_img->string("% SNP type 5");
	}
	else
	{
		$diagram_img->moveTo(5,$y_end);
        $diagram_img->setThickness(2);
        $diagram_img->fgcolor("red");
        $diagram_img->bgcolor("red");
        $diagram_img->lineTo(25,$y_end);
        $diagram_img->fgcolor("black");
        $diagram_img->moveTo(30,$y_end + 5);
        $diagram_img->fontsize(12);
        $diagram_img->font('Times');
        $diagram_img->string("% SNP Class 5 per gene (SNP Intra-Diploids = SNP Intra-Polyploid)");
	}

	$diagram_img->moveTo(5,$y_end + 30);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("Estimate of subgenomic contribution to the transcriptome for each gene (%G2)");


	$diagram_img->moveTo(25,$y_end + 45);
	$diagram_img->setThickness(10);
	$diagram_img->fgcolor($red);
	$diagram_img->bgcolor($red);
	$diagram_img->lineTo(30,$y_end + 45);
	$diagram_img->fgcolor("black");
	$diagram_img->moveTo(35,$y_end + 50);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("0-30%");

	$diagram_img->moveTo(95,$y_end + 45);
	$diagram_img->setThickness(10);
	$diagram_img->fgcolor($red_light);
	$diagram_img->bgcolor($red_light);
	$diagram_img->lineTo(100,$y_end + 45);
	$diagram_img->fgcolor("black");
	$diagram_img->moveTo(105,$y_end + 50);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("30-40%");

	$diagram_img->moveTo(165,$y_end + 45);
	$diagram_img->setThickness(10);
	$diagram_img->fgcolor($orange);
	$diagram_img->bgcolor($orange);
	$diagram_img->lineTo(170,$y_end + 45);
	$diagram_img->fgcolor("black");
	$diagram_img->moveTo(175,$y_end + 50);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("40-60%");

	$diagram_img->moveTo(235,$y_end + 45);
	$diagram_img->setThickness(10);
	$diagram_img->fgcolor($orange_light);
	$diagram_img->bgcolor($orange_light);
	$diagram_img->lineTo(240,$y_end + 45);
	$diagram_img->fgcolor("black");
	$diagram_img->moveTo(245,$y_end + 50);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("60-70%");

	$diagram_img->moveTo(305,$y_end + 45);
	$diagram_img->setThickness(10);
	$diagram_img->fgcolor($yellow);
	$diagram_img->bgcolor($yellow);
	$diagram_img->lineTo(310,$y_end + 45);
	$diagram_img->fgcolor("black");
	$diagram_img->moveTo(315,$y_end + 50);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("70-100%");

	$diagram_img->moveTo(25,$y_end + 60);
    $diagram_img->setThickness(10);
    $diagram_img->fgcolor("gray");
	$diagram_img->bgcolor("gray");
	$diagram_img->lineTo(30,$y_end + 60);
	$diagram_img->fgcolor("black");
	$diagram_img->moveTo(35,$y_end + 65);
	$diagram_img->fontsize(12);
	$diagram_img->font('Times');
	$diagram_img->string("No ratio information, no SNP class 5 in this gene");

}                             


open( DIAGRAM_PICT, ">$output_png" );
binmode(DIAGRAM_PICT);
print DIAGRAM_PICT $diagram_img->png;
close DIAGRAM_PICT;
