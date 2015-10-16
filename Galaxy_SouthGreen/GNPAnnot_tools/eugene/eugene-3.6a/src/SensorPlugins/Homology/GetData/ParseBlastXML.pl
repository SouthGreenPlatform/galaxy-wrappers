#!/usr/bin/perl

use XML::Parser;
use Data::Dumper;

my $p1 = new XML::Parser((Style => 'Tree',
		       ErrorContext => 2,
		       Handlers => {Start => \&handle_start,
				    Char  => \&handle_char_data,
				    End   => \&handle_end,
				    Final => \&handle_doc_end}));

my $file = $ARGV[0];

my %hits;      #Hash pour la gestion des hits
my $hit;       #Hash pour le hit courrant
my $current;   #Tag courrant
my $query_id;  #Id de la séquence query
my $query_len; #Taille de la séquence query
my $query_def; #Header de la séquence query

my $nb_seq;    #Nombre de séquence dans la base
my $nb_aa;     #Nombre d'aa dans la base
my $lambda;    #Le lambda (gap)
my $kappa;     #Le kappa  (gap)
my $entropy;   #L'entropie (gap)

#Variables locales pour stocker les infos de la HSP courrante
my $hsp_num;       #Num de la Hsp dans le couple query/subject courrant
my $hsp_bit_score; #Bit-score de la Hsp courrante
my $hsp_score;     #Score de la Hsp courrante
my $hsp_evalue;	   #E-value de la Hsp courrante
my $hsp_query_from;#Début de la hsp sur la query
my $hsp_query_to;  #Fin de la hsp sur la query
my $hsp_hit_from;  #Début de la hsp sur la subject
my $hsp_hit_to;	   #Fin de la hsp sur la subject
my $hsp_identity;  #nb d'aa identiques dans la hsp
my $hsp_positive;  #nb d'aa similaires dans la hsp
my $hsp_gaps;	   #nb de gaps dans la hsp
my $hsp_align_len; #taille de la hsp (avec gap)
my $hsp_density;   #densité (?)
my $hsp_qseq;      #hit sequence
my $hsp_query_frame;

my $hit_def;

my $tree = $p1->parsefile($file);

#####################################################################################
sub handle_start
  {
    my ($expat,$elem,%attrs) = @_;
    $current = $elem;
    if($elem eq "Hit") 
      {
	$hit = {};$hit->{'Hsp'} = {} ;
	$hit_def        = "";
      }
    if($elem eq "Hsp")
      {
	#Initialisation des variables des Hsp 
	#Ces variables permettent de stocker le texte des éléments
	$hsp_num        = "";
	$hsp_bit_score  = "";
	$hsp_score      = "";	
	$hsp_evalue     = "";	
	$hsp_query_from = "";
	$hsp_query_to   = "";
	$hsp_hit_from   = "";
	$hsp_hit_to     = "";	
	$hsp_identity   = "";
	$hsp_positive   = "";
	$hsp_gaps       = "";	
	$hsp_align_len  = "";
	$hsp_density    = "";
	$hsp_qseq       = "";
	$hsp_query_frame= "";
      }
  }
sub handle_char_data
  {
    my ($expat,$text) = @_;
    ## Récupération des info de la query
    if($current eq "BlastOutput_query-ID"  ) {                  $query_id   .= $text;}
    if($current eq "BlastOutput_query-def" ) {                  $query_def  .= $text;}
    if($current eq "BlastOutput_query-len" ) {$text =~ s/\s*//g;$query_len  .= $text;}
    ## Récupération des info du hit
    if($current eq "Hit_id"        ) {$text =~ s/\s*//g;$hit->{ $current } .= $text;}
    if($current eq "Hit_def"       ) {                  $hit->{ $current } .= $text;}
    if($current eq "Hit_num"       ) {$text =~ s/\s*//g;$hit->{ $current } .= $text;}
    if($current eq "Hit_accession" ) {$text =~ s/\s*//g;$hit->{ $current } .= $text;}
    if($current eq "Hit_len"       ) {$text =~ s/\s*//g;$hit->{ $current } .= $text;}
    ## Récupértion des statistiques du run
    if($current eq "Statistics_db-num"  ) {$text =~ s/\s*//g;$nb_seq  .= $text;}
    if($current eq "Statistics_db-len"  ) {$text =~ s/\s*//g;$nb_aa   .= $text;}
    if($current eq "Statistics_kappa"   ) {$text =~ s/\s*//g;$kappa   .= $text;}
    if($current eq "Statistics_lambda"  ) {$text =~ s/\s*//g;$lambda  .= $text;}
    if($current eq "Statistics_entropy" ) {$text =~ s/\s*//g;$entropy .= $text;}
    ## Gestion des Hsp multiples
    if($current eq "Hsp_num")          {$text=~s/\s*//g;$hsp_num .= $text;       }
    if($current eq "Hsp_bit-score")    {$text=~s/\s*//g;$hsp_bit_score .= $text; }
    if($current eq "Hsp_score")        {$text=~s/\s*//g;$hsp_score .= $text;     }
    if($current eq "Hsp_evalue")       {$text=~s/\s*//g;$hsp_evalue .= $text;    }
    if($current eq "Hsp_query-from")   {$text=~s/\s*//g;$hsp_query_from .= $text;}
    if($current eq "Hsp_query-to")     {$text=~s/\s*//g;$hsp_query_to .= $text;  }
    if($current eq "Hsp_hit-from")     {$text=~s/\s*//g;$hsp_hit_from .= $text;  }
    if($current eq "Hsp_hit-to")       {$text=~s/\s*//g;$hsp_hit_to .= $text;    }
    if($current eq "Hsp_identity")     {$text=~s/\s*//g;$hsp_identity .= $text;  }
    if($current eq "Hsp_positive")     {$text=~s/\s*//g;$hsp_positive .= $text;  }
    if($current eq "Hsp_gaps")         {$text=~s/\s*//g;$hsp_gaps .= $text;      }
    if($current eq "Hsp_align-len")    {$text=~s/\s*//g;$hsp_align_len .= $text; }
    if($current eq "Hsp_density")      {$text=~s/\s*//g;$hsp_density .= $text;   }
    if($current eq "Hsp_qseq")         {$text=~s/\s*//g;$hsp_qseq    .= $text;   }
    if($current eq "Hsp_query-frame")  {$text=~s/\s*//g;$hsp_query_frame .=$text;}
  }
sub handle_end
  {
    my ($expat,$elem) = @_;
    if($elem eq "Hit")
      {
	my $ident = $hit->{Hit_num} . "_" . $hit->{Hit_id};
	$hits{ $ident } = $hit;
      }
    if($elem eq "Hsp")
      {
	$hit->{'Hsp'}->{$hsp_num}                     = {};
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_num"}        = $hsp_num;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_bit-score"}  = $hsp_bit_score;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_score"}      = $hsp_score;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_evalue"}     = $hsp_evalue;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_query-from"} = $hsp_query_from;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_query-to"}   = $hsp_query_to;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_hit-from"}   = $hsp_hit_from;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_hit-to"}     = $hsp_hit_to;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_identity"}   = $hsp_identity;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_positive"}   = $hsp_positive ;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_gaps"}       = $hsp_gaps;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_align-len"}  = $hsp_align_len;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_density"}    = $hsp_density;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_qseq"}       = $hsp_qseq;
	$hit->{'Hsp'}->{$hsp_num}->{"Hsp_query-frame"}= $hsp_query_frame;
      }
  }
sub handle_doc_end
  {
    #Pour chaque hit
    foreach my $key (keys(%hits)) 
      {
	#Pour Chaque HSP du couple query/hit
	foreach my $hspkey (keys %{$hits{$key}->{'Hsp'}})
	  {
	    my @tmp;
	    #Création de l'identificateur du sommet source
	    my $Q      = parse_header($query_def); 
	    #Création de l'identificateur du sommet cible
	    my $S      = parse_header($hits{$key}->{"Hit_def"});

	    #Données 
	    my $Ql     = $query_len;
	    my $Sl     = $hits{$key}->{"Hit_len"};

	    my $s      = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_score"};
	    my $bs     = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_bit-score"};

	    my $Qsh    = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_query-from"};
	    my $Qeh    = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_query-to"};
	    my $Ssh    = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_hit-from"};
	    my $Seh    = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_hit-to"};
	    my $qs     = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_qseq"};

	    my $lq     = ($Qeh - $Qsh);       # taille du HSP sur la query
	    my $ls     = ($Seh - $Ssh);       # taille du HSP sur la subject
	    my $max    = 1+($lq>$ls?$lq:$ls); # la plus grande des 2 HSP
	    my $ev     = sprintf("%.0e",$hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_evalue"});
	    if ($ev==0) { $ev="0.0"};
	    my $qf     = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_query-frame"};
	    $qf = ( ($qf>0) ? "+"."$qf" : $qf );

	    my $tta    = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_identity"};
	    my $ttb    = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_positive"};

	    my $id     = sprintf("%.2f",($tta / $max * 100));
	    my $sm     = sprintf("%.2f",($ttb / $max * 100));
	    my $alen   = $hits{$key}->{'Hsp'}->{$hspkey}->{"Hsp_align-len"};
	    #Affiche les informations du HSP sur une seule ligne
#	    print $Q," ",$S," ",$Ql," ",$Sl," ",$s," ",$bs," ",$Qsh," ",$Qeh,
#	          " ",$Ssh," ",$Seh," ",$id," ",$sm," ",$alen,"\n";
	    if ( ("$Q" ne "$S") ||
		 ($Qsh != $Ssh) ||
		 ($Qeh != $Seh) )
	      {
		print $Qsh," ",$Qeh," ",$s," ",$ev," ",$qf," ",$S,"_",$Ssh,"_",$Seh," ",$Ssh," ",$Seh," ",$qs,"\n";
	      }
	  }
      }
  }

sub parse_header
  {
    my $head = shift;
#    $head ~= s/,//g;
    my @tok = split /\s+/, $head;
    return $tok[0];
  }
