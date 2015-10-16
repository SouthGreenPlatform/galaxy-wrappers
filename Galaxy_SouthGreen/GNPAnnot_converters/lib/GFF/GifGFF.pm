# Perl Module for GifGFF
#
# Cared for by Tim Hubbard <th@sanger.ac.uk> and 
# also watered by Richard Bruskiewich <rbsk@sanger.ac.uk>
#
#Copyright Genome Research Limited (1997). Please see information on licensing in LICENSE

package GFF::GifGFF;

use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use Exporter;
use Carp;
use GFF;
use GD;
use strict;

my $VERSION = '1.14' ; 

@EXPORT_OK = qw();

# @ISA has our inheritance.
@ISA = ( 'Exporter' );

my %colours ; # Global Declaration of the Colour Palette (see sub _colour_init() and _colour())

sub _colour_init ;
sub _colour ;

# Default font
my $FONT = gdSmallFont ; # 6x12

1;  # says use was ok

# creates a transformation matrix function by closure
# Vertical plots are reflected on a diagonal axis, running through (0,0)
# Either horizontal and vertical plots are then translated to 
# the specified GIF graph (h==x_origin,k==y_origin) point of origin
sub _transformation_matrix {
    my ($v)=shift ;
    my $h = shift ; $h=0 if !defined($h);
    my $k = shift ; $k=0 if !defined($k);
    if(defined($v) && $v) {
	return sub {
	    my($x,$y)=@_;
	    if(!defined($x)){
		1;# signal vertical mode
	    }else{
		($y+$h,$x+$k);
	    }
	} ;
    } else {
	return sub {
	    my($x,$y)=@_;
	    if(!defined($x)){
		0;# signal horizontal mode
	    }else{
		($x+$h,$y+$k); # identity matrix
	    }
	} ;
    }
}

sub gff2gif{
    my %args=(
	      x_origin=>0,
	      y_origin=>0,
	      width=>400,
	      margin=>20,
	      scale=>500,
	      yseparation=>20,
	      expandmap=>2,
	      mapA_offset => 0,
	      mapB_offset => 0,
	      compress=>1,
	      bump=>1,
	      bumpspace=>2,
	      intronstyle=>'straight',
	      barwidth=>3,
	      barcolour=>'black',
	      gridcolour=>'gray',
	      bumplabel=>0,
	      bumplabelseparation=>1,
	      resolution=>100,
	      font=>'small',
	      @_,
	      );
    # required
    my $gff=$args{'gff'};
    my $rhgifgffdata=$args{'layout'};
    my $group_features=$args{'filter'};

    # One of:

    my $OUT=$args{'io'};

    # *OR*

    my $im       = $args{'image'} ;   # existing user-provided GD::Image?
    my $palette  = $args{'palette'} ; # GD::Image $im->colorAllocate() values indexed with string labels

    # (h,k) plot origin (for axis of reflection during 'vertical' mode)
    my $x_origin = $args{'x_origin'} ;
    $x_origin = 0 if !defined($x_origin) ;
    my $y_origin = $args{'y_origin'} ;
    $y_origin = 0 if !defined($y_origin) ;

    # fonts
    if($args{'font'} eq 'tiny'){
	$FONT=gdTinyFont;
    }elsif($args{'font'} eq 'small'){
	$FONT=gdSmallFont;
    }elsif($args{'font'} eq 'medium'){
	$FONT=gdMediumBoldFont;
    }elsif($args{'font'} eq 'giant'){
	$FONT=gdGiantFont;
    }

    # from defaults
    my $maxx=$args{'width'};

    my $margin=$args{'margin'};
    my ($leftmargin,$rightmargin,$topmargin,$bottommargin) ;
    if(defined($args{'leftmargin'})) {
	$leftmargin = $args{'leftmargin'};
    } else {
	$leftmargin = $margin ;
    }

    # number of $FONT characters the left margin can hold
    my $lmchars = int($leftmargin/$FONT->width) ; 

    if(defined($args{'rightmargin'})) {
	$rightmargin = $args{'rightmargin'};
    } else {
	$rightmargin = $margin ;
    }
    if(defined($args{'topmargin'})) {
	$topmargin = $args{'topmargin'};
    } else {
	$topmargin = $margin ;
    }
    if(defined($args{'bottommargin'})) {
	$bottommargin = $args{'bottommargin'};
    } else {
	$bottommargin = $margin ;
    }

    my $yseparation=$args{'yseparation'};
    my $intronstyle=$args{'intronstyle'};

    # expand the map of exon features
    my $expandmap=$args{'expandmap'};

    # bulk offsets for coordinats of feature rows on
    # a pair of clone maps specified 
    # by 'A' or 'B' suffixes on feature labels
    my $mapA_offset=$args{'mapA_offset'};
    my $mapB_offset=$args{'mapB_offset'};

    # optional: filehandle to write out map information
    my $MOUT=$args{'mio'};

    my $compress=$args{'compress'};
    my $bump=$args{'bump'};
    my $bump_space=$args{'bumpspace'};

    # horizontal versus vertical plot?
    # select transformation matrix accordingly
    my $vertical=$args{'vertical'};
    $vertical=0 if !defined($vertical);
    my $T=_transformation_matrix($vertical,$x_origin,$y_origin);

    # optional: left margin row labels
    my $rowlabels=$args{'rowlabels'};
    
    # optional: gridlines
    my $gridlines=$args{'gridlines'};

    # get range in gff
    my($end,$start);
    if($args{'start'} && $args{'end'}){
	$start=$args{'start'};
	$end=$args{'end'};
	# filter gff
	my $gffold=$gff;
	$gff=new GFF::GeneFeatureSet;
	foreach my $gf ($gffold->eachGeneFeature()){
	    my $gfs=$gf->start;
	    my $gfe=$gf->end;
	    # remove completely out of range features
	    next if $gfs>$end;
	    next if $gfe<$start; 
	    # truncate those partly out of bounds
	    if($gfs<$start){$gf->start($start);}
	    if($gfe>$end){$gf->start($end);}
	    $gff->addGeneFeature($gf);
	}
    }else{
	($start,$end)=$gff->min_max_range;
    }
    my $len=($end-$start+1);

    # apply filter
    my %group=$gff->group($group_features);

    # extend groups for direction
    # label x becomes xR or xF if this type of labelling exists
    my %labelorder;
    # if bumping subsets, can't bump globalset
    my %non_bump;
    for(my $i=0;$i<scalar(@$rhgifgffdata);$i++){
	my($order,$height,$label,$colour1,$colour2,$fbox,$layout_filter)=@{$rhgifgffdata->[$i]};
	next unless $label;
	$labelorder{$label}=[$order,$layout_filter,$fbox];
	if($colour2=~/^(VB?)?([FR])/){
	    my $d=$2;
	    my $newlabel=$label.$d;
	    $non_bump{$label}=1;
	    if($group{$label} && !$group{$newlabel}){
		my $gff=$group{$label};
		my $gffnew=new GFF::GeneFeatureSet;
		foreach my $gf ($gff->eachGeneFeature()){
		    if($gf->strand eq '+'){
			next if($d eq 'R');
		    }else{
			next if($d eq 'F');
		    }
		    $gffnew->addGeneFeature($gf);
		}
		$group{$newlabel}=$gffnew;
		$labelorder{$newlabel}=[$order,$layout_filter,$fbox];
#		print "used new label $newlabel\n";
	    }
	    $rhgifgffdata->[$i]->[2]=$newlabel;
	}
    }

    # apply bumping within each group
    my %bump_level;
    # loop over groups
    # if groups are to be displayed split (F/R) then need to split at this level
    foreach my $group (keys %group){
	next if $non_bump{$group};
	my @mask;
	my $gff=$group{$group};
	my %group2;
#	print "$group\n";
	my $bump_filter;
	if($bump_filter=$labelorder{$group}->[1]){
	    next if $bump_filter eq 'none';
	    %group2=$gff->group($bump_filter);
	    foreach my $group2 (reverse sort keys %group2){
		my $gff2=$group2{$group2};
		my($start2,$end2)=$gff2->min_max_range;
		my $i=&_mask(\@mask,$start2,$end2,$end,$args{'resolution'});
		if(!defined($bump_level{$group}) or $bump_level{$group}<$i){$bump_level{$group}=$i;}
		# save bump level for this groups
		foreach my $gf ($gff2->eachGeneFeature()){
		    $gf->addMember($i,'bump_level');
		}
	    }
	}elsif($bump || ($group=~/^label/ && $args{'bumplabel'})){
	    my $fboxoff=$labelorder{$group}->[2];
	    # either global bumping or label bumping
	    my $i=1;
	    foreach my $gf ($gff->eachGeneFeature()){
		# fix for labels
		my($start2,$end2)=($gf->start,$gf->end);
		if($args{'bumplabel'} && $group=~/^label/){
		    my $text;
		    if($fboxoff){
			$text=&$fboxoff($gf);
		    }else{
			$text=$gf->seqname;
		    }
		    $end2=$start2+length($text)*$FONT->width*$args{'scale'}*$args{'bumplabelseparation'};
#		    print STDERR "got here: $start2, $end2, ".length($text)." $args{'scale'}\n";
		}
		my $i=&_mask(\@mask,$start2,$end2,$end,$args{'resolution'});
		# increment the max bump_level for this group
		if(!$bump_level{$group} or $bump_level{$group}<$i){$bump_level{$group}=$i;}
		$gf->addMember($i,'bump_level');
	    }
	}
	# now fix direction
	# if order is -ve, should bump in -ve direction
	# if order is +ve, should bump in +ve direction
	if($bump || $bump_filter || $args{'bumplabel'}){
	    # print "<P>### Bumping1! \$group: $group, \$labelorder\{$group\}: $labelorder{$group}, ",
	    #        "\$labelorder\{$group\}\-\>\[0\]: $labelorder{$group}->[0]</P>\n" ;
            my $labelorder = $labelorder{$group}->[0] ;
	    next if defined($labelorder) and $labelorder>=0;
	    my $maxbl=$bump_level{$group};
	    foreach my $gf ($gff->eachGeneFeature()){
		if(defined($gf->getMember('bump_level'))){
		    my $j=$gf->getMember('bump_level');
		    $j=$maxbl-$j;
		    $gf->addMember($j,'bump_level');
		}
	    }
	}
    }

    my $off=$start;
    my $xmin=$leftmargin;
    my $xmax=$maxx-$rightmargin;
    my $xblock=$xmax-$xmin;

    # work out y coordinates
    # 1. get biggest for this position
    my %order;
    foreach my $type (@$rhgifgffdata){
	my($order,$height,$label,$colour1,$colour2)=@$type;

	# if no label, this is a spacer, so include
	my $bump_level = 0 ;
	if(defined($label) && $label && $compress){
	    next unless defined($group{$label}) ; # excludes 'bars'
	    $bump_level = $bump_level{$label} if defined($bump_level{$label}) ;
	}
	# allocate space, including bumping
	#print "<P>### Bumping2! label: $label, height: $height,bump_space: ",
	#      "$bump_space, bump_level\{$label\}: $bump_level{$label},  ",
	#      "\$order: $order, \$order\{$order\}->[0]: $order{$order}->[0]</P>\n" ;
	my $bump_height=$height+($height+$bump_space)*$bump_level;
	
	if(!defined($order{$order}) or $bump_height>$order{$order}->[0]){
	    $order{$order}=[$bump_height,$label];
	    push(@{$order{$order}},1) if(defined($colour2) and $colour2=~/G$/) ;
	}
    }
    # 2. step through positions
    my %y;
    my $theight = 0;
    my $gridheight = 0 ;
    my $grid=1 ;
    my $llabel;
    foreach my $order (sort {$a<=>$b} keys %order){
	my($height,$label,$gridtoggle)=@{$order{$order}};
	# skip spacer if previous was spacer
	next unless((defined($label) and $label) or 
		    (defined($llabel) and $llabel));
	$y{$order}=$theight+$topmargin;
	$theight+=$height;
	$grid=!$grid if(defined($gridtoggle)) ;
	$gridheight+=$height if($grid);
	$llabel=$label;
#	print "$height,$theight,$label\n";
    }

    # decide scale/maxy
    my($maxy,$scale,$yblock,$nblock,$maxheight);
    $yblock=$yseparation+$theight;
    if($args{'height'}){
	# 1. height parameter was supplied - calculate scale to fit this
	$maxheight=$args{'height'};
	$nblock=int(($maxheight-($topmargin+$bottommargin))/$yblock);
	$maxy=$nblock*$yblock+($topmargin+$bottommargin);
	$scale=$len/(($maxx-($leftmargin+$rightmargin))*$nblock);

    } else {
	# 2. base on scale parameter
	$scale=$args{'scale'};
	$nblock=int($len/(($maxx-($leftmargin+$rightmargin))*$scale))+1;
	$maxy=$nblock*$yblock+($topmargin+$bottommargin);
    }
#    print "$maxheight,$maxy,$nblock,$yblock,$scale\n";
    # setup image, if none already provided
    if(!defined($im)) {
	$im = GD::Image->new(&$T($maxx,$maxy))
	            || die "cannot create GD object";
	# initialize palette
	&_colour_init($im,$palette) ;

	# background + interlacing
	$im->transparent(&_colour('white')) if defined($args{'transparent'});
	$im->interlaced('true');
	$im->filledRectangle(0,0,&$T($maxx,$maxy), &_colour('white'));
    } else {
	# initialize palette
	&_colour_init($im,$palette) ;
    }

# draw gridlines on canvas

    if(defined($gridlines)) {
	my $spacing= ($gridlines*1000)/$scale ; # in pixels
	for(my $pos=$leftmargin;$pos<$xmax;$pos+=$spacing) {    
	    $im->line(&$T($pos,$topmargin),&$T($pos,$topmargin+$gridheight),&_colour($args{'gridcolour'}));
	}
    }


# draw bars on canvas
    my %dummy;
    if($group{'bar'}){
	my $col = $args{'barcolour'} ;
	&_draw_set($im,\$group{'bar'},$topmargin,$theight,0,
		   $col,$col,
		   $xblock,$yblock,$off,$scale,$leftmargin,$xmin,$xmax,
		   $MOUT,0,$expandmap,
		   $bump_space,
		   $intronstyle,'bar',\%args,
		   $T,$maxy-($topmargin+$bottommargin),
		   $mapA_offset,$mapB_offset);
    }

# will compute map positions on canvas
    if($group{'map_position'}){
	# I exploit drawset for convenience sake but I'm not drawing anything ;-)
	&_draw_set($im,\$group{'map_position'},$topmargin,$theight,0,
		   'black','black',
		   $xblock,$yblock,$off,$scale,$leftmargin,$xmin,$xmax,
		   $MOUT,0,$expandmap,
		   $bump_space,
		   $intronstyle,'map_position',\%args,
		   $T,$maxy-($topmargin+$bottommargin),
		   $mapA_offset,$mapB_offset);
    }

    # draw features on canvas
    foreach my $type (@$rhgifgffdata){
	my($order,$height,$label,$colour1,$colour2,$fboxoff)=@$type;
	# skip empty sets
	next unless $label;
	next unless $group{$label};
	# print "<P>$label,$height,$colour1,$colour2\n</P>";
	if(!(defined(&_colour($colour1)) or $colour2 eq 'P')){
	    print "<P>colour 1: $colour1 not defined - using black</P>\n";
	    $colour1='black';
	}
	if(!($colour2 =~ /^((VB?)?([FR])?G?|P)$/ or defined(&_colour($colour2)))){
	    print "<P>colour 2: $colour2 not defined - using black</P>\n";
	    $colour2='black';
	}
# FIXME
# fs switch (for alternating contigs) not connected to API
#	print "$label,$order,$height\n";
	&_draw_set($im,\$group{$label},$y{$order},$height,0,
		   $colour1,$colour2,
		   $xblock,$yblock,$off,$scale,$leftmargin,$xmin,$xmax,
		   $MOUT,$fboxoff,$expandmap,
		   $bump_space,
		   $intronstyle,$label,\%args,
		   $T,$maxy-($topmargin+$bottommargin)
		   ,$mapA_offset,$mapB_offset);
#
# Draw optional $label in left margin
#
        if( defined($rowlabels) and $rowlabels and 
	    (ref($rowlabels) =~ /HASH/) and exists($rowlabels->{$label}) ) {

	       my $rspec = $rowlabels->{$label} ;

	       # $rspec is of form label<url>, where <url> is optional
	       my ($rlabel, $colour, $rurl) = split ',', $rspec ;
	       #print "<P>rowlabel: label($label:'$rlabel'[',length($rlabel),']), colour($colour?$colour1), URL($rurl)<P>\n";
	       $colour = $colour1 if !(defined($colour) && $colour);

               # print the label
	       $rlabel = substr($rlabel,0,$lmchars) if(length($rlabel)>$lmchars) ;
               my $tpos = $leftmargin ; 
	       my $mid=($height-$FONT->height)/2;
	       my $y1 = $y{$order}+($mid>0?$mid:0) ;
	       my $y2 = $y1+$FONT->height ;
	       my ($acx1,$acy1,$acx2,$acy2);
	       if($vertical) {
		   $tpos -= $FONT->width ; 
		   $tpos = ($tpos >= 0 ? $tpos : 0) ;
		   ($acx1,$acy1)=&$T($tpos,$y1);
		   $im->stringUp($FONT,$acx1,$acy1,$rlabel,&_colour($colour));
	       } else {
		   $tpos -= (length($rlabel)+1)*$FONT->width ; 
		   $tpos = ($tpos >= 0 ? $tpos : 0) ;
		   ($acx1,$acy1)=&$T($tpos,$y1);
		   $im->string($FONT,$acx1,$acy1,$rlabel,&_colour($colour));
	       }
	       ($acx2,$acy2)=&$T($tpos+length($rlabel)*$FONT->width,$y2+1);
	       unless($vertical) {
		   if($MOUT) {
                       # optional url
		       my $out="<area coords=\"$acx1,$acy1,$acx2,$acy2\"\n" ;
		       # link
		       if(defined($rurl)){
			   $out.="\t href=\"$rurl\"\n";
			   $im->line($acx1,$acy2,$acx2,$acy2, &_colour($colour)) ; # simulate URL link underline
		       } else {
			   $out.="\t href=\"#\"\n";
		       }
		       # and mouseover label
		       if(defined($rlabel)) {
			   $out.="\t onMouseOver=\"window.status=\'$rlabel\'; return true\"\n";
			   $out.="\t onMouseOut=\"window.status=\'\'; return true\"\n";
			   $out.="\t name=\"$rlabel\"\n";
			   $out.="\t alt=\"$rlabel\"\n";
		       }
		       $out.=">\n";
		       print $MOUT $out;
		   }
	       }
	   }
    }

# output gif
    print $OUT $im->gif if $OUT ;

# return maximum extent of map (in pixels), 
# assumed to be anchored on (x_origin,y_origin)
    return &$T($maxx,$maxy);
}

# returns bump level + max bump level
sub _mask{
    my($ramask,$start2,$end2,$end,$resolution)=@_;
    if($resolution){
	$start2=int($start2/$resolution)+1;
	$end2=int($end2/$resolution)-1;
	if($start2>$end2){$end2=$start2;}
	$end=int($end/$resolution)+1;
    }
    my $i=0;
    my $len=$end2-$start2+1;
    {
	if($$ramask[$i]){
	    # increment index of mask array, until find some space
	    if(substr($$ramask[$i],$start2,$len)=~/x+/){
		$i++;
		redo;
	    }
	}else{
	    # if new $i, setup blank string
	    $$ramask[$i]='-' x $end;
	}
    }
    # at this gene to mask
    substr($$ramask[$i],$start2,$len)='x' x $len;
    return $i;
}

# draw clones as tiled
sub _draw_set{
    my($im,$rgff,$ys,$height,$fs,
       $col1,$col2,
       $xblock,$yblock,$off,$scale,$leftmargin,$xmin,$xmax,
       $MOUT,$fboxoff,$expandmap,
       $bump_space,$intronstyle,$label2,
       $rhargs,
       $T,$plotheight,
       $mapA_offset,$mapB_offset)=@_;

    my $shift=0;
    foreach my $gf ($$rgff->eachGeneFeature()){
	my $bump_level;
	my $y=$ys;
	if($bump_level=$gf->getMember('bump_level')){
	    $y+=($height+$bump_space)*$bump_level;
	}
#	print "  $ys,$y,".$$rhbump_level{($gf->seqname)}.",".$gf->seqname."\n";

	my ($offset,$offset2) ;
	my $grp = $gf->group_value('MetaGroup') ; # may be undefined
	if(defined($grp)) {
	    # $map?_offset opposes effect of $off, 
	    # to shift map rightwards?
	    if($grp eq 'A') {
		$offset  = $off-$mapA_offset ;
		$offset2 = $off-$mapB_offset ; # complementary for Homol's
	    } else {
		$offset  = $off-$mapB_offset ;
		$offset2 = $off-$mapA_offset ; # complementary for Homol's
	    }
	} else {
	    $offset = $offset2 = $off ;
	}

	#print "<p>off($off), offset($offset), scale($scale), leftmargin($leftmargin)</p>\n";

# convert coords to pixels
	my $start = &_scale($gf->start,$offset,$scale,$leftmargin);
	my $end   = &_scale($gf->end,$offset,$scale,$leftmargin);
# constant thickness for bar
	if($label2 eq 'bar'){
	    my $barwidth=($$rhargs{'barwidth'});
	    if($$rhargs{'barwidth'}){
		if($end>=$xmax){
		    $start=$end-$barwidth;
		}else{
		    $end=$start+$barwidth;
		}
	    }
	}
	my($fhom,$start2,$end2);
	if( ref($gf) eq 'GFF::HomolGeneFeature' or

            # I use the 'Homology' tag rather than 'Target' tag here
            # to avoid problems with features like 'RepeatMasker' repeats,
            # which are also dumped with 'Target' tags from ACEDB
            # In this way, callers of GifGFF must explicitly mark their GFF to be
            # displayed as an homol feature
	    defined($gf->group_value_list('Homology')) ){ 
	    $fhom   = 1;
	    $start2 = &_scale($gf->start2,$offset2,$scale,$leftmargin);
	    $end2   = &_scale($gf->end2,$offset2,$scale,$leftmargin);
	}
	my $group  = $gf->group;
#	my $group  = $gf->dump_group;
	my $label  = $gf->getMember('label');
	my $href   = $gf->getMember('href');
	my $intron = $gf->getMember('intron');
#	print $gf->start.' '.$gf->end." $group $start,$end,$xmax\n";
	my $ys = 0 ;
	my $col ;
	my $vtext ; # undefined unless 'V'

	# drawing rules
	if($col2=~/^V/){ # directive for vertically oriented text
	    if($col2=~/B/){$vtext=-1;}else{$vtext=1;}
	}

	if($col2 eq 'P'){
	    # colour by percentage
	    unless($fhom){
		die "GF ".$gf->seqname." is not an homology GeneFeature - can't use 'P' shading\n";
	    }
	    my $percent = $gf->percentid;
	    my $gb = int((100-$percent)/100*200);
	    if($gb > 200) {$gb = 200;}
	    if(&_colour($col1,$gf) eq 'red'){
		$col = $im->colorAllocate(255,$gb,$gb);
	    }elsif(&_colour($col1,$gf) eq 'green'){
		$col = $im->colorAllocate($gb,255,$gb);
	    }elsif(&_colour($col1,$gf) eq 'blue'){
		$col = $im->colorAllocate($gb,$gb,255);
	    }else{
		print "Colour $col unknown for percentage shading\n";
	    }
	} elsif ($fhom) {
	    # hom's without 'P' drawing - always colour1
	    $col = &_colour($col1,$gf) ;
	} elsif ($gf->strand eq '+'){
	    next if $col2 =~ /^(VB?)?RG?$/ ;
	    $col = &_colour($col1,$gf) ; # works with 'V' alone too
	} elsif ($gf->strand eq '-'){
	    next if $col2  =~ /^(VB?)?FG?$/ ;
	    if($col2 =~ /^((VB?)?R|(VB?))?G?$/ ){ # R possibly with V(B) qualifier, or 'V(B)' alone
		$col = &_colour($col1,$gf) ;
	    } else {
		$col = &_colour($col2,$gf) ; # no specific directive... just a 'reverse' colour
	    }
	} else { # non-specific strand - just use colour #1
	    $col = &_colour($col1,$gf) ;
	}

	# feature range is $start-$end
	# will be drawn in blocks
	{
	    # only call if in range
	    if($start <= $xmax || ($fhom && $start2 <= $xmax) ) {
		&_draw_partial_feature($gf,$im,$fhom,
				       $y,$shift,$ys,$height,$xmax,$xmin,
				       $start,$end,$start2,$end2,
				       $col, &_colour('black'), &_colour('white'), $fboxoff,
				       $intron, $intronstyle, $MOUT, $expandmap,
				       $group, $label, $href, $label2, 
				       $vtext,$T,$plotheight);
	    }
	    # get ready for next block for this gf
	    $start -= $xblock;
	    $end   -= $xblock;
	    $ys    += $yblock;
	    if($fhom) {
		$start2 -= $xblock;
		$end2   -= $xblock;
	    }
	    redo if($end >= $xmin || ($fhom && $end2 >= $xmin));
	}
# special feature to alternate position according to order of gf objects
# useful for contig maps
	if($shift == 0 && $fs) {
	    $shift = $height;
	} else {
	    $shift = 0;
	}
    }
}

sub _draw_partial_feature{
    my($gf,$im,$fhom,
       $y,$shift,$ys,$height,$xmax,$xmin,
       $start,$end,$start2,$end2,
       $col, $black, $white, $fboxoff,
       $intron, $intronstyle, $MOUT, $expandmap,
       $group,$label,$href,$label2,$vtext,$T,$plotheight)=@_;

    my $vertical=&$T(); # without args, T returns transformation mode

# y1<$y2, but y1 is higher than y2 (canvas coordinates from upper left)

#    print "<P>_draw_partial_feature(feature: ",$gf->feature(),", T(\$y): ",(join ',',&$T(0,$y))," \$shift: $shift, \$ys: $ys)\n</P>\n" ;

    my $y1=$y+$shift+$ys;
    my $y2=$y1+$height;
#    print "Draw: $y1,$y2,$label2,$y,$ys,$shift,$height\n";
# special mode for drawing Homology blocks
# note: must draw boundaries correctly
# y2 is start/end; y1 is start2/end2
    if($fhom){
	my $poly = new GD::Polygon;
	$poly->addPt(&$T($start,$y1));
	$poly->addPt(&$T($start2,$y2));
	$poly->addPt(&$T($end2,$y2));
	$poly->addPt(&$T($end,$y1));
	$im->filledPolygon($poly,$col);
#	print "$start,$start2,$end,$end2\n";
# cover up writing in margins
	$im->filledRectangle(&$T(0,$y1),&$T($xmin,$y2),
			     $white);
	$im->filledRectangle(&$T($xmax,$y1),&$T($xmax+100,$y2),
			     $white);
# prepare for drawing boxes
	if($start2<$start){$start=$start2;}
	if($end2>$end){$end=$end2;}
	if($start<$xmin){$start=$xmin;}
	if($end>$xmax){$end=$xmax;}
    } else {
	return if($start>$xmax || $end<$xmin);
	if($start<$xmin){$start=$xmin;}
	if($end>$xmax){$end=$xmax;}

	# special mode for drawing labels
	if($label2=~/^label/){
	    my $text=&$fboxoff($gf);
	    my $xsu=$start;
	    my $ysu=$y1;
	    if($vertical^defined($vtext)) {
		if($vertical) {
		    $xsu+=length($text)*$FONT->width; 
		} else {
		    if($vtext<0) {
			$ysu+=$height ; 
		    } else {
			$ysu+=length($text)*$FONT->width ; 
		    }
		}
		my ($x,$y)=&$T($xsu,$ysu) ;
		$im->stringUp($FONT,$x,$y,$text,$col);
	    } else {
		if($vertical) {
		    # also vtext
		    if($vtext<0) {
			$ysu+=$height-length($text)*$FONT->width ; # bottom aligned 
		    }
		} # else, not vertical and not vtext
		my ($x,$y)=&$T($xsu,$ysu) ;
		$im->string($FONT,$x,$y,$text,$col);
	    }

	} elsif($label2=~/^map_position$/){ 
            # annotate gf with computed pixel position, without drawing
	    my $pixel_pos = [&$T($start,$y1),&$T($end,$y2)] ;
	    $gf->group_value_list('Pixel_Position',$pixel_pos) ;

	} elsif($intron) {
	    # special modes for drawing introns
	    # if style is not recognised, nothing is drawn
	    if($intronstyle eq 'straight'){
		my $ym  = int(($y1+$y2)/2);
		my $y11 = $ym-1;
		my $y21 = $ym+1;
		if($y11 < $y1+1) {$y11 = $y1+1;}
		if($y21 > $y2-1){$y21 = $y2-1;}
		if($intron == 1){
		    $im->filledRectangle(&$T($start,$y11),&$T($end,$y21),
					 $col);
		} elsif ($intron == 2){
		    $im->dashedLine(&$T($start,$y11),&$T($end,$y11),
				    $col);
		    $im->dashedLine(&$T($start,$ym),&$T($end,$ym),
				    $col);
		    $im->dashedLine(&$T($start,$y21),&$T($end,$y21),
				    $col);
		}
	    } elsif ($intronstyle eq 'angle'){
		my $y11 = $y1+1;
		my $y21 = $y2-1;
		if($intron == 1){
		    my $mid = int(($start+$end)/2);
		    $im->line(&$T($start,$y21),&$T($mid,$y11),
			      $col);
		    $im->line(&$T($mid,$y11),&$T($end,$y21),
			      $col);
		} elsif ($intron == 2){
		    my $mid1 = int(($start*2+$end)/3);
		    my $mid2 = int(($start+$end*2)/3);
		    $im->line(&$T($start,$y21),&$T($mid1,$y11),
			      $col);
		    $im->dashedLine(&$T($mid1,$y11),&$T($mid2,$y11),
				    $col);
		    $im->line(&$T($mid2,$y11),&$T($end,$y21),
			      $col);
		}
	    }
	} else {
	    # minimum pixel size for feature
	    my $fend ;
	    if($fboxoff and $fboxoff>1 and
	       $end-$start+1 < $fboxoff ) {
		$fend = $start+$fboxoff-1 ;
	    } else { # no change
		$fend = $end ;
	    }
	    $im->filledRectangle(&$T($start,$y1),&$T($fend,$y2),
				 $col);
	    if(!$fboxoff && $label2 ne 'bar' ){
		$im->rectangle(&$T($start,$y1),&$T($end,$y2),
			       $black);
	    }
	}
    }
# write html map
# since introns tend to be large and exons small then modify
# the map slightly so as to make it easier to hold the mouse over the exon
# FIXME
# this will make a mess if introns are very small - should be checked for
    if($MOUT){
	my($mstart,$mend);
# introns's contract
	if($intron){
	    $mstart = $start + $expandmap;
	    if($mstart > $xmax){$mstart = $xmax;}
	    $mend = $end - $expandmap;
	    if($mend < 1){$mend = 1;}
	    if($mend < $mstart){$mend = $mstart;}
# exon's expand
	}else{
	    $mstart = $start - $expandmap;
	    if($mstart < 1){$mstart = 1;}
	    $mend = $end + $expandmap;
	    if($mend > $xmax){$mend = $xmax;}
	}
	my ($acx1,$acy1)=&$T($mstart,$y1); 
	my ($acx2,$acy2)=&$T($mend,$y2) ;
	my $out="<area coords=\"$acx1,$acy1,$acx2,$acy2\"\n";
# link
	if($href){
	    $out.="\t href=\"$href\"\n";
	}else{
	    $out.="\t href=\"#\"\n";
	}
# label
	if($label){
	    $out.="\t onMouseOver=\"window.status=\'$label\'; return true\"\n";
	    $out.="\t onMouseOut=\"window.status=\'\'; return true\"\n";
	    $out.="\t name=\"$label\"\n";
	}
	$out.=">\n";
	print $MOUT $out;
    }
}

sub _draw_gridline {
    my ($im,$pos,$top,$bottom,$T)=@_ ;
    $im->line(&$T($pos,$top),&$T($pos,$bottom),gdStyled);
}

my %gifgffcolours = (
    'white'       => [255,255,255],
    'gray'        => [192,192,192],
    'black'       => [0,0,0],
    'red'         => [255,0,0],
    'lightred'    => [255,128,128],
    'darkred'     => [128,0,0],
    'green'       => [0,255,0],
    'lightgreen'  => [64,192,64],
    'darkgreen'   => [0,128,0],
    'blue'        => [0,0,255],
    'lightblue'   => [170,216,255],
    'darkblue'    => [0,0,128],
    'yellow'      => [255,215,0],
    'brown'       => [218,165,32],
    'magenta'     => [255,0,255],
    'darkmagenta' => [128,0,128],
    'cyan'        => [60,220,200],
    'darkcyan'    => [0,128,128],
) ;

###################################
# Initialize colour palette ##############
###################################
sub _colour_init {
    my $im           = shift ;
    my $user_palette = shift ;

    my $colour ;
    # record pre-allocated user colours first, if specified
    if(defined($user_palette) and ref($user_palette) =~ /HASH/) {
	foreach $colour (keys %{$user_palette}) { 
	    my $ucol = $user_palette->{$colour} ;
	    if(ref($ucol) =~ /ARRAY/) { 
                # An RGB ARRAY description? Allocate!
		$colours{$colour} = $im->colorAllocate(@{$ucol})  ;
	    } else {
		# I assume that the user has created pre-allocated colours for
		# the GD::Image they precreated and handed to me
		# so $ucol is a GD::Image colorAllocate index
		$colours{$colour} = $ucol  ;
	    }
        }
    } 

    # then allocate missing gifgff colours
    foreach $colour (keys %gifgffcolours) {
	if(!exists($colours{$colour})) {
	    $colours{$colour} = $im->colorAllocate(@{$gifgffcolours{$colour}}) ;
	}
    }
}

sub _colour {
    my $colour = shift ;
    if(ref($colour) =~ /CODE/) {
        my $gf = shift ;  # current $gf being drawn
        return &{$colour}(\%colours,$gf) ;
    } elsif(exists($colours{$colour})) {
        return $colours{"$colour"} ;
    } else {
	return undef; 
    }
}

# draw real genes
sub _scale{
    my($in,$off,$scale,$leftmargin)=@_;
    return int((($in-$off)/$scale)+$leftmargin);
}

__DATA__

=head1 NAME

GFF::GifGFF.pm - Perl extension for graphing GFF Gene Features

=head1 SYNOPSIS

use GFF::GifGFF ;

=head1 DESCRIPTION

GifGFF.pm contains a single non OO routine to create a gif image of
data in GFF objects based on a few global parameters and a drawing
specification array which selects colours etc.

Note: C<beta code> this module is evolving rapidly so I don't
guarantee that future versions will be backwards compatible.

=head1 AUTHOR

Created by B<Tim Hubbard> email th@sanger.ac.uk

Minor revisions by B<Richard Bruskiewich> email rbsk@sanger.ac.uk

=head1 SUBROUTINES

=head2 B<gff2gif>

This call write gif data to a filehandle.  Parameters and data are
passed as hash elements.  The default requirements are (1) a gff file,
(2) a gff filter subroutine that can be applied to group the features
that should be drawn under different labels, (3) a layout array that
specifies how to draw/colour different labels, (4) a reference to an
open filehandle for gif file.  Additional parameters control global
features such as scaling, canvas size, display range (to display only
a region of the gff) etc.  Information embedded within each
GeneFeature object specifies optional map html for mouseover labels
and href links.  Special functionality to link exons, draw bars to
separate contigs, bump overlapping features etc. is embedded.

=head2 Example

Here is an example of the use of GifGFF.pm This example constructs a
gff object containing 2 dna contigs and 2 overlapping genes composed
of 2 and 3 exons.

Extra data controlling drawing is added to the individual GeneFeature
objects using the addMember method.  It is this which allows the
definition of (1) mouseover text (2) href links (3) drawing as an
intron.  There are currently 2 types of introns: type 1 are drawn as
solid lines/boxes and type 2 are drawn as dashed lines/boxes.  One use
of type 2 drawing is to allow 'introns' to be drawn to join up a gene
structure, but indicate where there is likely to be a missing exon in
the structure.  This control is independent of the intronstyle
parameter which controls the selection of straight lines or angle
links.

Example script

    #!/usr/local/bin/perl

    use strict;
    use GFF::GifGFF;
    use GFF;
    
    # get GFF (normally this would be read from a file, but here it is
    # hard coded just to show you exactly what is needed)
    my $gff = new GFF::GeneFeatureSet;
    &get_gff($gff);
    
    # grouping function (used for bumping, to keep features with same name on same line)
    my $group_by_seqname=sub{
        my $self=shift;
        return $self->seqname;
    };
    
    # label function (specifies what is used for text written to gif)
    my $make_label=sub{
        my $self=shift;
        return $self->getMember('label');
    };
    
    # to map gff features->layout labels
    my $group_features=sub{
        my $self=shift;
        my $label=$self->feature;
        # exon,intron->gene
        if($label=~/^exon|intron$/){
    	return 'gene';
        }
        return $label;
    };
    
    # differential colour function based upon 'Sequence_by' [group] and strand
    my $seq_colour = sub {
        my $coltab = shift ; # a reference to a hash table of colours
        my $gf     = shift ;
	return $coltab->{'black'} if !(defined($gf) and $gf) ;
        my $source = $gf->group_value('Sequenced_by') ;
        my $strand = $gf->strand ;
        if($source =~ /Sanger/) {
	    if($strand eq '+') {
		return $coltab->{'blue'} ;
	    } else {
		return $coltab->{'darkblue'} ;
	    }
	} elsif($source =~ /WUSTL/) {
	    if($strand eq '+') {
		return $coltab->{'green'} ;
	    } else {
		return $coltab->{'darkgreen'} ;
	    }
	} else {
	    if($strand eq '+') {
		return $coltab->{'red'} ;
	    } else {
		return $coltab->{'darkred'} ;
	    }
	}
    };

    # this specifies how the draw each label and their relative positions
    # label names starting with /^label/ have special behaviour (see below)
    # 5th field can be either a colour (for reverse direction) or F/R/P (see doc)
    # (reverse genes are duplicated, to show the effect of this)
    # 6th field can be flag to indicate if feature should be drawn with a box unless 
    # label, where this field should contain a function to generate text to be drawn
    my @layout=(
    	    [-4,15,'label','black','black',$make_label,$group_by_seqname],
    	    [-3,15,'labelg','red','darkred',$group_by_seqname,$group_by_seqname],
    	    [-2,6,'gene','red','darkred',1,$group_by_seqname],
    	    [-1,15],
    	    [0,10,'contig',$seq_colour,$seq_colour],
    	    [1,15],
    	    [2,10,'hom','blue','P',1,$group_by_seqname],
    	    [3,6,'gene','red','R',1,$group_by_seqname],
    	    );

    # this specifies row labels in the left hand margin of the plot
    my %rowlabels = (
		    'gene'   => 'Gene,red,http://www.sanger.ac.uk/HGP/Genes/',
		    'contig' => 'Contig Map',
    ) ;

    my $file1='clone1.html';
    my $file2='clone1.gif';
    open(OUT1,">$file1") || die "cannot open $file1";
    print OUT1 <<ENDOFTEXT;
    <html>
    <head>
      <title>Output from example1.pl</title>
    </head>
    <body bgcolor="#FFFFFF">
    <map name="gifmap">
    ENDOFTEXT
    
    open(OUT2,">$file2") || die "cannot open $file2";
    binmode(OUT2);
    my($x,$y)=&GFF::GifGFF::gff2gif(gff => $gff,
    				layout  => \@layout,
				rowlabels => \%rowlabels,
    				filter  => $group_features,
    				io      => \*OUT2,
    				width   => 600,
				leftmargin => 100,
    				scale   => 200,
    				mio     =>  \*OUT1,
    				intronstyle => 'straight',
    				margin => 10,
    				);
    
    print OUT1 <<ENDOFTEXT;
    </map>
    <h2>Output from example.pl</h2>
    <img src="$file2" usemap="#gifmap" width=$x height=$y border=0>
    <pre>
    ENDOFTEXT
    $gff->dump(\*OUT1);
    print OUT1 <<ENDOFTEXT;
    </pre>
    </body>
    </html>
    ENDOFTEXT
    
    close(OUT1);
    close(OUT2);
    
    sub get_gff{
        my($gff)=@_;
    
        # make an example gff with a contig + a 2 exon gene
    
        # make a contig gf 
        my $gf=new GFF::GeneFeature;
        my $name='contig1';
        $gf->seqname($name);
        $gf->feature('contig');
        my $start=1;
        my $end=74999;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('+');

        # Flag source of sequence, for differential colouring
        $gf->group_value('Sequenced_by',0,'Sanger') ;

        # add a mouseover label
        my $contig1_label="$name:$start-$end";
        $gf->addMember($contig1_label,'label');
        # add a href to another file (dummy)
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make a contig gf 
        my $gf=new GFF::GeneFeature;
        my $name='contig1';
        $gf->seqname($name);
        $gf->feature('contig');
        my $start=75000;
        my $end=123456;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('+');

        # Flag source of sequence, for differential colouring
        $gf->group_value('Sequenced_by',0,'Whitehead') ;

        # add a mouseover label
        my $contig1_label="$name[Whitehead]:$start-$end";
        $gf->addMember($contig1_label,'label');
        # add a href to another file (dummy)
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make a bar (to show that contigs are not linked)
        $gf=new GFF::GeneFeature;
        $name='bar';
        $gf->seqname($name);
        $gf->feature('bar');
        $start=124456;
        $end=124456;
        $gf->start($start);
        $gf->end($end);
        $gff->addGeneFeature($gf);
    
        # another contig gf (reversed)
        $gf=new GFF::GeneFeature;
        $name='contig2';
        $gf->seqname($name);
        $gf->feature('contig');
        $start=125456;
        $end=234567;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('-');

        # Flag source of sequence, for differential colouring
        $gf->group_value('Sequenced_by',0,'Wustl') ;

        # add a mouseover label
        my $contig2_label="$name:$start-$end";
        $gf->addMember($contig2_label,'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make labels for each contig
        $gf=new GFF::GeneFeature;
        $name='contig';
        $gf->seqname($name);
        $gf->feature('label');
        $start=1;
        $end=1;
        $gf->start($start);
        $gf->end($end);
        $gf->addMember($contig1_label,'label');
        $gff->addGeneFeature($gf);
        $gf=new GFF::GeneFeature;
        $name='contig';
        $gf->seqname($name);
        $gf->feature('label');
        $start=125456;
        $end=125456;
        $gf->start($start);
        $gf->end($end);
        $gf->addMember($contig2_label,'label');
        $gff->addGeneFeature($gf);
    
        # make and label some genes
        # label this gene
        $gf=new GFF::GeneFeature;
        $gf->seqname('note');
        $gf->feature('label');
        $start=10000;
        $end=10000;
        $gf->strand('+');
        $gf->start($start);
        $gf->end($end);
        $gf->addMember('[genes overlap]','label');
        $gff->addGeneFeature($gf);
    
        # gene1 (2 exons - features grouped using name)
        $name='gene1';
        # make 1st exon gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('exon');
        $start=100;
        $end=1150;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('+');
        # add a mouseover label
        $gf->addMember("Gene $name, Exon 1:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make intron gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('intron');
        $start=1151;
        $end=10299;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('+');
        # identify as an intron
        $gf->addMember(1,'intron');
        # add a mouseover label
        $gf->addMember("Gene $name, Intron 1:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make 2nd exon gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('exon');
        $start=10300;
        $end=11450;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('+');
        # add a mouseover label
        $gf->addMember("Gene $name, Exon 2:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # label this gene
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('labelg');
        $start=100;
        $end=100;
        $gf->strand('+');
        $gf->start($start);
        $gf->end($end);
        $gff->addGeneFeature($gf);
    
        # gene2 (2 exons)
        $name='gene2';
        # make 2nd exon gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('exon');
        $start=10000;
        $end=10400;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('-');
        # add a mouseover label
        $gf->addMember("Gene $name, Exon 2:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make intron gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('intron');
        $start=10401;
        $end=30000;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('-');
        # identify as an intron
        $gf->addMember(1,'intron');
        # add a mouseover label
        $gf->addMember("Gene $name, Intron 1:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make 1st exon gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('exon');
        $start=30001;
        $end=30300;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('-');
        # add a mouseover label
        $gf->addMember("Gene $name, Exon 1:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make intron gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('intron');
        $start=30001;
        $end=60000;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('-');
        # identify as an intron
        $gf->addMember(2,'intron');
        # add a mouseover label
        $gf->addMember("Gene $name, Intron 2:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make 1st exon gf 
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('exon');
        $start=60001;
        $end=61000;
        $gf->start($start);
        $gf->end($end);
        $gf->strand('-');
        # add a mouseover label
        $gf->addMember("Gene $name, Exon 3:$start-$end",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    
        # make a label
        $gf=new GFF::GeneFeature;
        $gf->seqname($name);
        $gf->feature('labelg');
        $start=10000;
        $end=10000;
        $gf->strand('-');
        $gf->start($start);
        $gf->end($end);
        $gff->addGeneFeature($gf);
    
        # make hom features
    
        $gf=new GFF::HomolGeneFeature;
        $name='hom1';
        $gf->seqname($name);
        $gf->feature('hom');
        $start=10300;
        $end=13450;
        $gf->start($start);
        $gf->end($end);
        $gf->start2(14000);
        $gf->end2(18000);
        $gf->percentid(90);
        $gf->strand('+');
        # add a mouseover label
        $gf->addMember("Homol $name (90%)",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
    
        $gff->addGeneFeature($gf);
        $gf=new GFF::HomolGeneFeature;
        $name='hom2';
        $gf->seqname($name);
        $gf->feature('hom');
        $gf->start(80000);
        $gf->end(90000);
        $gf->start2(93456);
        $gf->end2(123456);
        $gf->percentid(20);
        $gf->strand('+');
        # add a mouseover label
        $gf->addMember("Homol $name (20%)",'label');
        # add a href to another file
        $gf->addMember("$name.html",'href');
        $gff->addGeneFeature($gf);
    }
    
    
There is a special behaviour for any GeneFeature object that is
grouped into a set labelled 'bar', which does not need any entry in
the layout array since it causes a vertical bar at that coordinate
point to be drawn.  Another special drawing mechanism is for any
layout array entry with a label of form /^label/ where text is written
at the coordinate of the start of the feature rather than a box of the
feature itself.

=over 12

=item C<'bar'>

Does not have a layout ruleset.  Causes a vertical black bar to be
drawn at the coordinates specified in the GeneFeature.  A bar is
specified as follows:

    my $gf=new GFF::GeneFeature;
    $gf->seqname('bar');
    $gf->feature('bar');
    $gf->start(123460);
    $gf->end(123500);
    $gff->addGeneFeature($gf);

=item C<'barwidth'>

By default the parameter 'barwidth' is set to 3, which means the C<bar>
elements are 3 pixels wide.  If you wish the have the width defined by
start,end (remember this will depend in terms of the pixels on the
coordinate scaling for that particular gif) then set 'barwidth' to zero.

=item C<'barcolour'>

By default the parameter 'barwidth' is set to 'black'. An alternate 
colour (string label) may be given to 'barcolour'.

=item C<'map_position'>

This feature does not have a layout ruleset nor triggers any drawing. For such GFF,
GifGFF merely computes the corresponding X,Y pixel coordinates of the feature
on the GD canvas and annotates the GeneFeature record accordingly, as the 
'Pixel_Position <x1> <y1> <x2> <y2>' [attribute] tag-value, such that the
<x1> and <x2> correspond to the <start> and <end> pixel coordinates of the feature, 
whereas the <y1> and <y2> are the top and bottom of the map (conceptually equivalent
to a 'bar' feature top and bottom).

Since the script invoking gifgff would know the reference pointer to such features (having
created them), this annotation is thus available upon return from the gifgff method call.

A typical use of this functionality is to permit the calling routine to better align
items subsequently drawn upon the canvas, to map features (previously) drawn by GifGFF.

=item C<'label'>

Does have a layout ruleset. Causes the text returned by the function
passed in the layout definition to be written at the start coordinates
specified in the GeneFeature at the top of the image.  A label might be
specified as follows:

    # layout line
    [-4,15,'label','black','black',$make_label,$group_by_seqname],

    # function for label (take text from seqname)
    my $make_label=sub{
	my $self=shift;
	return $self->seqname;
    };

    my $gf=new GFF::GeneFeature;
    my $text='Title Text';
    $gf->seqname($text);
    $gf->feature('label');
    $gf->start(0);
    $gf->end(0);
    $gff->addGeneFeature($gf);

=back

=head2 Required Parameters

=over 12

=item C<gff>

Supply a gff object as C<gff=$gff>.

=item C<io>

Supply a filehandle for gif file as C<io=$TOUT> where $TOUT is set by:

    my $TOUT;
    open(TOUT,">$html_dir/$clone.gif") || die "cannot open $clone.gif";
    $TOUT=\*TOUT;
    binmode(TOUT);

Image will be written to C<$html_dir/$clone.gif>.

=item C<filter>

Supply a gff filter as C<filter=$group_features> where an example of a
filter is:

    my $group_features=sub{
	my $self=shift;
	return $self->feature;
    };

Features will be displayed according labels derived from the
C<feature> element of each GeneFeature object in the GFF object.

=item C<layout>

Supply a layout array as C<layout=\@layout> where an example of a
layout is:

    my @layout=(
		[-3,5,'label','black','black',$make_label,$group_by_name],
		[-2,5,'gene','red','F',1,$group_by_name],
		[-1,5],
		[0,5,'contig','green','darkgreen'],
		[1,5],
		[2,5,'gene','darkred','R',1,$group_by_name],
		[3,15,'hom','red','P',1,$group_by_name],
		);

Each line describes a set of rules for drawing the set of features
which have been labelled, in this case 'gene' or 'contig', by the
filter function that was passed to the call.  The fields are: (1)
B<real> vertical position of feature line, (2) B<int> height of
feature in pixels, (3) B<text> label identifying this drawing line, (4)
B<text/ref(sub)> colour1 [see below] (5) B<text/ref(sub)> colour2 
[see below for exceptions], (6) no_box_outline/minimum pixel size flag 
[see below for exceptions], (7) bump filter [see
description in discussion of bumping].

Only features labelled C<gene> or C<contig> from output of
group_features filter will be displayed.

Contigs will be drawn in a single y axis region of the gif with green
for a forward contig and darkgreen for a reverse contig (indicated by
the C<strand> element of each GeneFeature object.

Genes will be drawn on two separate y axis regions of the gif
depending on the direction (C<strand>).  Forward will be drawn above
the contig, reverse will be drawn below the contig, in red and darkred
accordingly.

Both Contigs and Genes will be 5 pixels high and will be separated
from each other by 5 pixels.

Features labelled 'gene' will be first grouped according to the
function $group_by_name before being bumped, to allow all the elements
of a gene (exons, introns) to be kept together on a single line.

Contigs will be drawn as a box coloured green with a black outline.
Gene will be coloured in plain red (see that 'no_box_outline' flag is
set).

Field 6 of a @layout row is the 'no_box_outline' box flag, which specifies
whether an rectangular border is drawn around a feature box. It is
important to note that if this flag is set, and your feature widths
scale down to some small pixel value, that the colour of your features
will generally become 'all border' (generally black!).

Note that from version 1.12, this flag value can have any positive value
greater than 1. If this is the case, then feature in the given row are are 
coerced to have that number of minimal pixel width when rendered. This
can be useful for enhancing the print visibility of features which
resolve (after scaling) to only one pixel in thickness. Setting the
'no_box_outline' to 2 ensures that such features will always be at 
least 2 pixels wide.

Functions with a label starting with 'label' have special behaviour
(not used here).  In this case it is assumed that 'text' should be
written to the canvas rather than drawing a box.  The position of the
text is defined by $gf->start.  Obviously field 6 ('no_box_outline')
has no meaning in the context of writing text, so in this case it is
expected to be another filter that will return the text string to be
written.  This allows text to be built on the basis of values in the
$gf or simply the content of $gf->seqname etc. 

Field 5 can take a second colour for reversed features (if forward and
reverse are to be drawn on the same line) or it can be 'F' or 'R' to
indicate that this line should only contain forward or reverse
features.  This allows two alternative displaying modes for example of
forward and reverse features on opposite sides of the DNA or all
grouped together.

In addition, under these circumstances, field '5' can also take the
special directive 'V' (optionally with 'B') which direct that the 
specified line of row of text is drawn vertically, (i.e. at 90 degrees),
in the text direction reading from bottom to top. By default, the label
is written 'top' aligned, i.e. with the gene feature start coordinate 
positioned at the upper left hand corner of the text box (or upper 
right of the text itself).  If the 'B' directive is additionally given,
then the alignment becomes 'lower left hand corner' that is, 'Bottom' 
aligned.  Since the text characters are 6 x 12 pixels (gdGiantFont) 
in size, a feature height of n x 6 pixels (where n == size of the
longest text string, for the given layout row, should ideally 
be provided. The colour of the text is that designated by field 4.
The 'V' directive (with or without 'B') may be given alongside a 'F' 
or 'R' directive, in which case, strand specific labels are drawn vertically.
When 'V' or 'VB' is given, it must be the prefix directive (i.e. it must
*preceed' the 'F' or 'R' directive, i.e. 'VF','VR', 'VBF' or 'VBR').

An additional directive is the 'G' directive, which toggles gridlines off
and on in columns in a sequential manner, when the 'gridlines' argument
is set.  Note that this directive may be combined with 'VB(F|R)' but must
be the last directive (e.g. 'VBFG', not 'GVBF' or similar).

A further exception to Field 5 is when it takes the value 'P' to
indicate shading by percentage homology.  For this value to work, the
corresponding gene features must be of the 'GFF::HomolGeneFeature' class
or have the appropriate Version 2 tag-value tags set, most importantly,
a 'Homology <target> <start2> <end2>' tag value (functionally equivalent
to the ACEDB 'Target' dumped tag, just renamed, to avoid problems
handling 'Target' tagged GFF which is not meant to be drawn as 
percentage homology features, e.g. RepeatMasker repeats). In the case 
of the 'P' colour directive, colour1 (field 4) must be either 'red', 
'green' or 'blue'.  In this case a shape is drawn such that the
coordinates of the top line correspond to the values of 'start','end'
in the object and the bottom line to 'start2','end2'.  The 'percentid'
field is used to set the colour such that 100% homol matches are the
raw colour fading to a light colour as the homology decreases.

The colour argument to field 4 and/or field 5 (subject to the constraints
outlined above) may be a (Perl 'CODE') reference to a user defined colour 
function  expecting up to two arguments: the first argument is a module provided 
reference to a hash of allocated colours ('the palette'), keyed by 
string labels for each (system or user defined) allocated colour, from
which a GD::Image colour handle may be selected and returned to the 
module point calling the function.   Note: such a function may also be
defined for homology 'P' colour1 values, but rather, should return
one of 'red', 'green' or 'blue', *NOT* GD::Image colour handles.
 
The second argument, if present, is a reference to the current 
GFF::GeneFeature being drawn. The function should be designed to always 
return a colour value from the palette, including a default colour in
the absence of a defined Gene Feature.

It is allowed to have multiple drawing rules for the same label, so
complex gifs can be constructed with the same features repeated in different
places but drawn in the same or different ways.

With the exception of the constraints outlined above, the colours available
for use in gifgff C<layout> specifications (or user defined colour functions)
are specified as string scalars taken from the following list:

    white (presumably with a black border!)
    gray
    black
    red
    lightred
    darkred
    green
    lightgreen
    darkgreen
    blue
    lightblue
    darkblue
    yellow
    brown
    magenta (a.k.a. purple)
    darkmagenta
    cyan
    darkcyan

=back

=head2 Arguments Modifying Default Parameter Values

=over 12

=item C<image>

By default, a new blank GD::Image() is created by gff2gif() function;
however, an existing GD::Image() reference handle (i.e. one created
by a 'GD::Image->new($width,$height)' call) can be passed to the function,
in which case the gff2gif() function writes its features onto the
canvas of this existing image. See also C<xorigin> and C<xorigin> below.

=item C<palette>

A reference to a hash of values already colorAllocate()'d by the
caller, into the GD::Image() 'image', and indexed by string labels 
denoting custom colours.

=item C<width>

Overrides the default width (400) of the gif generated, (when a 
user-defined C<image> is not provided).

=item (C<x_origin>,<y_origin>)

Plot origins (h,k) in absolute GIF coordinates. Used to calculated
actual left and top borders, and also, in vertical plot mode,
as a point on the diagonal of reflection of coordinates (i.e. to
transform a horizontal map into a vertical one).

=item C<margin>

Overrides the default margin space (20) for all margins, around the 
objects on the gif image. Note that *all* margins are relative sizes,
to which the specified (x_origin,y_origin) is added to obtain the
true borders of the plot.

See also C<leftmargin>, C<rightmargin>, C<topmargin> and C<bottommargin>.

=item C<leftmargin>

Overrides the default left margin space equal to C<margin>.

=item C<rightmargin>

Overrides the default right margin space equal to C<margin>.

=item C<topmargin>

Overrides the default top margin space equal to C<margin>.

=item C<bottommargin>

Overrides the default bottom margin space equal to C<margin>.

=item C<scale>

Overrides the default scale (500) in bases per pixel of the DNA drawn
on the gif image.

=item C<yseparation>

Overrides the default separation (20) in pixels between the blocks in
the wrapped display.

=item C<expandmap>

Overrides default expansion (2) in pixels of the map box, along the x
axis, for features other than introns (which are contracted by same
amount).  The idea is that this makes it easier to hold the mouse over
a small feature.

=item C<compress>

C<compress=0> turns off (default on) mode to remove space for empty
label types in the gif image.  Make gif smaller.

=item C<bump>

C<bump=0> turns off bumping (default on).  If there are multiple
features (such as gene predictions) that overlap, then normally you
don't want these to overwrite each other.  Bumping fixes this.
Bumping is away from the centre of the display, which is defined as
index 0 in the layout array (such that genes and other features lie
preferentially close to the clone when displayed on either size.

=item C<bump_filter>

By default features are not grouped before being bumped, however there
are situations where you may want grouped bumping, e.g.  if you have
multiple blast hits to a protein, you would like them all on the same
line; all exons and introns of a gene should be on the same line etc.
This can be controlled by defining a filter and adding it as a
parameter to the layout array.  This allows different rules to be
defined for different lines in the layout.

    my $group_by_name=sub{
        my $self=shift;
        return $self->seqname;
    };

=item C<bumpspace>

Allows default spacing (2) in pixels between bumped features to be overriden.

=item C<intronstyle>

Overrides default intronstyle (straight).  Valid intronstyles are
'straight' and 'angle'.

=item C<font>

GD font to use in text labels of the display. Defaults to 'small' (gdSmallFont)
but can also be set to 'tiny' (GD::gdTinyFont) or 'giant'(GD::gdGiantFont)  instead.

=item C<bumplabel>

Turns on bumping for labels in 'label' type layout rows. Default 0 (off).

=item C<bumplabelseparation>

Scaling factor to modulate label row bumping threshhold. Defaults to 1.
Setting to a smaller value allows for more crowded labels prior to 
bumping being triggered.

=item C<resolution>

This affects the resolution of the 'bump' calculation.  The bumping
calculation is implemented by masking a string (see code) and 'resolution'
is a parameter which defines how many bases correspond to one character in
that string.  If resolution=1 then each base would have own character,
which would give the most precision, however if very long (megabase)
fragments were drawn the script would run out of memory.  If resolution is
set too high then features will be bumped when they don't need to be.  The
default of 100 is probably a good compromise for most situations.

=back

=head2 Optional Parameters

=over 12

=item C<mio>

filehandle for html file.  If set, map coordinates are generated and
written out to this filehandle.  <map> and </map> tags are not
written.

=item C<start>

By default, entire range of GFF object is displayed.  If start and end
are set then only that portion of GFF object is displayed.

=item C<end>

See C<start>.

=item C<height>

By default gif created is that has width C<width> and height
determined by the range of the GeneFeatures in the GFF object and
C<scale>.  If you set C<height> then C<scale> will be overridden so
that the scale on the image is that required to fit the GFF (or range
selected by C<start> and C<end>) into the size available.

=item C<transparent>

When defined, background colour (white) is made GIF transparent.

=item C<vertical>

When defined, dictates that the GifGFF map is plotted vertically,
rather than horizontally. To achieve this, the vertical plot is
a diagonally mirrored version of the horizontal plot. The
horizontal left to right scrolling of the map is thus read
from bottom to top, and the top to bottom stacking of rows of the
map is reoriented left to right.  All text in the map is adjusted 
accordingly such that horizontal text becomes vertical and
vertical text becomes horizontal, plotted in the same
location of the map.

=item C<rowlabels>

Setting this parameter triggers the printing of row labels
in the lefthand margin of the plot. The value of the parameter
must be a reference to a hash of such labels C<rowlabels=\%rowlabels> , 
indexed upon the C<@layout> row identification field tag(*). Rows 
without defined rowlabels are left blank. 

Optionally, a specified colour and URL may be appended to the label, 
comma delimited, for example:

   $rowlabel{'GFF'} = 'GFF,red,http://www.sanger.ac.uk/Software/GFF/' ;

If no colour is provided, then the field 4 colour is used. Note, however,
that if the field 4 colour is a reference to a user specified colour mapping
function, that no gene feature is available to the function, hence
only the default colour returned by the function will be used.

If a URL is provided, the row label box become a map coord clickable box.

(*) with strand decoration: namely, where the 'reverse' colour for a given C<@layout>
    identification 'tag' is 'F', 'R', the actual row identification 
    field tag is 'tagF' (or 'tagR').

=item C<gridlines>

Setting this parameter triggers the drawing of light gray background gridlines to 
the GifGFF map. The value of the argument is the number of kilobases between 
gridlines.

=item C<gridcolour>

Colour of C<gridlines> if turned on (default: gray)

=item C<mapA_offset>/C<mapB_offset>

Bulk (basepair) coordinate offsets for feature rows belonging to two independent
coordinate systems. Useful to align a pair of maps for two clones being compared 
(e.g. for comparative human/mouse genomics). Default both to zero.

Note: these offsets are added to the drawing coordinates of the rows, not to the
feature coordinates themselves. Moreover, these offsets are only applied to feature
rows which have a 'A' or 'B' suffix in the feature labels, e.g. contigA or geneB rows.

=back

=head2 Return values

Since width can be a default and height may be calculated by the
routine it is useful to return these values as C<x>, C<y> (e.g. so
they can be written in the <img> tag to speed loading).

=head1 REVISIONS

    1.16 (16/07/2000) rbsk: _colour_init: let user palette also be [RGB] vlaues to be allocated

    1.15 (05/04/2000) rbsk: rather, use 'Homology' tag (which would be a renamed 'Target' tag) for HomolGeneFeature's
			    labelling, to avoid problems with 'Target's which ought not be be treated as 
			    HomolGeneFeatures (like RepeatMasker similarity hits).

    1.14 (22/03/2000) rbsk: treat 'Target' tags as homology information, even if not a HomolGeneFeature

    1.13 (17/03/2000) rbsk: added provisions for two clone coordinates systems ('A' and 'B') in one map, such that 
			    all the rows belonging to one can be offset relative the other (although
			    the displayed coordinates remain relative to the given 'A' or 'B' clone map).

    1.12 (28/02/2000) rbsk: no_box_outline @layout flag generalized to signal 'minimum pixel size' if >1 in size
			    added (optional) map grid lines C<gridlines>, C<gridcolour>
			    Gridlines can also be toggled on and off by the 'G' colour2 directive
			    added 'map_pos' GFF feature type, which 
			    like 'label' and 'bar' has special behaviour (see above)

    1.11 (21/02/2000) rbsk: repaired logical bug in vertical drawing mode

    1.10 (21/11/99): rbsk:
    - 'transparent' flag provide for optional transparency in the GIF background
    - 'vertical' map drawing mode added!!
    - fixed _colour_init() bug (palette not properly allocated if no external palette)

    1.09 (11/11/99): rbsk:
    - colour allocation bug fixed in _colour(): returns undef if undefined, so test for that
    - C<font> argument added to provide flexibility in default base font definition

    1.08 (10/11/99): th/rbsk
    - several enhancements in label bumping, 'bar' drawing, etc.

    1.07 (26/10/99): rbsk
    - 'V' may now be combined with 'R' or 'F' ie. 'VR' or 'VF'
    - added 'B' specifier to direct 'V' alignment

    1.06 (26/10/99): rbsk
    - (bug fix) 'bar' object drawing needed to be offset by topmargin
    - added 'barcolour' argument

    1.05 (18/10/99): rbsk

    - User defined C<palette> of GD::Image() colours.
    - 'label*' text 'V' field 5 directive, for vertically drawn labels

    1.04 (11/10/99): rbsk

    - generalized colour specification mechanism in C<@Layout> fields 4 and 5
      to allow the use of an reference to a user defined feature discriminator 
      function so that colour used in drawing specific features is gene feature 
      content specific (e.g. different clone sources). See example script.

    1.03 (30/9/99): rbsk

     - added more colours to the GifGFF palette of colours: light versions of the primary colours,
       dark versions of the alternate colours, plus 'gray'; 'magenta' considered equivalent to 'purple'

    1.02 (28/9/99): rbsk

     - added arguments: 'image', 'leftmargin', 'rightmargin', 'topmargin', 'bottommargin' and 'rowlabels'

     - because pre-existing GD 'image' handle may now be provided, the 'io' argument may be
       undefined in such situations, leaving the user the responsibility of printing out the Gif image.

     - the 'rowlabels' argument provides a means of labeling,within the left margin space, any row 
       generated by a given @Layout entry. See example script.

    1.01 (07/99):   th - Creation

=cut

