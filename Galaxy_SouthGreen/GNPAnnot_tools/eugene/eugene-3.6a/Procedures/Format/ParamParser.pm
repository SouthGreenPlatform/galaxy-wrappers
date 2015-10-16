# 
# 	Jerome.Gouzy@toulouse.inra.fr
#	Last updated: October 29, 2004
#

package ParamParser;

=pod

=head1 NAME

ParamParser - parse parameters from different sources (CGI.pm, GetOpt, cgi-lib, configuration file, ARGV, ENV)

=head1 SYNOPSIS


	1. parameter source defined from a configuration file
	
	use ParamParser;
	$rh_param = New ParamParser($filename);

		------ example.cfg -------
		# lines starting with # are ignored
		OPTION=value of the option
		--------------------------

	2. from ARGV 
	
	use ParamParser;
	$rh_param = New ParamParser('ARGV');

		% program OPTION1="value of the option" OPTION2=value
 
	3. from environment variables
	 
	use ParamParser;
	$rh_param = New ParamParser('ENV');
	or
	$rh_param = New ParamParser('ENV','prefix'); to add a tag to environment variables

	4. from CGI object
	 
	use CGI;
	use ParamParser;
	$rh_param = New ParamParser('CGIPM');

	5. from CGI-LIB data structure (version 2)
	
	require "cgi-lib2.pl";
	use ParamParser;
	$rh_param = New ParamParser('CGILIB');


	6. from Getopt::Std object
	 
	use Getopt::Std;
	use ParamParser;
	$rh_param = New ParamParser('GETOPTSTD',"list_of_singlet-character_switches"); 

	run the command man Getopt::Std to see what is "list_of_singlet-character_switches"
	to use  the same options with the current module you must write
	$rh_param = New ParamParser('GETOPTSTD',"oif:");
	$rh_param = New ParamParser('GETOPTSTD',"oDI");

	7. from Getopt::Long object
	 
	use Getopt::Long;
	use ParamParser;
	$rh_param = New ParamParser('GETOPTLONG',(list_of_getoptlong_option)); 

	run the command man Getopt::Long to see what is a "list_of_getoptlong_option"
	to use the same options with the current module you must write
	$rh_param = New ParamParser('GETOPTLONG',("length=i","file=s","verbose"));

    8. automatic detection of the source
	$rh_param = New ParamParser('AUTO');

			1. ARGV
			2. CGIPM
			3. CGILIB
			4. ENV

=head1 DESCRIPTION

24-Jun-2004
add two new functions: IsDefined and SetUnlessDefined
26-Oct-2004
creation of Dump methods in order to export paramparser to a file or %ENV
29-Oct-2004
add the capability to define/select a NameSpace : SelectNameSpace
add a fonction to get all keys matching a given pattern : GetKeys
28-Dec-2004
remove a bug added to the Get function in the 29-Oct-2004 release
remove a bug in the __ToFile method where namespace was not used
the automatic detection is now an explicit parameter of the constructor and not the defaut

=cut


use strict;
use warnings;

use Carp;


BEGIN
{
    our $VERSION="0.36"; 
}

=item New

	see SYNOPSIS

=cut

sub New
{
my ($pkg,$source,@a_opt)=@_;

	my $self = {
		__h_opt => {},
		__h_behaviour => {},
		__nb => 0,
		__mode => "",
		__name_space => "",
		__usage => sub {print "\nWarning: something is wrong but the program usage is not yet described\n" },
		__possible_sources => "",
		__last_source => $source
	};

	&SetBehaviour($self,'assert_strict');
	&UnsetBehaviour($self,'assert_empty_file_allowed');
	&UnsetBehaviour($self,'ignore_space');
	&UnsetBehaviour($self,'exit_on_getopt_error');

	&__InitPossibleSources($self);
	&Update($self,$source,'I',@a_opt);

	bless($self,$pkg);
	return($self);

}

=item Update

	$rh_param->Update(source,mode,GetOpt_Std_or_Long_list_of_option);

	source: CGIPM|CGILIB|GetOptStd|GetOptLong|ARGV|$filename|ENV
	mode:
		I: init : clean the data structure first
		A: append mode : preserve the previous value of duplicate keys
		O: overwrite mode : replace the value of a duplicate key

	Update the data structure with a new parameter source.

=cut

sub Update
{
my ($self,$source,$mode,@a_opt)=@_;
my $opt = ( defined($a_opt[0]) ) ? $a_opt[0] :  "";

	$$self{'__mode'}=$mode;

	if ( defined($source) && -e $source && ! -z $source  )
	{
		&__FromFile($self,$source);
		$$self{'__last_source'}="$source";
		return;
	}

	if ( $source =~ /AUTO/i )
	{
		# the module tries to find automaticaly the source of parameter
		# (the source cannot be neither GetOpt* nor filename)
		if ( $$self{'__possible_sources'} =~ /ARGV/ )
		{
			$source="ARGV";
		}
		elsif ( $$self{'__possible_sources'} =~ /CGIPM/ )
		{
			$source="CGIPM";
		}
		elsif (  $$self{'__possible_sources'} =~ /CGILIB/ )
		{
			$source="CGILIB";
		}
		else
		{
			$source="ENV";
		}
	}

	if ( ! defined($source) || $$self{'__possible_sources'} !~ / $source / )
	{
		$$self{'__last_source'}="undef";
		return;
	}

	if ( $source =~ /CGILIB/i )
	{
		my(@a_backup)=@_; # this backup is needed because cgi-lib uses @_ as parameter input source
		&__FromCGILIB($self,@a_backup);
	}
	elsif ( $source =~ /CGIPM/i )
	{
		&__FromCGIPM($self);
	}
	elsif ( $source =~ /ENV/i )
	{
		&__FromENV($self);
	}
	elsif ( $source =~ /GETOPTSTD/i )
	{
		&__FromGetOptStd($self,$opt);
	}
	elsif ( $source =~ /GETOPTLONG/i )
	{
		&__FromGetOptLong($self,@a_opt);
	}
	elsif ( $source =~ /ARGV/i )
	{
		&__FromARGV($self);
	}
	$$self{'__last_source'}="\U$source";
}


=item Dump

	$rh_param->Dump(target[,prefix]);

	source: $filename|ENV
	prefix: add the prefix 'prefix' to %ENV keys

=cut

sub Dump
{
my ($self,$target,@a_opt)=@_;
my $opt = ( defined($a_opt[0]) ) ? $a_opt[0] :  "";

	if ( $$self{'__possible_sources'} =~ / $target / )
	{
		if ( $target =~ /ENV/ )
		{
			&__ToENV($self,$opt);
		}
	}
	else	# the parameter is assumed to be a filename
	{
			&__ToFile($self,$target,$opt);
	}
}

=item SelectNameSpace

	$rh_param->SelectNameSpace('NS');	#  create the namespace NS (in fact a prefix to all parameters)
	$rh_param->SelectNameSpace();		#  select the namespace which contains all parameters 

	Select/Init working NameSpace of parameters

=cut

sub SelectNameSpace
{
my($self,$opt) = @_;

	$opt = "" if ( ! defined($opt) );

	$$self{'__name_space'}=$opt;
}


=item Init

	$rh_param->Init();

	Initialise the data structure

=cut

sub Init
{
my($self) = @_;

	$$self{'__nb'} = 0;
	$$self{'__last_source'}="";
	$$self{'__mode'}  ="";
	$$self{'__usage'} = "";
	foreach my $key (keys (%{$$self{'__h_opt'}}))
	{
		delete($$self{'__h_opt'}{$key});
	}
}

=item Set

	$rh_param->Set($opt,$value);

	Associate a new value to $opt

=cut

sub Set
{
my($self,$opt,$value) = @_;

	$$self{'__last_source'}= "INLINE";
	$$self{'__nb'}++ if ( ! defined($$self{'__h_opt'}{$opt}) );
	my $key = $$self{'__name_space'}.$opt;
	$$self{'__h_opt'}{$key}=$value;
}

=item SetUnlessDefined

	$rh_param->SetUnlessDefined($opt,$value);

	Associate a new value to $opt ONLY if the key is not yet defined

=cut

sub SetUnlessDefined
{
my($self,$opt,$value) = @_;

	if ( ! defined($$self{'__h_opt'}{$opt}) ) 
	{
		$$self{'__last_source'}= "INLINE";
		$$self{'__nb'}++; 
		my $key = $$self{'__name_space'}.$opt;
		$$self{'__h_opt'}{$key}=$value;
	}
}

=item Delete

	$rh_param->Delete($opt);

	Delete the $opt key

=cut

sub Delete
{
my($self,$opt,$value) = @_;

	$$self{'__nb'}--;
	my $key = $$self{'__name_space'}.$opt;
	if ( defined($$self{'__h_opt'}{$key}) )
	{
		delete($$self{'__h_opt'}{$key});
	}
}

=item Get

	$value = $rh_param->Get($opt);

	Return the value of $opt key

=cut

sub Get
{
my($self,$opt) = @_;

	my $key = $$self{'__name_space'}.$opt;
	if ( defined($$self{'__h_opt'}{$key})  )
	{
		return $$self{'__h_opt'}{$key};
	}
	else
	{
		return "";
	}

}

=item GetKeys

	@a_keys = $rh_param->GetKeys(pattern);

	Return a list of parameters matching the given pattern

=cut

sub GetKeys
{
my($self,$pattern) = @_;
my @a_keys=();
my $cpt=0;

	my $ns = $$self{'__name_space'};
	foreach my $key (sort keys (%{$$self{'__h_opt'}}))
	{
		if ( $key =~ /^$ns/ )
		{
			my $nkey = $key;
			$nkey =~ s/^$ns//;
			if ( $nkey =~ /$pattern/ )
			{
				$a_keys[$cpt++]=$nkey;
			}
		}
	}
	return(@a_keys);
}

=item IsDefined

	$value = $rh_param->IsDefined($opt);

	boolean, TRUE if the key is defined

=cut

sub IsDefined
{
my($self,$opt) = @_;

	my $key = $$self{'__name_space'}.$opt;
	my($bool) = ( defined($$self{'__h_opt'}{$key}) ) ? 1 : 0;
	return $bool;
}

=item HowMany

	$value = $rh_param->HowMany();

	Return the number of parameters

=cut

sub HowMany
{
	my($self) = @_;

	return $$self{'__nb'};
}

=item GetSource

	$value = $rh_param->GetSource();

	Return the last parameter source

=cut

sub GetSource
{
my($self) = @_;

	return $$self{'_last_source'};
}

=item Print

	$rh_param->Print();
	$rh_param->Print('html');

	Print keys and associated values in text of html format

=cut

sub Print
{
my($self,$format) = @_;

	my($header) = "";
	my($tail) = "";
	my($sep) = ":";
	my($newline) = "\n";
	my($style) = "";
	if ( defined($format) && $format =~ /html/i ) 
	{
		$header = "<table>";
		$tail = "</table>";
		$sep = "<td>";
		$newline = "<tr><td>";
		$style = "<b>";
	}
	print "$header";
	foreach my $key (sort keys (%{$$self{'__h_opt'}}))
	{
		my $ns = $$self{'__name_space'};
		next if ( $key !~ /^$ns/ );
		if ( defined($key) && defined($$self{'__h_opt'}{$key}) )
		{
			print "$newline$style$key$sep ".$$self{'__h_opt'}{$key};
		}
	}
	print "${newline}Total number (all namespaces) of keys$sep ".$$self{'__nb'};
	print "${newline}Last source$sep ".$$self{'__last_source'};
	print "$tail";
}

=item SetBehaviour

	$rh_param->SetBehaviour('assert_strict'); # when set, the assertion will fail if the parameter is not defined (default)
	$rh_param->SetBehaviour('ignore_space');  # when set, the space between the '=' are ignored in the configuration file
	$rh_param->SetBehaviour('exit_on_getopt_error') # execute the usage function when GetOptions return an error code;
	$rh_param->SetBehaviour('assert_empty_file_allowed') # when set, no exit on empty files

	Control the behaviour of the parser

=cut

sub SetBehaviour
{
my($self,$key) = @_;

	$$self{'__h_behaviour'}{$key}= 1;
}

=item UnsetBehaviour

	$rh_param->UnsetBehaviour('assert_strict'); 		 # if unset, the assertion is true when the parameter is not defined
	$rh_param->UnsetBehaviour('ignore_space');  		 # if unset, the space between the '=' are not ignored in the configuration file (default)
	$rh_param->UnSetBehaviour('exit_on_getopt_error')    # ignore the value returned by GetOptions (default)
	$rh_param->UnSetBehaviour('assert_empty_file_allowed') # if unset, then the program exits on empty files (default)

	Control the behaviour of the parser

=cut

sub UnsetBehaviour
{
my($self,$key) = @_;

	$$self{'__h_behaviour'}{$key}= 0;
}

=item AssertFileExists

	$rh_param->AssertFileExists(@a_opt);

	The programs stops if the key $opt does not refer to a non empty file

=cut

sub AssertFileExists
{
my($self,@a_file)=@_;

	foreach my $file (@a_file)
	{
		my $key = $$self{'__name_space'}.$file;
		my($lfile)=$$self{'__h_opt'}{$key};
		next if( ! defined($lfile) && ! $$self{'__h_behaviour'}{'assert_strict'} );
		if ( ! defined($lfile) || ! -e $lfile  || ( -z $lfile &&  ! $$self{'__h_behaviour'}{'assert_empty_file_allowed'} ) )
		{
			&__PrintUsage($self);
			$lfile = &__DefinedIfNot($lfile);
			&Carp::croak("\n=>The value of the parameter $file is >$lfile< which is not a name of an existing and non empty file");
		}
	}

	return(1);
}

=item AssertDirExists

	$rh_param->AssertDirExists(@a_opt);

	The programs stops if the key $opt does not refer to a directory

=cut

sub AssertDirExists
{
my($self,@a_file)=@_;

	foreach my $file (@a_file)
	{
		my $key = $$self{'__name_space'}.$file;
		my($lfile)=$$self{'__h_opt'}{$key};
		next if( ! defined($lfile) && ! $$self{'__h_behaviour'}{'assert_strict'} );
		if ( ! defined($lfile) || ! -d $lfile )
		{
			&__PrintUsage($self);
			$lfile = &__DefinedIfNot($lfile);
			&Carp::croak("\n=>The value of the parameter $file is >$lfile< which is not a name of an existing directory");
		}
	}

	return(1);
}


=item AssertInteger

	$rh_param->AssertInteger(@a_opt);

	The programs stops if one of the key in the list does not refer to an integer

=cut

sub AssertInteger
{
my($self,@a_opt)=@_;

	foreach my $opt (@a_opt)
	{
		my $key = $$self{'__name_space'}.$opt;
		my($lopt)=$$self{'__h_opt'}{$key};
		next if( ! defined($lopt) &&  ! $$self{'__h_behaviour'}{'assert_strict'} );
		if ( ! defined($lopt) || $lopt !~ /^[\+\-]*\d+$/ )
		{
			&__PrintUsage($self);
			$lopt = &__DefinedIfNot($lopt);
			&Carp::croak("\n=>The value of the parameter $opt is >$lopt< which is not a valid integer value");
		}
	}
	return(1);
}

=item AssertDefined

	$rh_param->AssertDefined(@a_opt);

	The programs stop if one of the key in the list is not defined

=cut

sub AssertDefined
{
my($self,@a_opt)=@_;

	foreach my $opt (@a_opt)
	{
		my $key = $$self{'__name_space'}.$opt;
		my($lopt)=$$self{'__h_opt'}{$key};
		if ( ! defined($lopt) )
		{
			&__PrintUsage($self);
			&Carp::croak("=>The parameter $opt must be provided");
		}
	}
	return(1);
}

=item AssertAllowedValue

	$rh_param->AssertAllowedValue($value,@a_list_of_allowed_values);

	The program stop if the value of the key does not match one value of the list of allowed values

=cut

sub AssertAllowedValue
{
my($self,$value,@a_list_of_allowed_values)=@_;

	my $key = $$self{'__name_space'}.$value;
	my($lvalue) = $$self{'__h_opt'}{$key};
	if(  defined($lvalue) )
	{
		foreach my $one_value (@a_list_of_allowed_values)
		{
			if ( $lvalue =~ /^$one_value$/ )
			{
				return(1);
			}
		}
	}
	&__PrintUsage($self);
	my($allowed)=join(',',@a_list_of_allowed_values);
	$lvalue = &__DefinedIfNot($lvalue);
	&Carp::croak("=>The current value of the parameter $value is >$lvalue< which is not in the set of allowed values [$allowed]");
}

=item AssertNonEmptyFile

	$rh_param->AssertNonEmptyFile(@a_opt);

	The programs stops if the elements of the list does not refer to non empty files

=cut

sub AssertNonEmptyFile
{
my($self,@a_file)=@_;

	foreach my $file (@a_file)
	{
		my $file = $$self{'__name_space'}.$file;
		if ( ! defined($file) || ! -e $file || -z $file )
		{
			&__PrintUsage($self);
			$file = &__DefinedIfNot($file);
			&Carp::croak("AssertNonEmptyFile failed for $file");
		}
	}

	return(1);
}

=item Usage

	$rh_param->Usage();
	$rh_param->Usage('html');

	Print the usage of the program

=cut

sub Usage
{
my($self,$format)=@_;
my($head)="";
my($tail)="";

	if ( defined($format) && $format =~ /html/i ) 
	{
		$head="<html><head><title>$0</title></head><body><br><pre>";
		$tail="<br></pre></body></html>";
	}
	print $head;
	&__PrintUsage($self);
	print $tail;
	exit;
}

=item SetUsage

	$rh_param->SetUsage(my $usage= sub { &my_usage_fct();} )

	Attach an usage fonction to the ParamParser object

=cut

sub SetUsage
{
my($self,$r_fct_usage)=@_;

	$$self{'__usage'} = $r_fct_usage;	
	if ( defined($$self{'__h_opt'}{'help'}) || defined($$self{'__h_opt'}{'HELP'}) )
	{
		if ( $$self{'__last_source'} =~ /CGI/i )
		{
			&Usage($self,'html');
		}
		else
		{
			&Usage($self);
		}
	}
}

=head1 EXAMPLE1

	use CGI qw/:standard/;
	use ParamParser;

	my $rh_param =  New ParamParser("CGIPM");

	$rh_param->SetUsage(my $usage=sub { print "\nPlease read the documentation\n"; } ); # attach an usage fonction to the parser
	# the best way is to reference a real fonction $rh_param->SetUsage(my $usage=sub { &UsageFct(); } );
	$rh_param->Set('TIMEOUT','10000');  # add a single variable to the data structure
	$rh_param->Update('ENV',"O");	   # append all environment variables in overwrite mode (overwrite duplicates)
	$rh_param->AssertFileExists('CFG'); # check that the value of the parameter CFG is an existing file, print the usage and exit if it is not.
	$rh_param->Update($rh_param->Get('CFG'),"A");	   # add all variables contained in the configuration file in append mode (do not overwrite duplicates)
	print header;
	$rh_param->Print('html');

=cut

=head1 EXAMPLE2

	use Getopt::Long;
	use ParamParser;

	my $rh_param =  New ParamParser('GETOPTLONG',("help:s","min=i","max=i","inputfile=s","what=s"));

	$rh_param->SetUsage(my $usage=sub { print "\nPlease read the documentation\n"; } ); # attach an usage fonction to the parser
	# the best way is to reference a real fonction $rh_param->SetUsage(my $usage=sub { &UsageFct(); } );

	$rh_param->Update('ENV',"A");	   # append all environment variables in append mode (do not overwrite duplicates)

	$rh_param->AssertFileExists('inputfile'); 		# check that the value of the parameter inputfile is an existing file, print the usage and exit if it is not.
	$rh_param->AssertInteger('max','min');  # check that the value of the parameters are integers, print the usage and exit if one of them is not.
	$rh_param->AssertAllowedValue('what','yes','no','maybe');  # check that the value of the parameters is a correct value
	$rh_param->Print();

=cut


=head1 INTERNAL METHOD CALLS

=cut

=item __PrintUsage

	Print the usage of the program

=cut

sub __PrintUsage
{
my($self)=@_;

	&{$$self{'__usage'}}();
}

=item __UpdateIfPossible

	Update the value of the given key, depending on the selected insertion mode

=cut

sub __UpdateIfPossible
{
my($self,$item,$value)=@_;
	

	my $how = ( $$self{'__mode'} eq "" ) ? "A" : $$self{'__mode'} ;

	$item = $$self{'__name_space'}.$item;
	if (     ! defined($$self{'__h_opt'}{$item})			# the key doesn't already exist 
		||  (defined($$self{'__h_opt'}{$item}) && $how eq 'O') ) # or the key already exists but the mode is 'O'verwrite
	{
		$$self{'__nb'}++;
		$$self{'__h_opt'}{$item}=$value;
	}
}

=item __FromGetOptStd

	Initialize the ParamParser object using Getopt::Std style as source of param/values

=cut

sub __FromGetOptStd
{
my($self,$source,$optlist)=@_;

	use Getopt::Std;
	my @a_backup=@ARGV;

	our %options=();
	&getopts($optlist,\%options);
	#my $getopt_succeed = &getopts($optlist,\%options);
	#if ( ! $getopt_succeed && $$self{'__h_behaviour'}{'exit_on_getopt_error'} )
	#{
	#	&Usage();
	#}
	foreach my $key (keys(%options))
	{
		&__UpdateIfPossible($self,$key,$options{$key});
	}

	@ARGV = @a_backup; # restore original parameters 
			   #	-> can be parsed again is necessary
			   #	-> avoid site effect
}

=item __FromGetOptLong

	Initialize the ParamParser object using Getopt::Long style as source of param/values

=cut

sub __FromGetOptLong
{
my($self,@a_opt)=@_;

	use Getopt::Long;
	my @a_backup=@ARGV;
	my %h_options=();
	my %h_value=();

	foreach my $key (@a_opt)
	{
		$h_options{$key}=\$h_value{$key};
	}
	my $getopt_succeed = &GetOptions(%h_options);

	if ( ! $getopt_succeed && $$self{'__h_behaviour'}{'exit_on_getopt_error'} )
	{
		&Usage($self);
	}

	foreach my $key (keys(%h_value))
	{
		my(@F)=split(/[:=]/,$key);
		my($real_key)=$F[0]; 
		&__UpdateIfPossible($self,$real_key,$h_value{$key});
	}

	@ARGV = @a_backup; # restore original parameters 
			   #	-> can be parsed again is necessary
			   #	-> avoid site effect 
}


=item __FromCGILIB

	Initialize the ParamParser object using CGI-LIB2 as source of param/value

=cut

sub __FromCGILIB
{
my($self,@a_backup)=@_;

	@_=@a_backup; 
	my($keyin);

	if ( defined(ref(&main::ReadParse)) )
	{
		&main::ReadParse;

		foreach $keyin (keys(%main::in))
		{
			&__UpdateIfPossible($self,$keyin,$main::in{$keyin});
		}
	}
}

=item __FromCGIPM

	Initialize the ParamParser object using CGI.pm source

=cut

sub __FromCGIPM
{
	my($self)=@_;
	my($keyin);
	my($cgi)=new CGI;

	my (@a_all_params)=$cgi->param();

	foreach $keyin (@a_all_params)
	{
		&__UpdateIfPossible($self,$keyin,$cgi->param($keyin));
	}
}

=item __FromFile

	Initialize the ParamParser object using a configuration file.

=cut

sub __FromFile
{
my($self,$source)=@_;

	my($lign)="";

	open(PARAMPARSERCFG,"$source") or &Carp::croak("ERROR >$source<\n");
	while($lign=<PARAMPARSERCFG>)
	{
		next if ( $lign =~ /^#/ );
		chop($lign);
		my(@F);
		if ( $$self{'__h_behaviour'}{'ignore_space'} )
		{
			@F=split(/\s*=\s*/,$lign,2);
		}
		else
		{
			@F=split('=',$lign,2);
		}
		next if ( ! defined($F[0]) || ! defined($F[1]) );
		&__UpdateIfPossible($self,$F[0],$F[1]);
	}
	close(PARAMPARSERCFG);
}

=item __FromARGV

	Initialize the ParamParser object using the @ARGV array as source of param/values 


=cut

sub __FromARGV
{
my($self)=@_;

	foreach my $option (@ARGV)
	{
		my(@F)=split('=',$option,2);
		next if ( ! defined($F[0]) || ! defined($F[1]) );
		&__UpdateIfPossible($self,$F[0],$F[1]);
	}
}

=item __FromENV

	Initialize the ParamParser object using the %ENV hash as source of param/values 

=cut

sub __FromENV
{
my($self)=@_;

	foreach my $option (keys(%ENV))
	{
		next if ( ! defined($option) || ! defined($ENV{$option}) );
		&__UpdateIfPossible($self,$option,$ENV{$option});
	}
}

=item __ToFile

	Dump the paramparser into a file

=cut

sub __ToFile
{
my($self,$target,$prefix)=@_;
my $ns = $$self{'__name_space'};

    open(PARAMPARSERCFG,">$target") or &Carp::croak("ERROR >$target<\n");
	foreach my $key (sort keys (%{$$self{'__h_opt'}}))
	{
			if ( defined($key) && defined($$self{'__h_opt'}{$key}) && $key =~ /^$ns/ )
			{
					if ( $prefix ne "" && $key !~ /^$prefix/ )
					{
							my $nkey="$prefix$key";
							print PARAMPARSERCFG "$nkey=".$$self{'__h_opt'}{$key}."\n";
					}
					else
					{
							print PARAMPARSERCFG "$key=".$$self{'__h_opt'}{$key}."\n";
					}
			}
	}
	close(PARAMPARSERCFG);
}

=item __ToENV

	Dump the paramparser into a file

=cut

sub __ToENV
{
my($self,$prefix)=@_;
my $ns = $$self{'__name_space'};

	foreach my $key (sort keys (%{$$self{'__h_opt'}}))
	{
		next if ( $key !~ /^$ns/ );
		if ( defined($key) && defined($$self{'__h_opt'}{$key}) )
		{
			if ( $prefix ne "" && $key !~ /^$prefix/ )
			{
				my $nkey="$prefix$key";
				$ENV{$nkey}="$$self{'__h_opt'}{$key}";
			}
			else
			{
				$ENV{$key}="$$self{'__h_opt'}{$key}";
			}
		}
	}
}



=item __DefinedIfNot

	Init a variable if it is not defined (in order to avoir warnings)

=cut 

sub __DefinedIfNot
{
my($var)=@_;

	if ( !defined($var) || $var eq "" )
	{
		return "undef";
	}
	return$var;
}

=item __InitPossibleSources

	Build a list of possible sources depending on loaded modules

=cut

sub __InitPossibleSources
{
my($self)=@_;
my(%h_src) = (  "CGIPM" =>  defined($CGI::VERSION) ,
				"GETOPTSTD" =>  defined($Getopt::Std::VERSION) ,
				"GETOPTLONG" =>  defined($Getopt::Long::VERSION) ,
				"CGILIB" =>     defined($cgi_lib'version),
				"ARGV" => defined($ARGV[0])
				) ;

	$$self{'__possible_sources'} = " ENV "; 

	foreach my $key (keys(%h_src))
	{
		if ( $h_src{$key} ) 
		{
			$$self{'__possible_sources'} .= " $key ";
		}
	}
}


1;

