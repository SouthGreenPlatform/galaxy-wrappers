package Net::Domain::TLD;
use strict;
use base qw(Exporter);
use 5.006;
our @EXPORT_OK = qw(tlds tld_exists);
our $VERSION = '1.68';

use warnings;
use Carp;
use Storable qw ( dclone );

use constant TLD_TYPES => qw ( new_open new_restricted gtld_open gtld_restricted cc );

=head1 NAME

	Net::Domain::TLD - Work with TLD names 

=head1 SYNOPSIS

	use Net::Domain::TLD qw(tlds tld_exists);
	my @ccTLDs = tlds('cc');
	print "TLD ok\n" if tld_exists('ac','cc');

=head1 DESCRIPTION

	The purpose of this module is to provide user with current list of 
	available top level domain names including new ICANN additions and ccTLDs
	Currently TLD definitions have been acquired from the following sources:

	http://www.icann.org/tlds/
	http://www.dnso.org/constituency/gtld/gtld.html
	http://www.iana.org/cctld/cctld-whois.htm

=cut

my %tld_profile = (
	new_open => { 
		info => q{Unrestricted use},
	},
	new_restricted => { 
		aero => q{Air-transport industry},
		asia => q{Companies, organisations and individuals in the Asia-Pacific region},
		arpa => q{Address and Routing Parameter Area},
		biz => q{Businesses},
		cat => q{Catalan linguistic and cultural community},
		coop => q{Cooperatives},
		jobs => q{Human Resource Management},
		mobi => q{Mobile},
		museum => q{Museums},
		name => q{For registration by individuals},
		pro => q{Accountants, lawyers, and physicians},
		travel => q{Travel industry},
		tel => q{For businesses and individuals to publish contact data}
	},
	gtld_open => {
		com => q{Commercial organization},
		net => q{Network connection services provider},
		org => q{Non-profit organizations and industry standard groups}
	},
	gtld_restricted => {
		gov => q{United States Government},
		mil => q{United States Military},
		edu => q{Educational institution},
		int => q{International treaties/databases},
	},
	cc => {
		ac => q{Ascension Island},
		ad => q{Andorra},
		ae => q{United Arab Emirates},
		af => q{Afghanistan},
		ag => q{Antigua and Barbuda},
		ai => q{Anguilla},
		al => q{Albania},
		am => q{Armenia},
		an => q{Netherlands Antilles},
		ao => q{Angola},
		aq => q{Antartica},
		ar => q{Argentina},
		as => q{American Samoa},
		at => q{Austria},
		au => q{Australia},
		aw => q{Aruba},
	        ax => q(Aland Islands),
		az => q{Azerbaijan},
		ba => q{Bosnia and Herzegovina},
		bb => q{Barbados},
		bd => q{Bangladesh},
		be => q{Belgium},
		bf => q{Burkina Faso},
		bg => q{Bulgaria},
		bh => q{Bahrain},
		bi => q{Burundi},
		bj => q{Benin},
	        bl => q(Saint Barthelemy),
		bm => q{Bermuda},
		bn => q{Brunei Darussalam},
		bo => q{Bolivia},
		br => q{Brazil},
		bs => q{Bahamas},
		bt => q{Bhutan},
		bv => q{Bouvet Island},
		bw => q{Botswana},
		by => q{Belarus},
		bz => q{Belize},
		ca => q{Canada},
		cc => q{Cocos (Keeling) Islands},
		cd => q{Congo, Democratic Republic of the},
		cf => q{Central African Republic},
		cg => q{Congo, Republic of},
		ch => q{Switzerland},
		ci => q{Cote d'Ivoire},
		ck => q{Cook Islands},
		cl => q{Chile},
		cm => q{Cameroon},
		cn => q{China},
		co => q{Colombia},
		cr => q{Costa Rica},
		cu => q{Cuba},
		cv => q{Cap Verde},
		cx => q{Christmas Island},
		cy => q{Cyprus},
		cz => q{Czech Republic},
		de => q{Germany},
		dj => q{Djibouti},
		dk => q{Denmark},
		dm => q{Dominica},
		do => q{Dominican Republic},
		dz => q{Algeria},
		ec => q{Ecuador},
		ee => q{Estonia},
		eg => q{Egypt},
		eh => q{Western Sahara},
		er => q{Eritrea},
		es => q{Spain},
		et => q{Ethiopia},
		eu => q{European Union},
		fi => q{Finland},
		fj => q{Fiji},
		fk => q{Falkland Islands (Malvina)},
		fm => q{Micronesia, Federal State of},
		fo => q{Faroe Islands},
		fr => q{France},
		ga => q{Gabon},
		gb => q{United Kingdom},
		gd => q{Grenada},
		ge => q{Georgia},
		gf => q{French Guiana},
		gg => q{Guernsey},
		gh => q{Ghana},
		gi => q{Gibraltar},
		gl => q{Greenland},
		gm => q{Gambia},
		gn => q{Guinea},
		gp => q{Guadeloupe},
		gq => q{Equatorial Guinea},
		gr => q{Greece},
		gs => q{South Georgia and the South Sandwich Islands},
		gt => q{Guatemala},
		gu => q{Guam},
		gw => q{Guinea-Bissau},
		gy => q{Guyana},
		hk => q{Hong Kong},
		hm => q{Heard and McDonald Islands},
		hn => q{Honduras},
		hr => q{Croatia/Hrvatska},
		ht => q{Haiti},
		hu => q{Hungary},
		id => q{Indonesia},
		ie => q{Ireland},
		il => q{Israel},
		im => q{Isle of Man},
		in => q{India},
		io => q{British Indian Ocean Territory},
		iq => q{Iraq},
		ir => q{Iran (Islamic Republic of)},
		is => q{Iceland},
		it => q{Italy},
		je => q{Jersey},
		jm => q{Jamaica},
		jo => q{Jordan},
		jp => q{Japan},
		ke => q{Kenya},
		kg => q{Kyrgyzstan},
		kh => q{Cambodia},
		ki => q{Kiribati},
		km => q{Comoros},
		kn => q{Saint Kitts and Nevis},
		kp => q{Korea, Democratic People's Republic},
		kr => q{Korea, Republic of},
		kw => q{Kuwait},
		ky => q{Cayman Islands},
		kz => q{Kazakhstan},
		la => q{Lao People's Democratic Republic},
		lb => q{Lebanon},
		lc => q{Saint Lucia},
		li => q{Liechtenstein},
		lk => q{Sri Lanka},
		lr => q{Liberia},
		ls => q{Lesotho},
		lt => q{Lithuania},
		lu => q{Luxembourg},
		lv => q{Latvia},
		ly => q{Libyan Arab Jamahiriya},
		ma => q{Morocco},
		mc => q{Monaco},
		md => q{Moldova, Republic of},
	        me => q(Montenegro),
	        mf => q{Saint Martin (French part)},
		mg => q{Madagascar},
		mh => q{Marshall Islands},
		mk => q{Macedonia, Former Yugoslav Republic},
		ml => q{Mali},
		mm => q{Myanmar},
		mn => q{Mongolia},
		mo => q{Macau},
		mp => q{Northern Mariana Islands},
		mq => q{Martinique},
		mr => q{Mauritania},
		ms => q{Montserrat},
		mt => q{Malta},
		mu => q{Mauritius},
		mv => q{Maldives},
		mw => q{Malawi},
		mx => q{Mexico},
		my => q{Malaysia},
		mz => q{Mozambique},
		na => q{Namibia},
		nc => q{New Caledonia},
		ne => q{Niger},
		nf => q{Norfolk Island},
		ng => q{Nigeria},
		ni => q{Nicaragua},
		nl => q{Netherlands},
		no => q{Norway},
		np => q{Nepal},
		nr => q{Nauru},
		nu => q{Niue},
		nz => q{New Zealand},
		om => q{Oman},
		pa => q{Panama},
		pe => q{Peru},
		pf => q{French Polynesia},
		pg => q{Papua New Guinea},
		ph => q{Philippines},
		pk => q{Pakistan},
		pl => q{Poland},
		pm => q{St. Pierre and Miquelon},
		pn => q{Pitcairn Island},
		pr => q{Puerto Rico},
		ps => q{Palestinian Territories},
		pt => q{Portugal},
		pw => q{Palau},
		py => q{Paraguay},
		qa => q{Qatar},
		re => q{Reunion Island},
		ro => q{Romania},
	        rs => q(Serbia),
		ru => q{Russian Federation},
		rw => q{Rwanda},
		sa => q{Saudi Arabia},
		sb => q{Solomon Islands},
		sc => q{Seychelles},
		sd => q{Sudan},
		se => q{Sweden},
		sg => q{Singapore},
		sh => q{St. Helena},
		si => q{Slovenia},
		sj => q{Svalbard and Jan Mayen Islands},
		sk => q{Slovak Republic},
		sl => q{Sierra Leone},
		sm => q{San Marino},
		sn => q{Senegal},
		so => q{Somalia},
		sr => q{Suriname},
		st => q{Sao Tome and Principe},
		su => q{Soviet Union},
		sv => q{El Salvador},
		sy => q{Syrian Arab Republic},
		sz => q{Swaziland},
		tc => q{Turks and Caicos Islands},
		td => q{Chad},
		tf => q{French Southern Territories},
		tg => q{Togo},
		th => q{Thailand},
		tj => q{Tajikistan},
		tk => q{Tokelau},
		tl => q{Timor-Leste},
		tm => q{Turkmenistan},
		tn => q{Tunisia},
		to => q{Tonga},
		tp => q{East Timor},
		tr => q{Turkey},
		tt => q{Trinidad and Tobago},
		tv => q{Tuvalu},
		tw => q{Taiwan},
		tz => q{Tanzania},
		ua => q{Ukraine},
		ug => q{Uganda},
		uk => q{United Kingdom},
		um => q{US Minor Outlying Islands},
		us => q{United States},
		uy => q{Uruguay},
		uz => q{Uzbekistan},
		va => q{Holy See (City Vatican State)},
		vc => q{Saint Vincent and the Grenadines},
		ve => q{Venezuela},
		vg => q{Virgin Islands (British)},
		vi => q{Virgin Islands (USA)},
		vn => q{Vietnam},
		vu => q{Vanuatu},
		wf => q{Wallis and Futuna Islands},
		ws => q{Western Samoa},
		ye => q{Yemen},
		yt => q{Mayotte},
		yu => q{Yugoslavia},
		za => q{South Africa},
		zm => q{Zambia},
		zw => q{Zimbabwe}
	}
);

my $flat_profile = flatten ( \%tld_profile );

sub flatten {
	my $hashref = shift;
	my %results;
	@results{ keys %{ $hashref->{$_} } } = values % { $hashref->{$_} }
		for ( keys %$hashref );
	return \%results;
}

sub check_type {
	my $type = shift;
	croak "unknown TLD type: $type" unless grep { $type eq $_ } TLD_TYPES;
	return 1;
}

=head1 PUBLIC METHODS

	Each public function/method is described here.
	These are how you should interact with this module.

=head3 C<< tlds >>

	This routine returns the tlds requested.

	my @all_tlds = tlds; #array of tlds
	my $all_tlds = tlds; #hashref of tlds and their descriptions

	my @cc_tlds = tlds('cc'); #array of just 'cc' type tlds
	my $cc_tlds = tlds('cc'); #hashref of just 'cc' type tlds and their descriptions

	Valid types are:
		cc                 - country code domains
		gtld_open          - generic domains that anyone can register
		gtld_restricted    - generic restricted registration domains
		new_open           - recently added generic domains
		new_restricted     - new restricted registration domains

=cut

sub tlds {
	my $type = shift;
	check_type ( $type ) if $type;
	my $results = $type ? 
		wantarray ? [ keys %{ $tld_profile{$type} } ] : 
			dclone ( $tld_profile{$type} ) :
		wantarray ? [ map { keys %$_ } values %tld_profile ] : 
			$flat_profile;
	return wantarray ? @$results : $results;
}

=head3 C<< tld_exists >>

	This routine returns true if the given domain exists and false otherwise.

	die "no such domain" unless tld_exists($tld); #call without tld type 
	die "no such domain" unless tld_exists($tld, 'new_open'); #call with tld type

=cut

sub tld_exists {
	my ( $tld, $type )  = ( lc ( $_[0] ), $_[1] );
	check_type ( $type ) if $type;
	my $result = $type ? 
		$tld_profile{$type}{$tld} ? 1 : 0 :
		$flat_profile->{$tld} ? 1 : 0;
	return $result;
}

=head1 COPYRIGHT

	Copyright (c) 2003-2005 Alex Pavlovic, all rights reserved.  This program
	is free software; you can redistribute it and/or modify it under the same terms
	as Perl itself.

=head1 AUTHORS

	Alexander Pavlovic <alex.pavlovic@taskforce-1.com>
	Ricardo SIGNES <rjbs@cpan.org>

=cut

1;
