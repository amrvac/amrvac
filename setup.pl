#!/usr/bin/perl -i

use strict;
use warnings;
use Getopt::Long;

my $help_message =
"Usage: setup.pl [options]      To generate a new MPI-AMRVAC setup
       setup.pl -show          To show the current options

Options:

    -d=NM                       N is the the problem dimension (1 to 3)
                                M is the vector dimension (1 to 3)
    -p=<physics module>         Which physics module to use (rho, mhd, ...)
    -phi={1,2,3}                Index of vector phi-component (default: 3)
    -z={1,2,3}                  Index of vector z-component (default: 2)
    -eos=<equation of state>    The equation of state
    -nf=<number>                The number of fluid tracers (default: 0)
    -ndust=<number>             The number of dust species (default: 0)
    -arch=<name>                Use compilation flags from arch/<name>.defs

    -show                       Show current options
    -help                       Show this help message

Examples:

setup.pl -d=22 -p=mhd -arch=default
setup.pl -show\n";

# Locally define the variables that will hold the options
my $ndim;
my $ndir;
my $physics;
my $eos;
my $phi_dir;
my $z_dir;
my $arch;
my $ntracers;
my $ndust;
my $show;
my $help;

# Parse the options. Some are handled with a subroutine, which can do immediate
# error checking.
GetOptions(
    "d=i"     =>
    sub {
        my ($opt_name, $opt_value) = @_;
        $opt_value =~ /([123])([123])/ || die "-$opt_name flag incorrect\n";
        $ndim = $1;
        $ndir = $2;

        $ndim >= 1 && $ndim <= 3 ||
            die("1 <= ndim <= 3 does not hold\n");
        $ndim <= $ndir && $ndir <= 3 ||
            die("ndim <= ndir <= 3 does not hold\n");
    },
    "p=s"     => \$physics,
    "eos=s"   => \$eos,
    "phi=i"   => \$phi_dir,
    "z=i"     => \$z_dir,
    "arch=s"  => \$arch,
    "nf=i"    => \$ntracers,
    "ndust=i" => \$ndust,
    "show"    => \$show,
    "help"    => \$help)
    or die("Error in command line arguments\n");


# Show help if -help is given or if there are no other arguments
if ($help || !($ndim || $ndir || $show || $physics || $eos ||
               length($phi_dir) || length($z_dir) || $arch ||
               length($ntracers) || length($ndust))) {
    print STDERR $help_message;
    exit;
}

# Check if the environment variable AMRVAC_DIR is defined
if (!$ENV{AMRVAC_DIR}) {
    print STDERR "Error: \$AMRVAC_DIR variable undefined, set it with:\n";
    print STDERR "export \$AMRVAC_DIR=your/amrvac/dir\n";
    exit;
}

# Show parameters and quit when -s was given
if ($show) {
    show_current_parameters();
    exit;
}

# Get these files if they do not exist already
copy_if_not_present("makefile", "arch", "make_temp");
copy_if_not_present("definitions.h", "src");
copy_if_not_present('mod_indices.t', "src");
copy_if_not_present("definitions.h", "src");

if ($ndim) {
    replace_regexp_file("makefile", qr/ndim\s*=.*/, "ndim = $ndim");
}

if ($ndir) {
    replace_regexp_file("makefile", qr/ndir\s*=.*/, "ndir = $ndir");
}

if (length($ntracers)) {
    replace_regexp_file("makefile", qr/nf\s*=.*/, "nf = $ntracers");
}

if (length($ndust)) {
    replace_regexp_file("makefile", qr/ndust\s*=.*/, "ndust = $ndust");
}

if ($physics) {
    replace_regexp_file("makefile", qr/PHYSICS\s*=.*/, "PHYSICS = $physics");
}

if ($eos) {
    replace_regexp_file("makefile", qr/eos\s*=.*/, "eos = $eos");
}

if ($arch) {
    replace_regexp_file("makefile", qr/ARCH\s*=.*/, "ARCH = $arch.defs");
}

if (length($phi_dir)) {
    replace_regexp_file("makefile", qr/phi\s*=.*/, "phi = $phi_dir");
}

if (length($z_dir)) {
    replace_regexp_file("makefile", qr/z\s*=.*/, "z = $z_dir");
}

# Copy a file if it doesn't exist yet
# Usage: copy_if_not_present(filename, source directory)
# Optionally, a local filename can be specified as third argument
sub copy_if_not_present {
    my ( $filename, $location, $local_name ) = @_;

    if (!defined($local_name)) {
        $local_name = $filename;
    }

    # If the file does not exist, copy it
    unless (-e($filename)) {
        print "Getting $filename from $location/$local_name\n";
        `cp $ENV{AMRVAC_DIR}/$location/$local_name $filename`;
    }
}

# Replace lines in a file
# Usage: replace_regexp_file(filename, regexp, replacement_string)
# Note that the third argument is a **string**
sub replace_regexp_file {
    my ( $filename, $regexp, $replacement) = @_;
    my $do_replace = 1;

    @ARGV = ($filename);
    while(<>) {

        if ($do_replace) {
            s /$regexp/$replacement/; # Replace the regexp
            if (/SETVAC READS UP TO THIS POINT/) {
                $do_replace = 0; # Stop replacing
            }
        }
        print;
    }
}

# Print the current parameters
sub show_current_parameters {
    # Get a hash with the parameters
    my %params = get_current_parameters();
    my @param_names = sort keys %params;

    print " Setting         | Value\n";
    print " ----------------|------\n";

    for my $name (@param_names) {
        printf " %-15s | $params{$name}\n", $name;
    }

    print "\n Invocation of setup.pl:\n";
    printf " setup.pl -d=%d%d", $params{"ndim"}, $params{"ndir"};
    printf " -p=%s -phi=%d", $params{"physics"}, $params{"phi"};
    printf " -z=%d -eos=%s", $params{"z"}, $params{"eos"};
    printf " -nf=%d -ndust=%d", $params{"nf"}, $params{"ndust"};
    printf " -arch=%s\n", $params{"arch"};
}

# Return a hash object with the current parameters
sub get_current_parameters {
    my %params;

    if (!-e("makefile")) {
        die("Error: cannot read options; makefile is missing");
    }

    open(my $fh_makefile, "<", "makefile");
    while ($_ = <$fh_makefile>) {
        chop;

        # Note that $' returns the text after the match
        $params{"ndim"}    = $1 if /^ndim\s*=\s*(\d+)/ ;
        $params{"ndir"}    = $1 if /^ndir\s*=\s*(\d+)/;
        $params{"physics"} = $1 if /^PHYSICS\s*=\s*(\w+)/ ;
        $params{"arch"}    = $1 if /^ARCH\s*=\s*(\w+)/ ;
        $params{"phi"}     = $' if /^phi\s*=\s*/ ;
        $params{"z"}       = $' if /^z\s*=\s*/ ;
        $params{"nf"}      = $1 if /^nf\s*=\s*(\d+)/ ;
        $params{"ndust"}   = $1 if /^ndust\s*=\s*(\d+)*/ ;
        $params{"eos"}     = $1 if /^eos\s*=\s*(\w+)/ ;
        last if /SETVAC READS UP TO THIS POINT/;
    }
    close($fh_makefile);

    return %params;
}
