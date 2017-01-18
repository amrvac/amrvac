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
    -arch=<name>                Use compilation flags from arch/<name>.defs

    -show                       Show current options
    -help                       Show this help message

Examples:

setup.pl -d=22 -arch=default
setup.pl -show\n";

# Locally define the variables that will hold the options
my $ndim;
my $ndir;
my $arch;
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
    "arch=s"  => \$arch,
    "show"    => \$show,
    "help"    => \$help)
    or die("Error in command line arguments\n");


# Show help if -help is given or if there are no other arguments
if ($help || !($ndim || $ndir || $show || $arch )) {
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
copy_if_not_present("makefile", "arch", "template.make");

if ($ndim) {
    replace_regexp_file("makefile", qr/NDIM\s*[:?]?=.*/, "NDIM := $ndim");
}

if ($ndir) {
    replace_regexp_file("makefile", qr/NDIR\s*[:?]?=.*/, "NDIR = $ndir");
}

if ($arch) {
    replace_regexp_file("makefile", qr/ARCH\s*[:?]?=.*/, "ARCH = $arch.defs");
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
        $params{"arch"}    = $1 if /^ARCH\s*=\s*(\w+)/ ;
    }
    close($fh_makefile);

    return %params;
}
