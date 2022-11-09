#!/usr/bin/perl -i

use strict;
use warnings;
use Getopt::Long;

my $help_message =
"Usage: setup.pl [options]      To set up an AMRVAC problem

Options:

    -d=N                        N is the problem dimension (1 to 3)
    -arch=<name>                Use compilation flags from arch/<name>.defs
    -phys=<name>                For new setups: use this physics module
    -help                       Show this help message

Examples:

setup.pl -d=2\n";


# Locally define the variables that will hold the options
my $ndim;
my $arch;
my $phys;
my $help;

# Parse the options. Some are handled with a subroutine, which can do immediate
# error checking.
GetOptions(
    "d=i"     =>
    sub {
        my ($opt_name, $opt_value) = @_;
        $opt_value =~ /([123])/ || die "-$opt_name flag incorrect\n";
        $ndim = $1;
        $ndim >= 1 && $ndim <= 3 ||
            die("1 <= ndim <= 3 does not hold\n");
    },
    "arch=s"  => \$arch,
    "phys=s"  => \$phys,
    "help"    => \$help)
    or die("Error in command line arguments\n");


# Show help if -help is given or if there are no other arguments
if ($help || !($ndim || $arch )) {
    print STDERR $help_message;
    exit;
}

# Check if the environment variable AMRVAC_DIR is defined
if (!$ENV{AMRVAC_DIR}) {
    print STDERR "Error: \$AMRVAC_DIR variable undefined, set it with:\n";
    print STDERR "export \$AMRVAC_DIR=your/amrvac/dir\n";
    exit;
}

# Copy makefile
copy_file("makefile", "arch", "amrvac.make");

unless (-e 'amrvac.h') {
# create amrvac.h to use the std preprocessor
open(my $fh, '>', 'amrvac.h') or die "amrvac.h is not created $!";
close $fh;
}

if ($ndim) {
    replace_regexp_file("makefile", qr/NDIM\s*[:?]?=.*/, "NDIM := $ndim");
}

if ($arch) {
    replace_regexp_file("makefile", qr/ARCH\s*[:?]?=.*/, "ARCH = $arch");
}

# If mod_usr.t or mod_usr.f are not present, copy a default tempate
unless (-e("mod_usr.t") || -e("mod_usr.f")) {
    print "No user files found, copy the default template? [y/n] ";

    chomp(my $yn = <STDIN>);
    unless ($yn eq "y") {
        exit;
    }

    copy_if_not_present("mod_usr.t", "src", "mod_usr_template.t");

    unless ($phys) {
        print "Please enter the physics name (e.g., hd, mhd, rho): ";
        $phys = <STDIN>;
        chomp($phys);
    }

    replace_regexp_file("mod_usr.t", qr/hd_/, "$phys"."_");
    replace_regexp_file("mod_usr.t", qr/_hd/, "_"."$phys");
}

sub copy_file {
    my ( $filename, $location, $local_name ) = @_;

    if (!defined($local_name)) {
        $local_name = $filename;
    }

    print "Getting $filename from $location/$local_name\n";
    my $output = `cp $ENV{AMRVAC_DIR}/$location/$local_name $filename`;
}

# Copy a file if it doesn't exist yet
# Usage: copy_file(filename, source directory)
# Optionally, a local filename can be specified as third argument
sub copy_if_not_present {
    my ( $filename, $location, $local_name ) = @_;

    if (!defined($local_name)) {
        $local_name = $filename;
    }

    # If the file does not exist, copy it
    unless (-e($filename)) {
        print "Getting $filename from $location/$local_name\n";
        my $output = `cp $ENV{AMRVAC_DIR}/$location/$local_name $filename`;
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
