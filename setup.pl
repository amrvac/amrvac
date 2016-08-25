#!/usr/bin/perl -si

my $help_message =
"Usage: setup.pl [options]      To generate a new MPI-AMRVAC setup
        setup.pl -s             To show the current options

Options:

    -d=NM                       N is the the problem dimension (1 to 3)
                                M is the vector dimension (1 to 3)
    -g=<n_1, ..., n_N>          The size of a grid block (including ghostcells),
                                separated by commas, e.g.: -g=14,14 in 2D
    -p=<physics module>         Which physics module to use (rho, mhd, ...)
    -phi={1,2,3}                Index of vector phi-component (default: 3)
    -z={1,2,3}                  Index of vector z-component (default: 2)
    -eos=<equation of state>    The equation of state
    -nf=<number>                The number of fluid tracers (default: 0)
    -ndust=<number>             The number of dust species (default: 0)
    -arch=<name>                Use compilation flags from arch/<name>.defs

    -s                          Show current options

Examples:

setup.pl -d=22 -g=100,70 -phi=0 -z=2 -p=mhd -u=nul -arch=default
setup.pl -s\n";

# Check if at least one argument was defined
if (!($d || $g || $s || $p || $eos || length($phi)) &&
    !(length($z) || $arch || length($nf) || length($ndust))) {
    print STDERR "Error: no arguments have been specified\n\n";
    print STDERR $help_message;
    exit;
}

# Check if the environment variable AMRVAC_DIR is defined
if (!$ENV{AMRVAC_DIR}) {
    print STDERR "Error: AMRVAC_DIR environment variable undefined\n";
    print STDERR "You set the variable with:\n";
    print STDERR "export $AMRVAC_DIR=your/amrvac/dir";
    exit;
}

# Show parameters and quit when -s was given
if ($s) {
    show_current_parameters();
    exit;
}

# Get these files if they do not exist already
copy_if_not_present("makefile", "arch", "make_temp");
copy_if_not_present("definitions.h", "src");
copy_if_not_present('mod_indices.t', "src");
copy_if_not_present("definitions.h", "src");
copy_if_not_present("amrvacsettings.t", "src");

if ($d) {
    # Separate the -d=NM argument into ndim (N) and ndir (M)
    $c = $d - 10 * int($d/10);    # ndir is d mod 10
    $d = int($d/10);              # ndim is d/10

    replace_regexp_file("makefile", qr/ndim\s*=.*/, "ndim = $d");
    replace_regexp_file("makefile", qr/ndir\s*=.*/, "ndir = $c");
}

if ($g) {
    # Edit amrvacsettings.t specifying ixGhi[1,2,3] according to the $g flag
    # TODO: place these statements on separate lines
    @g = split(',', $g);        # Split $g into an array
    my $new_size = sprintf("ixGhi1 = %d", $g[0]);

    for ($i = 1; $i <= $#g; $i++) {
        # Concatenate other dimensions
        $new_size .= sprintf(", ixGhi%d = %d", $i+1, $g[i]);
    }
    replace_regexp_file("amrvacsettings.t", qr/ixGhi1\s*=.*/, $new_size);
}

if (length($nf)) {
    replace_regexp_file("makefile", qr/nf\s*=.*/, "nf = $nf");
}

if (length($ndust)) {
    replace_regexp_file("makefile", qr/ndust\s*=.*/, "ndust = $ndust");
}

if ($p) {
    replace_regexp_file("makefile", qr/PHYSICS\s*=.*/, "PHYSICS = $p");
}

if ($eos) {
    replace_regexp_file("makefile", qr/eos\s*=.*/, "eos = $eos");
}

if ($arch) {
    replace_regexp_file("makefile", qr/ARCH\s*=.*/, "ARCH = $arch.defs");
}

if (length($phi)) {
    replace_regexp_file("makefile", qr/phi\s*=.*/, "phi = $phi");
}

if (length($z)) {
    replace_regexp_file("makefile", qr/z\s*=.*/, "z = $z");
}

# Perform validity checks
check_validity();

# Check the validity of the current settings
sub check_validity {
    my %params = get_current_parameters();

    if ($params{"ndim"} < 1 || $params{"ndim"} > 3) {
        show_current_parameters();
        exit(print "Error: 1 <= ndim <= 3 does not hold\n");
    }

    if ($params{"ndim"} > $params{"ndir"} || $params{"ndir"} > 3) {
        show_current_parameters();
        exit(print "Error: ndim <= ndir <= 3 does not hold\n");
    }

    if ($params{"phi"} == $params{"z"} && $params{"phi"} != 0) {
        show_current_parameters();
        exit(print "Error: phi and z have to differ unless both are 0\n");
    }

    if ($params{"phi"} > 0 && $params{"phi"} !=2 && $params{"phi"} !=3) {
        show_current_parameters();
        exit(print "phi can only be 2, 3, or <= 0\n");
    }

    if ($params{"z"} > 0 && $params{"z"} !=2 && $params{"z"} !=3) {
        show_current_parameters();
        exit(print "z can only be 2, 3, or <= 0\n");
    }
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

    for my $name (@param_names) {
        printf " %-15s = $params{$name}\n", $name;
    }
}

# Return a hash object with the current parameters
sub get_current_parameters {
    my %params;

    open(makefile, "makefile");
    while ($_ = <makefile>) {
        chop;

        # Note that $' returns the text after the match
        $params{"ndim"}  = $1 if /^ndim\s*=\s*(\d+)/ ;
        $params{"ndir"}  = $1 if /^ndir\s*=\s*(\d+)/;
        $params{"phys"}  = $1 if /^PHYSICS\s*=\s*(\w+)/ ;
        $params{"arch"}  = $1 if /^ARCH\s*=\s*(\w+)/ ;
        $params{"phi"}   = $' if /^phi\s*=\s*/ ;
        $params{"z"}     = $' if /^z\s*=\s*/ ;
        $params{"nf"}    = $1 if /^nf\s*=\s*(\d+)/ ;
        $params{"ndust"} = $1 if /^ndust\s*=\s(\d+)*/ ;
        $params{"eos"}   = $1 if /^eos\s*=\s*(\w+)/ ;
        last if /SETVAC READS UP TO THIS POINT/;
    }
    close(makefile);

    # Read the grid size from amrvacsettings.t
    if (-e("amrvacsettings.t")) {
        open(vacdef, "amrvacsettings.t");
        while ($_ = <vacdef>) {
            if (/ixGhi1/) {
                my @block_size = ($_ =~ /ixGhi[123]\s*=\s*(\d+)/g);
                $params{"block_size"} = join(", ", @block_size);
            }
        }
        close(VACDEF);
    }
    return %params;
}
