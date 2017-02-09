#!/usr/bin/perl

# This script writes updated user source and par files, and warns when manual
# work is required. The updated files get a prefix $new_prefix (= new_)

use strict;
use warnings;
use File::Copy;

my $backup_suffix = ".bak";

# Replace words (pattern => replacement)
my %simple_replacements = (
    qr/\bamrlist\b/ => "meshlist",
    qr/dixB/ => "nghostcells",                   # Don't use word boundary here
    qr/\berrorestimate\b/ => "refine_criterion", # (TODO: of type string)
    qr/\bfilenameini\b/ => "restart_from_file",
    qr/\bfilenameout\b/ => "base_filename",
    qr/\bflags\(\b/ => "w_for_refine(", # Don't replace the word flags (too generic)
    qr/\bmxnest\b/ => "refine_max_level",
    qr/\bngridshi\b/ => "max_blocks",
    qr/\bnormt\b/ => "time_convert_factor",
    qr/\bnormvar(0)\b/ => "length_convert_factor",
    qr/\bnxblock/ => "block_nx",
    qr/\bnxlone/ => "domain_nx",
    qr/\bsmallrho\b/ => "small_density",
    qr/\bsmallp\b/ => "small_pressure",
    qr/\bsmallT\b/ => "small_temperature",
    qr/(?<!['.])\bt\b/ => "global_time", # Prevent 't and .t from matching
    qr/\btmax\b/ => "time_max",
    qr/\btreset\b/ => "time_reset",
    qr/\btmaxexact\b/ => "time_max_exact",
    qr/\btypeadvance\b/ => "time_integrator",
    qr/\btypefull1\b/ => "flux_scheme",
    qr/\btypegradlimiter1\b/ => "gradient_limiter",
    qr/\btypelimiter1\b/ => "limiter",
    qr/\bwflags\b/ => "w_refine_weight",
    qr/\bwnames\b/ => "w_names",
    qr/\bwritew\b/ => "w_write",
    qr/\btol\b/ => "refine_threshold",
    qr/\btolratio\b/ => "derefine_ratio",
    qr/\btypegridfill\b/ => "prolongation_method",
    qr/\b13\*/ => "20*",   # nlevelshi went from 13 to 20
    qr/ *useprimitive *=.*\n/ => "", # Remove useprimitive = ... lines
    qr/ *filenamelog *=.*\n/ => "", # Remove useprimitive = ... lines
    # qr/tsave\(1\)/ => "tsave_log",
    # qr/tsave\(2\)/ => "tsave_dat",
    # qr/tsave\(3\)/ => "tsave_slice",
    # qr/tsave\(4\)/ => "tsave_collapsed",
    # qr/tsave\(5\)/ => "tsave_analysis"
    # qr/\s*filenamelog.*=.*/ => "", # Remove filenamelog
    );

# List of regular expressions for which the associated functions will be called. They should return a pattern.
my %replacement_functions = (
    # qr/\bnormvar/ => sub {
    #     return qq{TODO: normvar(1:nw) is now w_convert_factor and
    #                  normvar(0) is now length_convert_factor\n@_};},
    # qr/\btypeB\b/ => sub {
    #     return qq{TODO: typeB is now split in typeboundary_min1,
    #                   typeboundary_max1, typeboundary_min2, etc.\n@_};}
    # qr/\btypeB\b/ => sub {
    #     return qq{TODO: typeB is now split in typeboundary_min1,
    #                   typeboundary_max1, typeboundary_min2, etc.\n@_};}
    );

# List of regular expressions and associated warnings
my %warnings = (
    qr/\buseprimitive\b/ => "useprimitive is now always .true.",
    qr/\bfilenamelog\b/ => "filenamelog has been removed, log files now have the\
same name as the simulation output (controlled by base_filename)",
    # qr/\bnormvar/ => "TODO: normvar(1:nw) is now w_convert_factor and
    #                  normvar(0) is now length_convert_factor",
    # qr/\btypeB/ => "TODO: typeB is now split in typeboundary_min1,
    #                   typeboundary_max1, typeboundary_min2, etc."
    qr/\btypeB\b/ => "typeB can now be split in typeboundary_min1,
typeboundary_max1, typeboundary_min2, etc."
    );


# Get list of files, but exclude previous output of script
my @in_files = (glob "*.par *.t");

foreach my $fname (@in_files) {
    print "Updating par file: $fname\n";

    my $backup = $fname.$backup_suffix;

    if (-f "$backup") {
        print "Already have $backup\n";
    } else {
        copy("$fname", "$backup") or die "Backup failed: $!\n";
    }

    my $data = read_file($fname);

    foreach my $pat (keys %simple_replacements) {
        my $repl = $simple_replacements{$pat};
        $data =~ s/$pat/$repl/g;
    }

    foreach my $pat (keys %replacement_functions) {
        my $func = $replacement_functions{$pat};
        $data =~ s/$pat/&$func($&)/ge;
    }

    foreach my $pat (keys %warnings) {
        my $msg = $warnings{$pat};
        print STDERR "Warning: $msg\n" if ($data =~ /$pat/g);
    }

    write_file($fname, $data);
}

sub read_file {
    my ($filename) = @_;
    open my $in, '<:encoding(UTF-8)', $filename or
        die "Could not open '$filename' for reading $!";
    local $/ = undef;           # To read file at once
    my $all = <$in>;
    close $in;

    return $all;
}

sub write_file {
    my ($filename, $content) = @_;

    open my $out, '>:encoding(UTF-8)', $filename or
        die "Could not open '$filename' for writing $!";;
    print $out $content;
    close $out;
}
