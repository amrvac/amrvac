#!/usr/bin/perl

# This script writes updated user source and par files, and warns when manual
# work is required. The updated files get a prefix $new_prefix (= new_)

use strict;
use warnings;
use Cwd;
use File::Copy;

my $backup_suffix = ".bak";

# Replace words (pattern => replacement)
my %simple_replacements = (
    qr/\bamrlist\b/ => "meshlist",
    qr/dixB/ => "nghostcells",                   # Don't use word boundary here
    qr/%w0/ => "%B0",                   # Don't use word boundary here
    qr/\berrorestimate\b/ => "refine_criterion", # (TODO: of type string)
    qr/\bfilenameini\b/ => "restart_from_file",
    qr/\bfilenameout\b/ => "base_filename",
    qr/\bflags\(\b/ => "w_for_refine(", # Don't replace the word flags (too generic)
    qr/\bmxnest\b/ => "refine_max_level",
    qr/\bngridshi\b/ => "max_blocks",
    qr/\bnormt\b/ => "time_convert_factor",
    qr/\bnormvar(0)\b/ => "length_convert_factor",
    qr/\bnormvar\b/ => "w_convert_factor",
    qr/\bnxblock/ => "block_nx",
    qr/\bnxlone/ => "domain_nx",
    qr/\bsmallrho\b/ => "small_density",
    qr/\bsmallp\b/ => "small_pressure",
    qr/\bsmallT\b/ => "small_temperature",
    qr/(?<!['.])\bt\b/ => "global_time", # Prevent 't and .t from matching
    qr/\btmax\b/ => "time_max",
    qr/\bitmax\b/ => "it_max",
    qr/\btreset\b/ => "reset_time",
    qr/\brestart_reset_time\b/ => "reset_time",
    qr/\bitreset\b/ => "reset_it",
    qr/\btmaxexact\b/ => "time_max_exact",
    qr/\bresetgrid\b/ => "reset_grid",
    qr/\btypeadvance\b/ => "time_integrator",
    qr/\btime_integrator\b/ => "time_stepper",
    qr/\btypefull1\b/ => "flux_scheme",
    qr/\btypegradlimiter1\b/ => "gradient_limiter",
    qr/\btypelimiter1\b/ => "limiter",
    qr/\btypetvdlf\b/ => "typeboundspeed",
    qr/\btypeB\b/ => "typeboundary",
    qr/\bwflags\b/ => "w_refine_weight",
    qr/\bwnames\b/ => "w_names",
    qr/\bwritew\b/ => "w_write",
    qr/\btol\b/ => "refine_threshold",
    qr/\btolratio\b/ => "derefine_ratio",
    qr/\btypegridfill\b/ => "prolongation_method",
    qr/\bb1_\b/ => "mag(1)",
    qr/\bb2_\b/ => "mag(2)",
    qr/\bb3_\b/ => "mag(3)",
    qr/(\bv1_\b|\bm1_)/ => "mom(1)",
    qr/(\bv2_\b|\bm2_)/ => "mom(2)",
    qr/(\bv3_\b|\bm3_)/ => "mom(3)",
    qr/mygeo/ => "block",
    qr/normvar\(0\)/ => "length_convert_factor",
    qr/normvar\(rho_\)/ => "w_convert_factor(rho_)",
    qr/normvar\(p_\)/ => "w_convert_factor(p_)",
    qr/normvar\(1:nw\)\b/ => "w_convert_factor(:)",
    qr/mhcgspar/ => "hydrogen_mass_cgs",
    qr/kbcgspar/ => "kboltzmann_cgs",
    qr/include 'amrvacdef.f'/ => "use mod_global_parameters",
    qr/INCLUDE:amrvacmodules\/cooling.t/ => "",
    qr/INCLUDE:amrvacmodules\/heatconduct.t/ => "",
    qr/INCLUDE:amrvacmodules\/viscosity.t/ => "",
    qr/INCLUDE:amrvacmodules\/gravity.t/ => "",
    qr/INCLUDE:amrvacmodules\/handle_particles.t/ => "",
    qr/INCLUDE:amrvacmodules\/integrate_particles.t/ => "",
    qr/INCLUDE:amrvacmodules\/magnetofriction.t/ => "INCLUDE:includes/magnetofriction.t",
    qr/INCLUDE:amrvacmodules\/fff.t/ => "INCLUDE:includes/fff.t",
    qr/INCLUDE:amrvacmodules\/pfss.t/ => "INCLUDE:includes/pfss.t",
    qr/INCLUDE:amrvacnul\/usrflags.t/ => "",
    qr/INCLUDE:amrvacnul\/specialbound.t/ => "",
    qr/INCLUDE:amrvacnul\/specialini.t/ => "",
    qr/INCLUDE:amrvacnul\/specialimpl.t/ => "",
    qr/INCLUDE:amrvacnul\/speciallog.t/ => "",
    qr/INCLUDE:amrvacnul\/specialsource.t/ => "",
    qr/\btypeaxial\b/ => "coordinate",
    );

# Replace words only used in par files (pattern => replacement)
my %par_file_replacements = (
    qr/\b13\*/ => "20*",   # nlevelshi went from 13 to 20
    qr/ *useprimitive *=.*\n/ => "", # Remove useprimitive = ... lines
    qr/ *filenamelog *=.*\n/ => "", # Remove useprimitive = ... lines
    qr/ *fileheadout *=.*\n/ => "", # Remove fileheadout = ... lines
    qr/ *snapshotini *=.*\n/ => "", # Remove snapshotini = ... lines
    qr/ *typeaxial *=.*\n/ => "", # Remove typeaxial = ... lines
    qr/ *ssplitdivb *=.*\n/ => "", # Remove ssplitdivb = ... lines
    qr/ *primnames *=.*\n/ => "", # Remove primnames = ... lines
    qr/ *w_names *=.*\n/ => "", # Remove w_names = ... lines
    qr/ *ssplituser *=.*\n/ => "", # Remove ssplituser = ... lines
    qr/ *fixsmall *=.*\n/ => "", # Remove fixsmall = ... lines
    qr/ *nghostcells *=.*\n/ => "", # Remove nghostcells = ... lines
    qr/ *typeparIO *=.*\n/ => "", # Remove typeparIO = ... lines
    qr/ *strictsmall *=.*\n/ => "", # Remove strictsmall = ... lines
    qr/ *strictgetaux *=.*\n/ => "", # Remove strictgetaux = ... lines
    qr/ *nflatgetaux *=.*\n/ => "", # Remove nflatgetaux = ... lines
    qr/ *\bconduction\b *=.*\n/ => "", # Remove conduction = ... lines
    qr/ *TCsaturate *=.*\n/ => "", # Remove TCsaturate = ... lines
    qr/ *fixprocess *=.*\n/ => "", # Remove fixprocess= ... lines
    qr/ *sliceascii *=.*\n/ => "", # Remove sliceascii= ... lines
    qr/ *restrictprimitive *=.*\n/ => "", # Remove restrictprimitive= ... lines
    qr/ *amrentropy*=.*\n/ => "", # Remove amrentropy= ... lines
    qr/ *w_for_refine\(1\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(2\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(3\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(4\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(5\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(6\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(7\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(8\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(9\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(10\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(11\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(12\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(13\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(14\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(15\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(16\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(17\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(18\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(19\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/ *w_for_refine\(20\) *=.*\n/ => "", # Remove w_for_refine*= ... lines
    qr/tsave\(1\)/ => "tsave_log",
    qr/tsave\(2\)/ => "tsave_dat",
    qr/tsave\(3\)/ => "tsave_slice",
    qr/tsave\(4\)/ => "tsave_collapsed",
    qr/tsave\(5\)/ => "tsave_custom",
    qr/\btypetvdlf\b/ => "typeboundspeed",
    qr/\bitmax\b/ => "it_max",
    qr/\bglm1\b/ => "glm",
    qr/ *typelimited *=.*\n/ => "", # Remove typelimited= ... lines
    );

# List of regular expressions and associated warnings
my %warnings = ();

# List of regular expressions and associated warnings in par files
my %par_file_warnings = (
    qr/\buseprimitive\b/ => "useprimitive is now always .true.",
    qr/\bfilenamelog\b/ => "filenamelog has been removed, log files now have the\
same name as the simulation output (controlled by base_filename)",
    qr/\bnormvar/ => "TODO: normvar(1:nw) is now w_convert_factor and
                     normvar(0) is now length_convert_factor",
    qr/\btypeboundary *=/ => "typeboundary (formerly typeB) can now be split in typeboundary_min1, typeboundary_max1, typeboundary_min2, etc.",
    );

# Get list of files, but exclude previous output of script
my @in_files = (glob "*.par *.t");

foreach my $fname (@in_files) {
    my $backup = $fname.$backup_suffix;
    my $data = read_file($fname);
    my $orig_data = $data;
    my $cwd = getcwd;

    foreach my $pat (keys %simple_replacements) {
        my $repl = $simple_replacements{$pat};
        $data =~ s/$pat/$repl/g;
    }

    # Only do these replacements in par files
    if ($fname =~ /\.par$/) {
        foreach my $pat (keys %par_file_replacements) {
        my $repl = $par_file_replacements{$pat};
        $data =~ s/$pat/$repl/g;
        }

        foreach my $pat (keys %par_file_warnings) {
            my $msg = $par_file_warnings{$pat};
            if ($data =~ /$pat/g) {
                my $cwd = getcwd;
                print STDERR "** Warning ** in file $cwd/$fname:\n";
                print STDERR "$msg\n";
            }
        }
    }

    foreach my $pat (keys %warnings) {
        my $msg = $warnings{$pat};
        if ($data =~ /$pat/g) {
            print STDERR "** Warning ** in file $cwd/$fname:\n";
            print STDERR "$msg\n";
        }
    }

    if ($data ne $orig_data) {
        print "Modifying $cwd/$fname\n";

        if (-f "$backup") {
            print "Already have a backup, not storing a new one\n";
        } else {
            copy("$fname", "$backup") or die "Backup failed: $!\n";
        }
        write_file($fname, $data);
    }
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
