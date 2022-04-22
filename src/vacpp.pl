#!/usr/bin/perl -s
#############################################################################
#
# AMRVAC PreProcessor (Based on vacpp.pl by Gabor Toth)
#                  June 2020: Simplified usage (Teunissen/Keppens)
#
#    Todo: here: remove IFDEF/IFNDEF constructions
#          verify and reduce the patterns in use
#    TODO: in code: remove RAY/D1/ENERGY/EVOLVINGBOUNDARY
#
# Translates the dimension independent notation to Fortran 90 by expanding
# the Loop Annotation SYntax (LASY).
#
# Usage: vacpp.pl -d=2                              # Interactive
#        vacpp.pl -d=2 file.t > file.f              # Translation
#        vacpp.pl -d=2 -maxlen=72 file.t > file.f   # Limit line length
#############################################################################
use Cwd;

my $help_message =
"Usage: vacpp.pl [options] file

Required options:

    -d=N                        N is the problem dimension (1 to 3)

Optional options:

    -maxlen=<number>            Maximum line length (default: 78)

Examples:

    # Interactive usage
    vacpp.pl -d=2 -

    # File translation
    vacpp.pl -d=2 file.t > file.f
";

# Check presence of required arguments
if (! ($d)) {
    print STDERR "Error: not all required arguments are present\n\n";
    print STDERR $help_message;
    exit
}

$AMRVAC_DIR=$ENV{AMRVAC_DIR};

# Store current working directory
my $cwd = getcwd();

# Include these directories for header files
# TODO: check whether $p has a valid value
@INC = {"$CWD","$AMRVAC_DIR/$p"};

# Default maximum length for lines unless set by eg. "vacpp.pl -maxlen=72 ..."
$maxlen=78 unless $maxlen;

# Check whether the problem and vector dimension lie between 1 and 3, by
# matching with a regexp.
$d =~ /([123])/ || die "Incorrect -d flag value\n";

# Store the problem dimension ($1 is input argument of the matching groups
# of the regexp)
$ndim = $1;

# Call these routines
&definepatterns;

# leftover switch D1, only used in amrvacio/mod_slice.t
# to allow for nesting in lasy syntax, acts as nooned
%defvars = ();
if ($ndim==1){$defvars{D1}=1;}


# Process the files given as arguments
foreach $file (@ARGV) {
   &processfile($file, 'fh00');
}

#============================================================================

sub processfile {
   local($filename, $input) = @_;
   my $partial_line = "";

   $input++;
   open($input, $filename) || die"Can't open $filename: $!\n";

   while (<$input>) {
      # Print empty and comment lines as they are
      if(/^$/ || /^ *![^S]/ && not /^ *!\$[^S]/) {
          print;
          next;
      }

      # Hide quoted text from translation
      $_=&quotation($_);

      # Stop if there are tabs
      if (/\t/) {
          die "vacpp.pl error: $filename contains TABs";
      }

      # Collect continuation lines
      ($_, $partial_line) = contline($_, $partial_line);
      next if ($partial_line);

      # Expand the line
      &processline($input);

      # Print the line formatted according to maxlen
      &printline if $_;
   }
}

# Join continuation lines into a single line
sub contline {
    my ($line, $partial_line) = @_;

    if ($line =~ /^ *!\$OMP/ || $line =~ /^#/ || $line =~ /^ *print *\*.*& *$/ ||
        $line =~ /\/\/ *& *$/) {
        # Don't break OpenMP, std preprocessor directives (lines starting with '#'), print and //&
        $line = $partial_line.$line;
        $partial_line = "";
    } elsif ($line =~ /[^A-Z]& *$/) {
        $partial_line .= $line;
    } else {
        $line = $partial_line.$line;
        $line =~ s/([^A-Z])& *\n */$1/g;
        $partial_line = "";
    }
    return ($line, $partial_line);
}
#============================================================================
sub processline{
   local($input)=@_;
   # ATTACH A NUL-SEPARATOR CHARACTER TO THE BEGINNING OF LINE FOR LEFTMATCH
   $_=$spc.$_;
   # DO SUBSTITUTIONS FOR BLOCKS "{....}"
   while(/$ope/){
      &block($input);
   }
   # DO SUBSTITUTION FOR EXPRESSIONS "...^pattern..."
   while(/$pat/){
       &expression;
      #print "Expr:$_";
   }
   # REMOVE NUL-SEPARATOR CHARACTERS
   s/$spc//g;
   # REPLACE BREAK CHARACTERS BY NEWLINES AND INDENTATION
   if(/\\/){/^ */; $indent=$&; s/\\/\n$indent/g};
}
#============================================================================
sub definepatterns{

   $ope='{'; $clo='}'; $bar='\|';
   $pre='#'; $ifdef='IFDEF'; $ifndef='IFNDEF';
   $pat='\^'; $patid='[A-Z]'; $patchr='[A-Z&%]'; $uni='%';
   $brk='\\'; $sep='[,\+\-\*/:; ]'; $spc='~';
   $rbound=" ,;~)\n"; $lbound=" ,;~(\n";

   $nsubdefault=$ndim;
   # Define number of substitutes and substitute strings for patterns
   # E.g. ^D -> 1,2 is defined by &patdef('D',2,'1','2','3')
   &patdef('ND'	,1	,$ndim			);
   &patdef('DE'	,$ndim-1,2	,3		);
   &patdef('DE&',$ndim-1			);
   &patdef('DE%',$ndim-1,'^%2'	,'^%3'		);
   &patdef('DM&',$ndim-1			);
   &patdef('DM'	,$ndim-1,1	,2		);
   &patdef('DMB',$ndim-1,$ndim-1,$ndim-2        );
   &patdef('D'	,$ndim	,1	,2	,3	);
   &patdef('D&'	,$ndim				);
   &patdef('DB'	,$ndim	,$ndim	,$ndim-1,$ndim-2);
   &patdef('DD'	,$ndim	,'^D' 	,'^D'	,'^D'	);
   &patdef('DDB',$ndim  ,'^DB'  ,'^DB'  ,'^DB'  );
   &patdef('DD&',$ndim	,'^D&'	,'^D&'	,'^D&'	);
   &patdef('D%'	,$ndim	,'^%1'	,'^%2'	,'^%3'	);
   &patdef('S'	,$ndim	,'^LIM1:','^LIM2:','^LIM3:');
   &patdef('T'	,$ndim	,'^LLIM1:','^LLIM2:','^LLIM3:');

   &patdef('LIM'	,2	,'min'	,'max'	);
   &patdef('LLIM'	,2	,'lo'	,'hi'	);
   &patdef('L'		,2	,'min^D','max^D');
   &patdef('LL'		,2	,'lo^D'	,'hi^D');
   &patdef('LSUB'	,2	,'+'	,'-'	);
   &patdef('LADD'	,2	,'-'	,'+'	);
   &patdef('LT'		,2	,'>'	,'<'	);

   &patdef('IFONED'	,$ndim==1		);
   &patdef('IFTWOD'	,$ndim==2		);
   &patdef('IFTHREED'	,$ndim==3		);
   &patdef('NOONED'	,$ndim!=1		);
}
#============================================================================
sub patdef{
   # Put pattern definitions into global %sub and %nsub associative arrays
   local($pattern,$nsub,@substitute)=@_; local($isub);
   for($isub=0;$isub<=$#substitute;$isub++){
       $sub{$pattern,$isub+1}="$substitute[$isub]";
   }
   $nsub{$pattern}=$nsub;
}
#============================================================================
sub block{
   # Find block opening by $ope and closed by matching $clo, and substitute
   local($input)=@_;
   local($rightpos,$leftpos,$line,$nline);
   my $partial_line = "";

   $leftpos=index($_,$ope);
   $level=0;
   $rightpos=&rightmatch($leftpos,$ope,$clo,$clo);

   # If no match, ie $rightpos>=length($_), read at most 2000 more lines
   while ($rightpos>=length($_) && $nline<2000 && ($line=<$input>)) {
       # Print empty and comment lines as they are
       if(/^$/ || /^ *![^S]/ && not /^ *!\$[^S]/) {
           print;
           next;
       }

       # Hide quoted text from translation
       $line = &quotation($line);

       # Stop if there are tabs
       if ($line =~ /\t/) {
          die "vacpp.pl error: $filename contains TABs";
       }

       # Collect continuation lines
       ($line, $partial_line) = contline($line, $partial_line);
       next if ($partial_line);

       # $_.=&unquote($line);
       $_.=$line;
       $nline++;
       $rightpos=&rightmatch($rightpos-1,$ope,$clo,$clo);
   }

   die "$ope without matching $clo within 2000 lines" if $rightpos>=length($_);

   &substitute($leftpos,$rightpos);
}
#============================================================================
sub expression{
   # Find expression surrounding a bounded pattern and do substitution
   local($patpos,$rightpos,$leftpos);
   $patpos=index($_,'^');
   $leftpos = &leftmatch($patpos,'(',')',$lbound);
   $rightpos=&rightmatch($patpos,'(',')',$rbound);
   &substitute($leftpos,$rightpos);
}
#============================================================================
sub rightmatch{
   local($p,$lparen,$rparen,$bound)=@_; local($i,$c);
   for($i=$p+1;$i<length($_);$i++){
      $c=substr($_,$i,1);
      last unless ($level!=0 || index($bound,$c)<0);
      $level++ if $c eq $lparen;
      $level-- if $c eq $rparen;
   }
   $i;
}
#============================================================================
sub leftmatch{
   local($p,$lparen,$rparen,$bound)=@_; local($i,$c,$level);
   for($i=$p-1;$i>0;$i--){
      $c=substr($_,$i,1);
      last unless ($level!=0 || index($bound,$c)<0);
      $level++ if $c eq $rparen;
      $level-- if $c eq $lparen;
   }
   $i;
}
#============================================================================
sub substitute{
   #Substitute patterns of the same ID as of the first pattern
   local($leftpos,$rightpos)=@_;
   local($head,$lchr,$expr,$rchr,$tail);
   local($id,$nsub,$expa,$result);

   $head=substr($_,0,$leftpos);
   $lchr=substr($_,$leftpos,1);
   $expr=substr($_,$leftpos+1,$rightpos-$leftpos-1);
   $rchr=substr($_,$rightpos,1);
   $tail=substr($_,$rightpos+1);
   $symbol=substr($_,$leftpos+1,1);

   if($symbol==$pre){
     # Handle #IFDEF ...
     if($expr=~/^$pre\s*$ifdef\s*(\S*)\s.*/){
       $var=$1;
       if (exists $defvars{$var}){
	 $expr=~ s/$pre\s*$ifdef\s*$var(.*)/$1/;
	 $lchr='' if $lchr eq $ope || ($nsub==0 && $lchr eq $separator);
	 $rchr='' if $rchr eq $clo;
	 $_=$head.$lchr.$expr.$rchr.$tail;
	 return;
       } else {
	 $lchr='' if $lchr eq $ope || ($nsub==0 && $lchr eq $separator);
	 $rchr='' if $rchr eq $clo;
	 $_=$head.$lchr.$rchr.$tail;
	 return;
       }
     }

     # Handle #IFNDEF ...
     if ($expr=~/^$pre\s*$ifndef\s*(\S*)\s.*/) {
       $var=$1;
       if (not exists $defvars{$var}){
	 $expr=~ s/$pre\s*$ifndef\s*$var(.*)/$1/;
	 $lchr='' if $lchr eq $ope || ($nsub==0 && $lchr eq $separator);
	 $rchr='' if $rchr eq $clo;
	 $_=$head.$lchr.$expr.$rchr.$tail;
	 return;
       } else{
	 $lchr='' if $lchr eq $ope || ($nsub==0 && $lchr eq $separator);
	 $rchr='' if $rchr eq $clo;
	 $_=$head.$lchr.$rchr.$tail;
	 return;
       }
     }
   }

   #Determine separator
   if($rchr eq ";")                 {$separator=";"}
   elsif($expr=~s/$bar([^$bar]*)$//){$separator="$spc$1$spc"}
   elsif($expr=~s/\\$//)            {$separator="$spc$brk$spc"}
   elsif($expr=~s/$sep$//)          {$separator="$&"}
   else                             {$separator=","};

   # Determine number of substitutions and ID
   if ($expr=~/$pat($patid)($patchr*)/) {
       # Example: ^D& gives $1 = D and $2 = &
       $id=$1;

       if (exists $nsub{"$1$2"}) {
           $nsub=$nsub{"$1$2"};     # Look up in $nsub array
       } else {
           die "Unknown pattern $_\n";
       }
   } else {
       # Patterns such as {end do\}
       $nsub=$nsubdefault
   }

   # Do the substitutions
   for($isub=1;$isub<=$nsub;$isub++){
     $expa=$expr;
     # If there is a unitvector pattern ^%N choose preceeding or subsequent
     # part depending on N==$isub or not.
     $expa=~s/$pat$uni$isub.*// || $expa=~s/.*$pat$uni.//;
     # Substitute patterns with $pat$id to their $isub-th substitute
     $expa=~s/$pat($id$patchr*)/$sub{$1,$isub}/g;
     $result.=$separator if $isub>1;
     $result.=$expa;
   }

   $lchr='' if $lchr eq $ope || ($nsub==0 && $lchr eq $separator);
   $rchr=';' if $rchr eq $clo && $separator eq ";";
   $rchr='' if $rchr eq $clo;
   $_=$head.$lchr.$result.$rchr.$tail;
}
#============================================================================
sub printline{
    local($line,$comment); $sss="   ";

    # PRINT FORMATTED OUTPUT LINE BY LINE
    while(s/(.*)\n//){
       $line=$1;
       # PUT TRAILING COMMENTS INTO $comment
       if($line=~s/ *!.*//){
           $comment=&unquote($&);
           # If line is longer than $maxlen, try reducing length of comment
           if(length("$line$comment")>$maxlen){$comment=~s/ *! */ !/}
       }else{$comment=""};

       # Print full line
       $line = &unquote(&format90($line));
       print $line, "$comment\n";
    }
}
#===========================================================================
# Break long lines into continuation lines and/or reduce indentation
sub format90 {
   local($line)=@_;

   local($bestlen,$goodlen,$len,$maxindent,$indent,$indentnow,$c,$answer);

   # If line is not too long return
   return($line) if length($line)<=$maxlen;

   if ($line =~ /^ *!\$OMP/ || $line =~ /^#/ || $line =~ /^ *print *\*.*& *$/ ||
        $line =~ /\/\/ *& *$/ || $line =~ /error stop/) {
       # Don't break lines with OpenMP,std preprocessor directives (lines starting with '#'), print, //&, error stop
       return ($line);
   }

   # Don't break lines with print statements
   return($line) if ($line =~ /^ *(print)/);

   # Determine line indentation. If too much, reduce it to maximum.
   $maxindent=int(0.6*$maxlen);
   $line=~s/^( {0,$maxindent}) */$1/; $indent=$1;

   # We are happy if the length of the line is between $goodlen and $maxlen
   $goodlen=$maxlen / 2;

   # Start breaking line
   while(length($line)>$maxlen){
       # Remove & followed by newlines
       $line =~ s/& *\n *//g;

       # Check for semicolon
       if(($len=rindex($line,';',$maxlen-1))>=0){
          $answer.=substr($line,0,$len)."\n";
          $line=substr($line,$len+1); $line=~s/^ */$indent/;
          next;
       }
       # Find best breakpoint after indentation but before $maxlen-2
       $line=~/^ */;
       $bestlen=$indentnow=length($&);

       foreach $c (',', ' ', '=>', '/=', '>=', '<=', '==',
                   '+', '-', ' / ', '**', ' * ', '.or.',
                   '.and.', ' > ', ' < ', '(', ')') {
           # Get the start position of the last occurrence of $c beginning at or
           # before $maxlen
           $len = rindex($line, $c, $maxlen);

           # If the line is long enough with the current break, exit
           if ($len > $bestlen) {
               $bestlen = $len;
               if ($bestlen > $goodlen) {
                   # Include operator on current line
                   $bestlen += length($c);
                   last;
               }
           }
       }

       if ($bestlen>$indentnow) {
           # Collect broken parts in $answer and break the continuation further
           if (substr($line,$bestlen) =~ /^ *\n*$/) {
               # If there is only whitespace on the new line, ignore it
               $answer.=substr($line,0,$bestlen)."\n";
               $line="";
           } else {
               $answer.=substr($line,0,$bestlen)."&\n";
               $line="$indent$sss".substr($line,$bestlen);
           }
       } else {
          # No break was found. Remove indentation if there is any or die.
          $line=~s/^ +// || die "Couldn't break line:".&unquote($line)
       }
   }
   $answer.=$line;
}
#===========================================================================
sub quotation{
   # Hide quoted strings by converting them from 0-127 to 128-255 ASCII.
   # Only matched quotation marks ( ' and " ) count.
   local($line)=@_; local($head,$q);

   # From left to right invert quoted text into $head. The rest is in $line.
   while($line=~/(['"])/){
      $head.=$`; $q=$1; $line=$';
      # Check for matching quotation mark. If found, convert quoted part.
      if($line=~/$q/){$q="$q$`$q"; $line=$'; $q=~tr/\x00-\x7f/\x80-\xff/;}
      $head.=$q;
   }
   # return result
   $head.$line;

}
#===========================================================================
sub unquote{
    # Uncover hidden quotations by deleting 8-th bit.
    local($line)=@_;
    $line=~tr/\x80-\xff/\x00-\x7f/;
    $line;
}
