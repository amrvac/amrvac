#!/usr/bin/perl -si
if(!($d||$g||$s||$p||$eos||length($phi)||$u||length($z)||$arch||length($nf)||length($ndust))){
   print STDERR <<HELPMSG ;

$0 [-d=NdimNdir] [-g=n1,n2,n3] [-p=PHYSICS] [-eos=EOS] [-phi=Idirphi] [-z=Idirz] \\
       [-u=USERFILE] [-nf=nf] [-ndust=nd] [-arch=arch] [-s]

Ndim and Ndir are number of spatial DIMENSIONS and vector components,
respectively, the GRID size is defined by Ndim integers n1,n2...
PHYSICS is the extension name for amrvacphys.t, and amrvacpar.t.
Idirphi and Idirz define the indexnames for cylindrical coordinates (2,3 or 0).
USERFILE is the extension for amrvacusr.t and the optional amrvacusrpar.t file.
To see the SETTINGS use the -s flag. Dont forget the "=" signs and commas!

Examples:

setup.pl -d=22 -g=100,70 -phi=0 -z=2 -p=mhd -u=nul -arch=default
setup.pl -s
HELPMSG
exit}

if(!$ENV{AMRVAC_DIR}){
    print STDERR <<ENVMSG ;

Set the environment variable $AMRVAC_DIR to the code main directory, e.g. \

export $AMRVAC_DIR=youramrvacdir

ENVMSG
exit}

$showonly= $s && !($d||$g||$p||$u||$eos||$arch||length($nf)||length($ndust)||length($phi)||length($z));

# Get the makefile:
    $filename="makefile";
    unless (-e($filename)){                                                                        
	print "Getting makefile from source directory\n";
	`cp $ENV{AMRVAC_DIR}/arch/make_temp $filename`;
    }

# Get the definitions.h:
    $filename="definitions.h";
    unless (-e($filename)){                                                                        
	print "Getting definitions.h from source directory\n";
	`cp $ENV{AMRVAC_DIR}/src/$filename $filename`;
    }

if($d){
    # Calculate ndim=$d and ndir=$c from the $d=10*ndim+ndir flag.
    $c=$d-10*int($d/10); $d=int($d/10);
    1<=$d  && $d<=3 || exit(print "1 <= ndim <= 3 should hold!\n");
    $d<=$c && $c<=3 || exit(print "ndim <= ndir <= 3 should hold!\n");
}

if($g){
    # Edit the line in amrvacsettings.t specifying ixGhi according to the $g flag.
    @g=split(',',$g);
    $filename="amrvacsettings.t";
    @ARGV=($filename);
    unless (-e($filename)){
	print "Getting ".$filename." from source directory\n";
	`cp $ENV{AMRVAC_DIR}/src/$filename .`;
    }
    while(<>){
        if(/ixGhi1=/){
            $_=$`;for($i=0;$i<=$#g;$i++){$_ .= "ixGhi".($i+1)."=".$g[$i].",";};
            s/,$/\n/;
        }
        print;
    }
}

if($d){
    # Set switches in makefile according to $c, $d
    $filename="makefile";                                                                                                                                             
    @ARGV=($filename);                                                                                                                                                       
    $doit=1;
    while(<>){
      if($doit){
         s/ndim=.*/ndim=$d/ if($d);
         s/ndir=.*/ndir=$c/ if($c);
      }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
      print;
    }
}

if(length($nf)){
    # Set switches in makefile according $nf
    $filename="makefile";                                                                                                                                             
    @ARGV=($filename);                                                                                                                                                       
    $doit=1;
    while(<>){
      if($doit){
         s/nf=.*/nf=$nf/ ;
      }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
      print;
    }
}

if(length($ndust)){
    # Set switches in makefile according $nd
    $filename="makefile";                                                                                                                                             
    @ARGV=($filename);                                                                                                                                                       
    $doit=1;
    while(<>){
      if($doit){
         s/ndust=.*/ndust=$ndust/ ;
      }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
      print;
    }
}

if($p){
    $filename="makefile";                                                                                                                                             
    @ARGV=($filename);   
    $doit=1;
    while(<>){
        if(/PHYSICS\s*=.*/ && $doit){
            $_="PHYSICS       = $p\n"
        }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
        print;    	
    }
}

if($u){
    # Set switches in makefile according $u
    $filename="makefile";                                                                                                                                             
    @ARGV=($filename);                                                                                                                                                       
    $doit=1;
    while(<>){
      if($doit){
         s/USER=.*/USER=$u/ ;
      }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
      print;
    }
}

if($eos){
    # Set switches in makefile according $eos
    $filename="makefile";                                                                                                                                             
    @ARGV=($filename);                                                                                                                                                       
    $doit=1;
    while(<>){
      if($doit){
         s/eos=.*/eos=$eos/ ;
      }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
      print;
    }
}

if($arch){
# Put the arch into makefile:
    $filename="makefile";
    unless (-e($filename)){                                                                        
	print "Getting $filename from source directory\n";
	`cp $ENV{AMRVAC_DIR}/arch/make_temp $filename`;
    }
    @ARGV=($filename);                                                                                                                 
    $doit=1;
    while(<>){
        if(/ARCH\s*=.*/ && $doit){
            $_="ARCH          = $arch.defs\n"
        }
      $doit=0 if /SETVAC READS UP TO THIS POINT/;
        print;    	
    }
}

if(length($phi) || length($z)){
    # Set phi_ and z_ in vacpp.pl
    # Read current values
# Read the parameters set in the makefile
 open(makefile,"makefile");
while($_=<makefile>){
   chop;
   $phi0="$'" if /phi\s*=\s*/ ;
   $z0="$'" if /z\s*=\s*/ ;
   last if /SETVAC READS UP TO THIS POINT/;
}
close(makefile);
    $phi0=0 if $phi0<0; $z0=0 if $z0<0;

    # Use current values if not defined by the switches
    $phi=$phi0 unless length($phi);
    $z=$z0 unless length($z);
    if($phi != $phi0 || $z != $z0){
        # Check for validity
        exit(print "-phi=$phi and -z=$z have to differ unless both are 0\n") 
           if $phi==$z && $phi!=0;
        exit(print "-phi can only be 0, 2, or 3\n") 
           if $phi!=0 && $phi!=2 && $phi!=3;
        exit(print "-z can only be 0, 2, or 3\n") 
           if $z!=0 && $z!=2 && $z!=3;
        # Change settings in amrvacdef.t. Big negative integers are used for 
        # the switched off directions so that mphi_ and bphi_ remain negative.
        $phi=-9 if $phi==0; $z=-8 if $z==0;
        @ARGV=("makefile");
	$doit=1;
        while(<>){
	    if ($doit){
		s/phi=.*/phi=$phi/;
		s/z=.*/z=$z/;
	    }
	    $doit=0 if /SETVAC READS UP TO THIS POINT/;
	    print;
        }
    }
}
 unless ($showonly){
   $u = "nul" unless $u;
    foreach $filename ('amrvacusr.t','amrvacusrpar.t'){
      unless (-e($filename)){
	print "Getting $filename.$u from user directory\n";
	`cp $ENV{AMRVAC_DIR}/src/usr/$filename.$u $filename`;
      }else{
	print STDERR "$filename is already there, doing nothing!\n";
      } 
      unless (-e($filename)){
	print "$filename.$u not found, getting $filename.nul instead\n";	  
	`cp $ENV{AMRVAC_DIR}/src/usr/$filename.nul $filename`;
      } 
}
   foreach $filename ('mod_indices.t'){
     unless (-e($filename)){                                                
       print "Getting $filename from source directory\n";
       `cp $ENV{AMRVAC_DIR}/src/$filename $filename`;
     } else{
       print STDERR "$filename is already there, doing nothing!\n";
     }
   }
 }
# Read the parameters set in the makefile
 open(makefile,"makefile");
# Find ndim and ndir and switches
while($_=<makefile>){
   chop;
   $ndim="$'" if /ndim\s*=\s*/ ;
   $ndir="$'" if /ndir\s*=\s*/;
   $phys="$'" if /PHYSICS\s*=\s*/ ;
   $arch="$'" if /ARCH\s*=\s*/ ;
   $arch =~ s/.defs$//;
   $phi="$'" if /phi\s*=\s*/ ;
   $z="$'" if /z\s*=\s*/ ;
   $nf="$'" if /nf\s*=\s*/ ;
   $ndust="$'" if /ndust\s*=\s*/ ;
   $eos="$'" if /eos\s*=\s*/ ;
   $u="$'" if /USER\s*=\s*/ ;
   last if /SETVAC READS UP TO THIS POINT/;
}
close(makefile);
$d=10*$ndim+$ndir;
$phi=0 if $phi<0; $z=0 if $z<0;

$switch=~s/;? *\$if_/,/g; $switch=~s/^,//; $switch=~s/;$//;

# Read the grid size from amrvacsettings.t
   $filename="amrvacsettings.t";
    unless (-e($filename)){                                                                        
	print "Getting ".$filename." from source directory\n";
	`cp $ENV{AMRVAC_DIR}/src/$filename .`;
    }  
 open(VACDEF,$filename);
while(($_=<VACDEF>)!~/ixGhi/){}; 
close(VACDEF);
/(ixGhi.*)/;$_=$1;s/ixGhi[123]=//g;$g=$_;
@g=split(',',$g);
exit(print"ndim=$ndim is greater than the number of defined sizes g=$g!\n")
        if $ndim > 1+$#g;

$settings="-d=$d -phi=$phi -z=$z -g=$g -p=$phys -eos=$eos -nf=$nf -ndust=$ndust -u=$u -arch=$arch";

if($s){
    # Show the current settings
    print "$0 $settings\n";
    exit if $showonly;
}
