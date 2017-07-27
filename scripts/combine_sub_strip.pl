#!/usr/bin/perl
#
# Mesfin Tsige
# 08/30/2006 
# SIUC
#  
# USAGE combine_sub_strip.pl <Substrate> <film> <gap>

#
sub read_masses{
   my @line;
   while(!/[0-9].*/){
     $_ = <FILE>;
   }
   while(/[0-9].*/){
     @line=split;
     $countM ++;
     $Masses[$countM] = @line[1];
     $_ = <FILE>;
   }
}
sub read_pair_coeff{
   my @line;
   while(!/[0-9].*/){
     $_ = <FILE>;
   }
   while(/[0-9].*/){
     @line=split;
     $countP ++;
     if($countP==1){$pinput=$#line; printf STDERR "pair length $pinput\n";}
     for($i=1; $i <= $pinput; $i++){
      $pairv[$countP][$i] = $line[$i];
     }
     $_ = <FILE>;
   }
}

sub read_bond_coeff{
   my @line;
   while(!/[0-9].*/){
     $_ = <FILE>;
   }    
   while(/[0-9].*/){
     @line=split;
     $countB ++;
     if($countB==1){$binput=$#line; printf STDERR "bond length $binput\n";}
     for($i=1; $i <= $binput; $i++){
      $bondv[$countB][$i] = $line[$i];
     }
   
     $_ = <FILE>;
   }
}

sub read_angle_coeff{
   my @line;
   while(!/[0-9].*/){
     $_ = <FILE>;
   }
   while(/[0-9].*/){
     @line=split;
     $countA ++;
     if($countA==1){$ainput=$#line; printf STDERR "angle length $ainput\n";}
     for($i=1; $i <= $ainput; $i++){
      $anglev[$countA][$i] = $line[$i];
     }
     $_ = <FILE>;
   }
}
sub read_dihed_coeff{
   my @line;
   while(!/[0-9].*/){
     $_ = <FILE>;
   }
   while(/[0-9].*/){
     @line=split;
     $countD ++;
     if($countD==1){$dinput=$#line; printf STDERR "dihedral length $dinput\n";}
     for($i=1; $i <= $dinput; $i++){
      $dihedv[$countD][$i] = $line[$i];
     }
     $_ = <FILE>;
   }
}
sub read_atom_info{
  my ($i, @line);
#  $Nmol1=0;
printf STDERR "Starting value of atoms $countatoms\n";
  while(!/[0-9].*/){
        $_=<FILE>;
  }
  while(/[0-9].*/){
   @line=split;	
   $countatoms++;
   if($k==2){
	$line[0]=$line[0]+$natoms1;
    $type[$line[0]]=@line[2]+$ntypes1;
    $molecule[$line[0]]=@line[1]+$Nmol1;
   }
   else{
    $type[$line[0]]=@line[2];
    $molecule[$line[0]]=@line[1];
   }
    $charge[$line[0]]=@line[3];
   for($i=4; $i < 7; $i++){
     $pos[$line[0]][$i-4]=@line[$i];
   }
   if(scalar(@line) > 8){
    $foldT=1;
    for($i=7; $i <= scalar(@line); $i++){
     $fold[$line[0]][$i-7]=@line[$i];
    }
   }
   if($k==1){  #substrate
	# This will allow us to combine both substrate and film or two bulks. But, we have remove the periodicity now.
	$pos[$line[0]][2]=$pos[$line[0]][2]+$fold[$line[0]][2]*$Zbox1;
    if($pos[$line[0]][2] < $zmin1){$zmin1=$pos[$line[0]][2];}
    if($pos[$line[0]][2] > $zmax1){$zmax1=$pos[$line[0]][2];}
    if($molecule[$line[0]] > $Nmol1){$Nmol1=$molecule[$line[0]];}
   }
   elsif($k==2){ #slab
	# This will allow us to combine both substrate and film or two bulks. But, we have remove the periodicity now.
	$pos[$line[0]][2]=$pos[$line[0]][2]+$fold[$line[0]][2]*$Zbox2;
    if($pos[$line[0]][2] < $zmin2){$zmin2=$pos[$line[0]][2];}
    if($pos[$line[0]][2] > $zmax2){$zmax2=$pos[$line[0]][2];}
   }
   $_ = <FILE>;
  }
printf STDERR "Counted so far is $countatoms\n";
printf STDERR "Counted no of chains of substrate is %6d\n",$Nmol1;
}

sub read_angles{
 my @line;
 while(!/[0-9].*/){
  $_=<FILE>;
 }
 while(/[0-9].*/){
  @line=split;
  $countangles++;
  if($k==2){
   $angleT[$countangles]=@line[1]+$atypes1;
   $angle1[$countangles]=@line[2]+$natoms1;
   $angle2[$countangles]=@line[3]+$natoms1;
   $angle3[$countangles]=@line[4]+$natoms1;
  }
  else{
   $angleT[$countangles]=@line[1];
   $angle1[$countangles]=@line[2];
   $angle2[$countangles]=@line[3];
   $angle3[$countangles]=@line[4];
  }
  $_=<FILE>;
 }
}

sub read_dihed{
  my @line;
  while(!/[0-9].*/){
        $_=<FILE>;
  }
  while(/[0-9].*/){
   @line=split;
   $countdihed++;
   if($k==2){ 
    $dihedT[$countdihed]=@line[1]+$dtypes1;
    $dihed1[$countdihed]=@line[2]+$natoms1;
    $dihed2[$countdihed]=@line[3]+$natoms1;
    $dihed3[$countdihed]=@line[4]+$natoms1;
    $dihed4[$countdihed]=@line[5]+$natoms1;
   }
   else{
    $dihedT[$countdihed]=@line[1];
    $dihed1[$countdihed]=@line[2];
    $dihed2[$countdihed]=@line[3];
    $dihed3[$countdihed]=@line[4];
    $dihed4[$countdihed]=@line[5];
   }
   $_=<FILE>;
  }
}

sub read_bonds{
  my @line;
  $maxz=1.0;
  $minz=0.0;
  while(!/[0-9].*/){
   $_=<FILE>;
  }
  while(/[0-9].*/){
    @line=split;
    $countbonds++;
    if($k==2){ 
     $bondT[$countbonds]=@line[1]+$btypes1;
     $Bond1[$countbonds]=@line[2]+$natoms1;
     $Bond2[$countbonds]=@line[3]+$natoms1;
    }
    else{ 
     $bondT[$countbonds]=@line[1];
     $Bond1[$countbonds]=@line[2];
     $Bond2[$countbonds]=@line[3]; 
    }
    $_=<FILE>
  }
}


sub print_atoms{
  my $i;
  my $xbox;
  my $ybox;
  my $zbox;
  if($countatoms != ($natoms1 + $natoms2)){
   printf STDERR "ERROR: Counted different number of atoms %8d, %8d\n",$countatoms,$natoms1 + $natoms2;
   exit;
  }
  $strlength[0]=$Box1_[1] - $Box1_[0];
  $length[0]=$Box2_[1] - $Box2_[0];
  $strlength[1]=$Box1_[3] - $Box1_[2];
  $length[1]=$Box2_[3] - $Box2_[2];
  if(abs($strlength[0] - $length[0]) > 0.0001 || abs($strlength[1] - $length[1]) > 0.0001){
   printf STDERR "ERROR: x and y dimensions of substrant and stripe are not the same!\n";
  printf STDERR "%11.5f != %11.5f\n",$strlength[0],$length[0];
  printf STDERR "%11.5f != %11.5f\n",$strlength[1],$length[1];
   exit;
  }
  $difference=$zmax1-$zmin2;
  $shift=$difference+$ARGV[2];
#  $shift=0.0;
  printf STDERR "Zmax, zmin, and shift are %11.4f %11.4f %11.4f\n",$zmax1,$zmin2,$shift;
  #open OUT2, "> lammps_noZperiod\_$ARGV[2]A\.dat" || die "Could not open out_real.\n";
  open OUT2, "> lammps_noZperiod.dat" || die "Could not open out_real.\n";
  printf OUT2 "\n\n    %8d  atoms\n",$countatoms;
  printf OUT2 "    %8d  bonds\n",$countbonds;
  printf OUT2 "    %8d  angles\n",$countangles;
  printf OUT2 "    %8d  dihedrals\n",$countdihed;
  printf OUT2 "    %8d  impropers\n\n",0;
  printf OUT2 "    %4d  atom types\n",$countP;
  printf OUT2 "    %4d  bond types\n",$countB;
  printf OUT2 "    %4d  angle types\n",$countA;
  printf OUT2 "    %4d  dihedral types\n\n",$countD;
  printf OUT2 "   %11.8f %11.8f xlo xhi\n",$Box1_[0],$Box1_[1];
  printf OUT2 "   %11.8f %11.8f ylo yhi\n",$Box1_[2],$Box1_[3];
  printf OUT2 "   %11.8f %11.8f zlo zhi\n\n",$zmin1,$zmax2+$shift;
  printf OUT2 " Masses\n\n";
  for($i=1; $i <= $countP; $i++){
   printf OUT2 " %4d %11.4f\n",$i,$Masses[$i];
  }
  printf OUT2 "\n Pair Coeffs\n\n";
  for($i=1; $i <= $countP; $i++){
   printf OUT2 " %4d",$i;
   for($j=1; $j <= $pinput; $j++){
    printf OUT2 "%11.4f",$pairv[$i][$j];
   }
   printf OUT2 "\n";
  }
  printf OUT2 "\n Bond Coeffs\n\n";
  for($i=1; $i <= $countB; $i++){
   printf OUT2 " %4d",$i;
   for($j=1; $j <= $binput; $j++){
    printf OUT2 "%11.4f",$bondv[$i][$j];
   }
   printf OUT2 "\n";
  }
  printf OUT2 "\n Angle Coeffs\n\n";
  for($i=1; $i <= $countA; $i++){
   printf OUT2 " %4d",$i;
   for($j=1; $j <= $ainput; $j++){
    printf OUT2 "%11.4f",$anglev[$i][$j];
   }
   printf OUT2 "\n";
  }
  printf OUT2 "\n Dihedral Coeffs\n\n";
  for($i=1; $i <= $countD; $i++){
   printf OUT2 " %4d",$i;
   for($j=1; $j <= $dinput; $j++){
    printf OUT2 "%11.4f",$dihedv[$i][$j];
   }
   printf OUT2 "\n";
  }
  printf OUT2 "\nAtoms\n\n";
printf STDERR "ATTENTION Nmol1 is $Nmol1\n";
  if($foldT){
   for($i=1; $i <= $countatoms; $i++){
    if($i > $natoms1){
     $pos[$i][2]=$pos[$i][2]+$shift;
    }
    printf OUT2 "%6d %3d %3d %7.4f %14.7f %14.7f %14.7f %3d %3d %3d\n",
     $i, $molecule[$i], $type[$i], $charge[$i],
     $pos[$i][0], $pos[$i][1], $pos[$i][2], $fold[$i][0],$fold[$i][1],0;
   }
  }
  else{
   for($i=1; $i <= $countatoms; $i++){
    if($i > $natoms1){
     $pos[$i][2]=$pos[$i][2]+$shift;
    }
    printf OUT2 "%6d %3d %3d %7.4f %14.7f %14.7f %14.7f\n",
     $i, $molecule[$i], $type[$i], $charge[$i],
     $pos[$i][0], $pos[$i][1], $pos[$i][2];
   }
  }
  printf OUT2 "\nBonds\n\n";
  for($j=1; $j <= $countbonds; $j++){
   printf OUT2 " %8d %3d %8d %8d\n",$j,$bondT[$j],$Bond1[$j],$Bond2[$j];
  } 
  printf OUT2 "\nAngles\n\n";
  for($j=1; $j <= $countangles; $j++){
   printf OUT2 " %8d %3d %8d %8d %8d\n",$j,$angleT[$j],$angle1[$j],$angle2[$j],$angle3[$j];
  }
  printf OUT2 "\nDihedrals\n\n";
  for($j=1; $j <= $countdihed; $j++){
   printf OUT2 " %8d %3d %8d %8d %8d %8d\n",$j,$dihedT[$j],$dihed1[$j],$dihed2[$j],$dihed3[$j],$dihed4[$j];
  }
  close OUT2;
}
#
# MAIN STARTS HERE
#
$countM=0;
$countP=0;
$countB=0;
$countA=0;
$countD=0;
$countatoms=0;
$coountbonds=0;
$countangles=0;
$countdihed=0;
if(scalar(@ARGV) < 3){
   printf STDERR "\nUSAGE:   remove_periodic_z_cpp_PE.pl <initial_data_file> <dump_file>  <time step>\n\n";
   exit;
 }
for($k = 1; $k <= scalar(@ARGV); $k++){
 open FILE, "<$ARGV[$k-1]" || die "Could not open $ARGV[$k-1]";
 while(<FILE>){
  if(/atoms/){
   split;
   if($k==1){ 
    $natoms1 = $_[0];
    printf STDERR "number of atoms in the substrate is $natoms1\n";
   }
   else{ 
    $natoms2 = $_[0];
    printf STDERR "number of atoms in the film is $natoms2\n";
   }
  }
  if(/bonds/){
   split;
   if($k==1){$nbonds1=$_[0];}
   else{$nbonds2=$_[0];}
  }
  if(/angles/){
   split;
   if($k==1){ $nangles1=$_[0];}
   else{ $nangles2=$_[0];}
  }
  if(/dihedrals/){
   split;
   if($k==1){$ndihedrals1=$_[0];}
   else{$ndihedrals2=$_[0];}
  }
  if(/atom types/){
   split;
   if($k==1){$ntypes1=$_[0];}
   else{$ntypes2=$_[0];}
  }
  if(/bond types/){
   split; 
   if($k==1){$btypes1=$_[0];}
   else{$btypes2=$_[0];}
  }
  if(/angle types/){
   split;
   if($k==1){$atypes1=$_[0];}
   else{$atypes2=$_[0];}
  }
  if(/dihedral types/){
   split; 
   if($k==1){$dtypes1=$_[0];}
   else{$dtypes2=$_[0];}
  }
  if(/xlo/){
   split;
   if($k==1){
    $Box1_[0]=$_[0];
    $Box1_[1]=$_[1];
	printf STDERR "Check this $Box1_[0] $Box1_[1]\n";
   }
   else{
    $Box2_[0]=$_[0];
    $Box2_[1]=$_[1];
   }
  }
  if(/ylo/){
   split;
   if($k==1){
    $Box1_[2]=$_[0];
    $Box1_[3]=$_[1];
	printf STDERR "Check this $Box1_[2] $Box1_[3]\n";
   }
   else{
    $Box2_[2]=$_[0];
    $Box2_[3]=$_[1];
   }
  }
  if(/zlo/){
   split;
   if($k==1){
    $Box1_[4]=$_[0];
    $Box1_[5]=$_[1];
	$Zbox1=$Box1_[5]-$Box1_[4];
	$zmin1=$zmin2=$_[1];
	$zmax1=$zmax2=$_[0];
	printf STDERR "Check this $Box1_[4] $Box1_[5]\n";
   }
   else{
    $Box2_[4]=$_[0];
    $Box2_[5]=$_[1];
	$Zbox2=$Box2_[5]-$Box2_[4];
   }
  }
  if(/Masses/){
   printf STDERR "Reading masses\n";
   read_masses();
  }
  if(/Pair Coeffs/){
   printf STDERR "Reading Pair coeffs\n",
   read_pair_coeff();
  }
  if(/Bond Coeffs/){
   printf STDERR "Reading Bond coeffs\n",
   read_bond_coeff();
  }
  if(/Angle Coeffs/){
   printf STDERR "Reading Angle coeffs\n",
   read_angle_coeff();
  }
  if(/Dihedral Coeffs/){
   printf STDERR "Reading Dihedral coeffs\n",
   read_dihed_coeff();
  }
  if(/Atoms/){
   printf STDERR "Reading charge info\n";
   read_atom_info();
  }
  if(/Bonds/){
   printf STDERR "Reading bonds\n";
   read_bonds();
  }
  if(/Angles/){
   printf STDERR "Reading angles\n";
   read_angles();
  }
  if(/Dihedrals/){
   printf STDERR "Reading dihedrals\n";
   read_dihed();
  }
 }
}
close FILE;
print_atoms();
printf STDERR "\n\nDONE SUCCESSFULLY\n\n";
