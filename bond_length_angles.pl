
#perl code to find bond length, bond angle and dihedral angles 
#Shrushti salunke - BIM-2022-23
#The output of the code is save in file 1YL0_result.csv

use strict;


my @x = ();
my @y = ();
my @z = ();
my @atom=();
my @aa = ();
my @anum = ();
my @xval = ();
my @yval = ();
my @zval =();
my @distance =();
my @al =();
my @am = ();
my @an = ();
my @line =();
my @plane = ();
my @degrees =();
my @anormal = ();
my @bnormal =();
my @cnormal =();
my @degrees_d =();

use Math::Trig;

open(F,"C:\\Users\\Shrushti\\Downloads\\4grl.pdb") or die "file not found\n";
open(H,">C:\\Users\\Shrushti\\Downloads\\4grl_torsion.csv") or die "fnf\n";

my @a =();
while(chomp(my $l=<F>))
{
	#to extract the atom and co-ordinate lines
	if($l=~/^(ATOM)\s+(\d+)\s+((N|CA|C)\s+(\w+)\s(\w)\s+(\d+)\s+(.*))/)
	{
		push@a,$3;
	}
}
foreach my $i(0..@a)
{
	my @m = split/\s+/,$a[$i];
	push@x, $m[4];			#to extract X co-ordinates of each atom
	push@y, $m[5];			#to extract y co-ordinates of each atom
	push@z, $m[6];			#to extract z co-ordinates of each atom
	push@atom,$m[0];		#to extract atom
	push@aa,$m[1];			#to extract amino acid three letter name
	push@anum,$m[3];		#to extract amino acid number
}
 my $xlen = @x;
foreach my $i(0..$xlen-3)
{
	my $x = $x[$i] - $x[$i+1];		#x1-x2
	push@xval,$x;
	my $y = $y[$i] - $y[$i+1];		#y1-y2
	push@yval,$y;
	my $z = $z[$i] - $z[$i+1];		#z1-z2
	push@zval,$z;
	my $distance = ($x**2+$y**2+$z**2)**(0.5);
	#print "The distance between $atom[$i] of $aa[$i]$anum[$i] and $atom[$i+1] of $aa[$i+1]$anum[$i+1]  = $distance\n";
	push@distance,$distance;
	#print H " x2-x1 = $x\ty2-y1 = $y\tz2-z1 = $z\n";
	#if($distance == 0){break;}
	#else{
	my $al= $x/$distance;			#l1
	push@al, $al;
	my $am = $y/$distance;			#m1
	push@am, $am;
	my $an = $z/$distance;			#n1
	push@an, $an;		
	my $line = $atom[$i] . $anum[$i] . $atom[$i+1] . $anum[$i+1];
	push@line,$line;

	my $plane = $atom[$i] . $anum[$i] . $atom[$i+1] . $anum[$i+1] . $atom[$i+2] . $anum[$i+2];
	push@plane,$plane;
	#}
		       
}

my $len = @al;
for(my $i=0;$i<$len-1;$i++)
{ 
	#bond angle calculations
	my $cos_a = $al[$i]*(-$al[$i+1]) + ($am[$i])*(-$am[$i+1]) + $an[$i]*(-$an[$i+1]);		#to find the dot product of l1 m1 n1 and l2 m2 n2
	#print H "$cos_a\n";
		
	my $angle = acos($cos_a);									#to find the cos inverse 
	my $degrees  = rad2deg($angle);									#to convert radian to degrees
	#print  "The angle between $line[$i] and $line[$i+1] = $degrees\n";
	push@degrees,$degrees;
	
}

my $xval = @xval;
for(my $j=0;$j<=$xval-1;$j++)
{
	#for normal vector points - a,b,c - cross product calculations
	my $a = $yval[$j]*$zval[$j+1] - $yval[$j+1]*$zval[$j];		 #i
	my $b = (-($xval[$j]*$zval[$j+1] - $xval[$j+1]*$zval[$j]));	 #j
	my $c = $xval[$j]*$yval[$j+1] - $yval[$j]*$xval[$j+1]; 		 #k
	#print "$a\t$b\t$c\n";
	
	
	my $t = ($a**2+$b**2+$c**2);
	my $dis = $t**(0.5);
	#print "$dis\n";
	if($dis != 0)
	#{break;}
	{
	my $l = $a/$dis;		#l1 for normal
	push@anormal,$l;
	my $m = $b/$dis;		#m1 for normal
	push@bnormal,$m;
	my $n = $c/$dis;		#n1 for normal
	push@cnormal,$n;	
	#print "$a\t$b\t$c\n";
	#print "l = $l, m = $m, n = $n\n";
	}
}
my $length = @anormal;

for(my $k=0;$k<=$length-1;$k++)
{
	#torsion angles calculations
	my $cos_d = ($anormal[$k]*$anormal[$k+1]) + ($bnormal[$k]*$bnormal[$k+1]) + ($cnormal[$k]*$cnormal[$k+1]);		#to find the dot product of l1 m1 n1 and l2 m2 n2
	#print  "$cos_d\n";
		
	my $angle_d = acos($cos_d);							#to find the cos inverse 
	my $degrees_d  = rad2deg($angle_d);						#to convert radian to degrees
	#print "The angle between plane $plane[$k] and $plane[$k+1] = $degrees_d\n";
	push@degrees_d,$degrees_d;
	
	
}	


print H "AAR,AAR number,atom1,atom2,bond length,\t,\t,atom1atom2,atom2atom3,bond angle,\t,\t,plane1,plane2,dihedral angle\n";
foreach my $i(0..$xlen-3)
{
	print H "$aa[$i],$anum[$i],$atom[$i],$atom[$i+1],$distance[$i],\t,\t,$line[$i],$line[$i+1],$degrees[$i],\t,\t,$plane[$i],$plane[$i+1],$degrees_d[$i]\n";
}

print "done!!\nThe result is in 1YL0_result.csv file.\n";


close F;
close H;









