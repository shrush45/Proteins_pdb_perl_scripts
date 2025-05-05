
#perl code to find the bond angle between main chain atoms of a peptide chain
#shrushti salunke - BIM-2022-23
use Math::Trig;				#perl module for trigonometric functions

open(F,"1YL0.pdb") or die "file not found\n";
open(H,">1YL0_output.txt") or die "file not found\n";

while(chomp($l=<F>))
{
	if($l=~/^(ATOM)\s+(\d+)\s+((N|CA|C)\s+(\w+)\s(\w)\s+(\d+)\s+(.*))/)
	{
		push@a,$3;
	}
}
foreach $i(0..@a)
{
	@m = split/\s+/,$a[$i];
	push@x, $m[4];			#to extract X co-ordinates of each atom
	push@y, $m[5];			#to extract y co-ordinates of each atom
	push@z, $m[6];			#to extract z co-ordinates of each atom
	push@atom,$m[0];		#to extract atom
	push@aa,$m[1];			#to extract amino acid three letter name
	push@anum,$m[3];		#to extract amino acid number
}
#print "@x1\n";

$alen = @a;
$xlen = @x;
for($i=0;$i<$xlen-1;$i++)					
{		
		$a = ($x[$i]-$x[$i+1]);		#to find x1 - x2
		$b = ($y[$i]-$y[$i+1]);		#to find y1 - y2 
		$c = ($z[$i]-$z[$i+1]);		#to find z1 - z2
		$s = $a**2 + $b**2 + $c**2;
		$dis = $s**(0.5);		#distance between two atoms
		if($dis<=1.55)
		{print "The distance between $atom[$i] of $aa[$i]$anum[$i] and $atom[$i+1] of $aa[$i+1]$anum[$i+1]= $dis\n";
		if($dis == 0){break;}
		else{
		$l = $a/$dis;			#l1
		push@l, $l;
		$m = $b/$dis;			#m1
		push@m, $m;
		$n = $c/$dis;			#n1
		push@n, $n;		
		$line = $atom[$i].$anum[$i].$atom[$i+1].$anum[$i+1];
		push@line,$line;
		}
		}
}
$len = @l;
for($i=0;$i<$len-1;$i++)
{ 
	$cos = $l[$i]*(-$l[$i+1]) + ($m[$i])*(-$m[$i+1]) + $n[$i]*(-$n[$i+1]);		#to find the dot product of l1 m1 n1 and l2 m2 n2
	print H "$cos\n";
		
	$angle = acos($cos);							#to find the cos inverse 
	$degrees  = rad2deg($angle);						#to convert radian to degrees
	print H "The angle between $line[$i] and $line[$i+1] = $degrees\n";
	
}		
close F;













		
=head
for($i=0;$i<$len-1;$i++)
{
	for($j=$i+1;$j<$len-1;$j++)
	{
		$cos = $l[$i]*$l[$j] + $m[$i]*$m[$j] + $n[$i]*$n[$j];
		#print "$cos\n";
		
		$angle = acos($cos);
		$degrees  = rad2deg($angle);
		print "The angle between $line[$i] and $line[$j] = $degrees\n";
	}
}
for($j=2;$j<=$xlen-1;$j+=4)
{	
	$a1 = $x[$j] - $x[$j+2];
	$b1 = $y[$j] - $y[$j+2];			
	$c1 = $z[$j] - $z[$j+2];
	$s = $a**2 + $b**2 + $c**2;
	$dis = $s**(0.5);		#distance between two atoms
	if($dis<=1.55)
	{#print "The distance between $atom[$i] of $aa[$i]$anum[$i] and $atom[$i+1] of $aa[$i+1]$anum[$i+1]= $dis\n";
	$l1 = $a1/$dis;			#l1
	push@l1, $l1;
	$m1= $b1/$dis;			#m1
	push@m1, $m1;
	$n1 = $c1/$dis;			#n1
	push@n1, $n1;
	$line = $atom[$j].$anum[$j].$atom[$j+1].$anum[$i+1];
	push@line,$line;
	}
}
=cut
		
		
		
		
		