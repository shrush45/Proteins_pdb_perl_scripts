
#perl code to calculate bond length of each atom present in the peptide chain
#shrushti salunke - BIM-2022-23

open(F,"1YL0.pdb") or die "file not found\n";

while(chomp($l=<F>))
{
	if($l=~/^(ATOM)\s+(\d+)\s+((N|CA|C|O)\s+(\w+)\s(\w)\s+(\d+)\s+(.*))/)
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
	for($j=$i;$j<$xlen-1;$j++)
	{		
		$t = (($x[$i]-$x[$j])**2 + ($y[$i]-$y[$j])**2 + ($z[$i]-$z[$j])**2);
		$dis = $t**(0.5);
		#print H "distance between $atom1[$i] of $aa1[$i]$anum1[$i] and $atom1[$j] of $aa1[$j]$anum1[$j] = $dis\n";
		if($dis == 0)
		{break;}
		elsif($dis<=2.00)
		{print  "The bond lenghth of $atom[$i] of $aa[$i]$anum[$i] and $atom[$j] of $aa[$j]$anum[$j] is = $dis\n";}
		
	}
}
print "done!!\n";







































































































close F;