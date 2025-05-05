

#to check for chainbreaks in protein structures(polypeptide chain)

open(F,"5ksa.pdb") or die "file not found\n";
open(G, ">output.txt") ;
@l = <F>;
$len = @l;
for($i=0;$i<$len;$i++)
{
	#print "$l\n";
	if($l[$i]=~/(^HEADER)\s+(.*)/)
	{
		print "protein name and PDB ID of protein is $2\n";
	}
	if($l[$i]=~/^(ATOM)\s+\d+\s+((C|N)\s+(\w+)\s((\w)\s+(\d+))\s+(.*))/)
	{
		push@a,$2;		#rest of the atom line is treated as one array element
	}
}
#print "@a\n";
$alen = @a;
for($i=0;$i<=$alen-1;$i++)
{
	@m = split/\s+/,$a[$i];
	push@x,$m[4];
	push@y,$m[5];
	push@z,$m[6];
	push@atom,$m[0];
	push@aar,$m[1];
	push@chain,$m[2];
	push@aar_num,$m[3];
	
}
#print "$atom[2]\t$aar[2]\t$chain[2]\t$aar_num[2]\n";

$xlen=@x;

for($i=0
;$i<=$xlen-2;$i++)
{
	if($aar_num[$i]!=$aar_num[$i+1])
	{
	$t = ($x[$i]-$x[$i+1])**2+($y[$i]-$y[$i+1])**2+($z[$i]-$z[$i+1])**2;
	$dis = $t**(0.5);
	if($dis>=1.35)
	{print G "There is a chainbreak between $atom[$i] of $aar[$i]$aar_num[$i] chain $chain[$i] and $atom[$i+1] of $aar[$i+1]$aar_num[$i+1] chain $chain[$i+1] separated by distance = $dis angs\n"; 
	}
	}
	
}
print "The result is store in the file output.txt\n";
print "done!!!\n";
close G;
close F;
