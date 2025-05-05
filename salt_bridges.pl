

#to check salt bridges in the protein backbone structure
#shrushti salunke - BIM-2022-23


open(F,"1abe.pdb") or die"file not found\n";
@l = <F>;
$len = @l;

for($i=0;$i<=$len;$i++)
{
	if($l[$i]=~/^(ATOM)\s+\d+\s+((NZ|NH1|NE2|OD1|OE1)\s+(LYS|ARG|HIS|ASP|GLU)\s\w\s+\d+(.*))/)
	{
		
		push@charge,$2;						#array of charge atoms of basic and acidic amino acids
	}
}
#print "@charge\n";
$len = @charge;

for($i=0;$i<=$len-1;$i++)
{
	@b = split/\s+/,$charge[$i];
	push@batom,$b[0];
	push@baar,$b[1];
	push@bchain,$b[2];
	push@baar_num,$b[3];
	push@xb,$b[4];
	push@yb,$b[5];
	push@zb,$b[6];
}
for($i=0;$i<=$len-1;$i++)
{
	for($j=$i;$j<=$len-1;$j++)
	{
		$x = ($xb[$i]-$xb[$j])**2;
		$y = ($yb[$i]-$yb[$j])**2;
		$z = ($zb[$i]-$zb[$j])**2;
		$distance = ($x+$y+$z)**(0.5);
	if($distance == 0)
	{break;}
	elsif($distance<4)
	{print "The distance between $batom[$i]  $baar[$i]$baar_num[$i] and $batom[$j] of $baar[$j]$baar_num[$j] = $distance\n";
	 print "The saltbridge is may present between above atoms\n";}
}
}

 close F;
 





