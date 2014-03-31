#!/usr/bin/awk -f
{
if (substr($0, 1, 1) == "#") print $0;
else if (match($8, /^DP=/)){	
	#split($8, a, ";");
	#if (match(a[1], /^DP=/)){
	if ($6 >= qlimit){
		print $0;
	}
}
}
