#!/usr/bin/perl

for ($i = 0; $i < 256; $i++) {
	$b = 1;
	$cnt = 0;
	for ($j = 0; $j < 8; $j++) {
		$cnt += (($i & $b)!=0);
		$b <<= 1;
	}
	print "," if ($i > 0);
	print "\n" if ($i % 30 == 0);
	print "$cnt";
}
