#!/usr/bin/perl

$PAT_BEGIN = "^\s*#ifdef *WITH_(NEIGHBOR|DOMCUT|CALNEWALI|DEBUG)";
$PAT_END = "^\s*#endif";
$PAT_COUNT = "^\s*#ifn*def";

while(<>){
	if(/$PAT_BEGIN/){
		$cut = 1;
		$count++;
	} elsif($cut && /$PAT_COUNT/){
		$cut++;
	} elsif($cut && /$PAT_END/){
		$cut--;
		next;
	}
	if (! $cut) {
		print;
	}
}
