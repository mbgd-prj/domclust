#!/usr/bin/perl
while (<>) {
	($sp, $name) = split;
	$name =~ s/\([0-9]+\)$//;
	$spname = "$sp:$name";
	if (! $NUM{$spname}) {
		$NUM{$spname} = ++$NEWNUM;
	}
	if ($oldspname ne $spname) {
		if ($oldspname) {
			print "$oldspname $NUM{$spname} $cnt\n";
		}
		$oldspname = $spname;
		$cnt = 0;
	}
	$cnt++;
}
