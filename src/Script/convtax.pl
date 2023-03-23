#!/usr/bin/perl -s

$taxfile = $ARGV[0];

#$ranklist = "superkingdom,phylum,class,order,family,genus";
$outrank = 'genome' if (! defined $outrank);

foreach $sp (split(/,/, $splist)) {
	$SpList{$sp} = 1;
}
foreach $rank (split(/,/, $ranklist)) {
	$RankList{$rank} = 1;
}

open(T, $taxfile) || die;
while (<T>) {
	my($level, $name, $spl, $rank) = split(/\t/);
	my($output);
	while ($level <= $PrevLevels[$#PrevLevels]) {
		if (! $output) {
			$output = join(', ',@outlist) . "\n";
		}
		$output .= ")";
		pop(@PrevLevels);
	}
	print "$output\n" if ($output);

	if (! defined %RankList || $RankList{$rank} || $name eq 'root') {
		my($flag);
		undef(@outlist);
		foreach $sp (split(/,/, $spl)) {
			if (! defined %SpList || $SpList{$sp}) {
				$flag = 1;
if (! $rank || $rank eq $outrank) {
				push(@outlist, $sp);
}
			}
		}
		if ($flag) {
			$name .= " <$rank>" if ($with_rank);
##			print "(#$name\n";
			print "([$name]\n";
			push(@PrevLevels, $level);
		}
	}
}
close(T);
if (! $output) {
	$output = join(', ',@outlist) . "\n";
}
while (pop(@PrevLevels)) {
	$output .= ")";
}
print "$output\n" if ($output);
