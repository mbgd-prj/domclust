#!/usr/bin/perl -s

$EVAL_CUT = 0.01;
$infile = $ARGV[0];

$BITUNIT = 1/3;    ## the default scoring system is in 1/3 bit units

if (! $skip_sort) {
	## eliminate comment lines (when using blastall -m 9)
	$infile = "egrep -v '^#' $infile | " .
	## excahnge name1 and name2, and also start and end positions
	q{awk ' BEGIN {OFS="\t"}
		$1<=$2 {print}
		$1>$2 {print $2,$1,$3,$4,$5,$6,$9,$10,$7,$8,$11,$12}' | } .
	## sort by name pair followed by E-value
	"sort -k 1,2 -k 11,11g | ";
}

open(IN, $infile) || die;

while (<IN>) {
	next if (/^#/);
	($qid, $sid, $ident, $alilen, $mismatch, $gap,
		$qstart, $qend, $sstart, $send, $eval, $score) = split;
	next if ($eval > $EVAL_CUT);
	if ($qid eq $prev_qid && $sid eq $prev_sid) {
		next;
	}
	$prev_qid = $qid; $prev_sid = $sid;

	my($dir) = 1;
	if ( $nucl ) {
		my($tmp);
		if ($qstart > $qend) {
			$tmp = $qstart; $qstart = $qend; $qend = $tmp;
			$dir *= -1;
		}
		if ($sstart > $send) {
			$tmp = $sstart; $sstart = $send; $send = $tmp;
			$dir *= -1;
		}
	}
	if ( $evalcorr) {
		# corrected E-value
		$eval *= $len
	}
	$score /= $BITUNIT;
	$dist = 100 - $ident;
	if ( $distconv && $dist > 0 ) {
	    # Kimura's correction formula for protein sequence distances
		$dist /= 100;
		$dist = - log( 1 - $dist - 0.2 * $dist * $dist);
		$dist *= 100;
	}
	$dist = sprintf("%.0f", $dist);
	$score = sprintf("%.0f", $score);
	print "$qid $sid $qstart $qend $sstart $send $dist $score";
	print " $dir" if ($nucl);
	print "\n";
}
close(IN);
