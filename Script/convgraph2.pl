#!/usr/bin/perl
while (<>) {
	next if (/^Cluster/);
	next if (/^$/);
	($n1, $n2) = split;
	$G{$n1}->{$n2} = 1;
}


print "[";
$cnt = 0;
foreach $n1 (keys %G) {
	next if ($Found{$n1});
	print "," if ($cnt > 0);
	&printNode($n1);
	$cnt++;
}
print "]\n";

sub printNode {
	local($n1) = @_;
	local($cnt);

	if ($Found{$n1}) {
		print qq{r("$n1")}
	} else {
		$Found{$n1} = 1;
		if ($n1 !~ /^[0-9]+$/) {
			$attr = qq{a("OBJECT", "$n1")};
		} else {
			$attr = '';
		}
		print qq{l("$n1", n("", [$attr], [};
		$cnt = 0;
		foreach $n2 (keys %{$G{$n1}}) {
			$Enum++;
			print "," if ($cnt > 0);
			print qq{l("$Enum", e("", [], };
			&printNode($n2);
			print "))\n";
			$cnt++;
		}
		print "]))\n";
	}
}
