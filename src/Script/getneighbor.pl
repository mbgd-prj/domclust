#!/usr/bin/perl
while (<>) {
	if (/Cluster ([0-9]+)/) {
		$clustid = $1;
	} elsif (/^ *\[([0-9]+)\]\([^,]+,([^,]+)\)/) {
		$hier = $1;
		$nodeid = $2;
		$nodeclust{$nodeid} = $clustid;
	} elsif (/Left:.*\(([^)]*)\)/) {
		if ($hier == 1) {
			$leftnode{$clustid} = $1;
		}
	} elsif (/Right:.*\(([^)]*)\)/) {
		if ($hier == 1) {
			$rightnode{$clustid} = $1;
		}
	} elsif (/^ *([a-z0-9]+:.*)/) {
		push(@{$clustent[$clustid]}, $1);
	}
}

$maxclust = $clustid;

foreach $clid (1..$maxclust) {
	$clid0 = $clid;
	while ($nodeclust{$leftnode{$clid0}} && $clid0 != $nodeclust{$leftnode{$clid0}}) {
		$clid0 = $nodeclust{$leftnode{$clid0}};
	}
	do {
		last if ($DONE{$clid0});
		&print_cluster($clid0);
		$DONE{$clid0} = 1;
	} while ($clid0 = $nodeclust{$rightnode{$clid0}});
}

sub print_cluster {
	my($clid) = @_;
	print "Cluster $clid";
	if ($leftnode{$clid} && $nodeclust{$leftnode{$clid}}) {
		print " L:$nodeclust{$leftnode{$clid}}";
	}
	if ($rightnode{$clid} && $nodeclust{$rightnode{$clid}}) {
		print " R:$nodeclust{$rightnode{$clid}}";
	}
	print "\n";
	foreach $ent (@{$clustent[$clid]}) {
		print "$ent\n";
	}
	print "\n";
}
