#!/usr/local/bin/perl -s

package ClustTab;
sub read {
	my($class, $file, %opt) = @_;
	my($this) = {};
	bless $this, $class;
	if (! $file) {
		$file = "<&STDIN";
	}
	if ($opt{nested}) {
		$this->{nested} = 1;
		$this->{tabclass} = "NestedCluster";
	} else {
		$this->{tabclass} = "Cluster";
	}
	$this->init;
	$this->readDefault($file, %opt);
	$this;
}

sub readDefault {
	my($this, $file) = @_;
	my($clustnum, $cluster, $sp, $ent, $dom, $entname, $from, $to);
	open(F, $file) || die;
	while (<F>) {
		chomp;
		if (/^Cluster\s*(\d+)/) {
			$clustnum = $1;
			$cluster = $this->{tabclass}->new($clustnum);
			$this->{Clusters}->{$clustnum} = $cluster;
			$subcluster = 0;
		} elsif (/^SubCluster\s*(\d+)/) {
			$subclustid = $1;
			$subcluster = Cluster->new($subclustid);
			$cluster->addSubCluster($subcluster);
		} elsif (/^OutGroup/) {
			$subcluster = Cluster->new();
			$cluster->addOutGroup($subcluster);
		} elsif (/^\/\// || /^\s*$/) {
		} else {
			($entname, $from, $to) = split;
			if ($subcluster){
				$subcluster->add($entname,$from,$to);
			} else {
				$cluster->add($entname,$from,$to);
			}
		}
	}
	close(F);
}
sub init {
	my($this) = @_;
	$this->{curr_clnum} = 0;
	my @list = sort {$a<=>$b} (keys %{$this->{Clusters}});
	$this->{clustnum_list} = \@list;
}
sub next {
	my($this) = @_;
	my $clnum = $this->{clustnum_list}->[$this->{curr_clnum}++];
	return $this->{Clusters}->{$clnum};
}
sub list_all_species {
	my($this) = @_;
	my($speccnt) = {};
	$this->init;
	while ($cl = $this->next) {
		$cl->count_species($speccnt);
	}
	$this->{tabclass}->list_species( $speccnt );
}
sub printTab {
	my($this, %opt) = @_;
	my $species = $this->list_all_species;
	$this->init;
	if ($opt{format} eq 'HTML') {
		print "<TABLE BORDER=1>";
	}
	while ($cl = $this->next) {
		if ($opt{format} eq 'HTML') {
			print "<TR>";
			$cl->printHTML($species);
			print "</TR>";
		} else {
			$cl->printTab($species);
		}
	}
	if ($opt{format} eq 'HTML') {
		print "</TABLE>";
	}
}

package Cluster;
sub new {
	my($class, $clustid) = @_;
	my($this) = {};
	$this->{clustid} = $clustid;
	bless $this, $class;
}
sub add {
	my($this,$entname,$from,$to) = @_;
	my($sp, $ent, $dom) = &splitDom($entname);
	push(@{ $this->{members} }, 
		{sp=>$sp, ent=>$ent, dom=>$dom, from=>$from, to=>$to} );
}
sub get_init {
	my($this) = @_;
	$this->{getnum} = 0;
}
sub get_next {
	my($this) = @_;
	$this->{members}->[ $this->{getnum}++ ];
}
sub set_species {
	my($class, $species) = @_;
	${$class::species} = $species;
}
sub list_species {
	my($this, $speccnt) = @_;
	join(',', sort keys %{ $speccnt });
}
sub count_species {
	my($this, $speccnt) = @_;
	if (! $speccnt) {
		$speccnt = {};
	}
	foreach $m (@{ $this->{members} }) {
		$speccnt->{$m->{sp}}++;
	}
	$speccnt;
}
sub make_spectab {
	my($this) = @_;
	foreach $m (@{ $this->{members} }) {
		push(@{$this->{spectab}->{$m->{sp}}}, $m);
	}
}
sub phylopat {
	my($this) = @_;
	my($pat);
	foreach $sp (sort (keys %{$this->{species}})) {
		if ($this->{species}->{$sp}) {
			$pat .= '1';
		} else {
			$pat .= '0';
		}
	}
	$pat;
}
sub splitDom {
	my($entname) = @_;
	my($sp,$name) = split(/:/, $entname);
	my($dom);
	if ($name =~ /^(.*)\((\d+)\)$/) {
		$name = $1;
		$dom = $2;
	}
	($sp, $name , $dom);
}
sub printList {
	my($this) = @_;
	$this->get_init;
	while ($d = $this->get_next) {
		print join(' ', "$d->{sp}:$d->{ent}",
				$d->{from},$d->{to}),"\n";
	}
	print "\n";
}

sub printTab {
	my($this, $species) = @_;
	my($flag);
	$this->make_spectab;
	print $this->{clustid};
	foreach $sp (split(/,/, $species)) {
		print "\t";
		foreach $d (@{$this->{spectab}->{$sp}}) {
			print " " if ($flag);
			print "$d->{sp}:$d->{ent}";
			print "($d->{dom})" if ($d->{dom});
			$flag ++;
		}
	}
	print "\n";
}

sub printHTML {
	my($this, $species, %opt) = @_;
	my($flag);
	my($OPT_rowspan) = "rowspan=$opt{rowspan}" if ($opt{rowspan});
	$this->make_spectab;
	print "<TH $OPT_rowspan>$this->{clustid}</TH>";
	foreach $sp (split(/,/, $species)) {
		print "<TD $OPT_rowspan>";
		foreach $d (@{$this->{spectab}->{$sp}}) {
			print " " if ($flag);
			print "$d->{sp}:$d->{ent}";
			print "($d->{dom})" if ($d->{dom});
			$flag ++;
		}
		print "</TD>";
	}
	print "\n";
}

package NestedCluster;
@ISA = (Cluster);
sub addSubCluster {
	my($this,$subcluster) = @_;
	push(@{ $this->{subclusters} }, $subcluster);
}
sub addOutGroup {
	my($this,$subcluster) = @_;
	push(@{ $this->{outgroup} }, $subcluster);
}
sub count_species {
	my($this, $speccnt) = @_;
	if (! $speccnt->{outgroup}) {
		$speccnt->{outgroup} = {};
	}
	if (! $speccnt->{ingroup}) {
		$speccnt->{ingroup} = {};
	}
	foreach $c (@{ $this->{subclusters} }) {
		$c->count_species($speccnt->{ingroup});
	}
	foreach $c (@{ $this->{outgroup} }) {
		$c->count_species($speccnt->{outgroup});
	}
	$speccnt;
}
sub list_species {
	my($this, $speccnt) = @_;
	my @tmp1 = sort (keys %{ $speccnt->{ingroup} });
	my @tmp2 = sort (keys %{ $speccnt->{outgroup} });
	join(';', join(',', @tmp1), join(',', @tmp2));
}
sub printTab {
	my($this, $species) = @_;
	my($subclusters) = $this->{subclusters};
	my($num_ingroup) = 0+@{$subclusters};
	my($clustid) = $this->{clustid};
	my($in_sp,$out_sp) = split(/;/, $species);

	if (@{$subclusters}) {
		print "$clustid\t"; $clustid = '';
		$subclusters->[0]->printTab($in_sp, nobreak=>1, noheader=>1);
	}
	if (@{$this->{outgroup}}) {
		print "$clustid\t"; $clustid = '';
		$this->{outgroup}->[0]->printTab($out_sp);
	}
	if (@{$subclusters}) {
		foreach $s (@{$subclusters}[1..$#{$subclusters}]) {
			print "$clustid\t"; $clustid = '';
			$s->printTab($in_sp);
		}
	}
}
sub printHTML {
	my($this, $species) = @_;
	my($subclusters) = $this->{subclusters};
	my($outgroup) = $this->{outgroup};
	my($n_ingrp) = 0+@{$subclusters};
	my($n_outgrp) = 0+@{$outgroup};
	my($n_group) = ($n_ingrp > $n_outgrp ? $n_ingrp : $n_outgrp);
	my($clustid) = $this->{clustid};
	my($in_sp,$out_sp) = split(/;/, $species);

	print "<TH rowspan=$n_group>$clustid</TH>";
	if (@{$subclusters}) {
		$subclusters->[0]->printHTML($in_sp);
	}
	if (@{$this->{outgroup}}) {
		$this->{outgroup}->[0]->printHTML($out_sp, rowspan=>$n_group, noheader=>1);
	}
	if (@{$subclusters}) {
		foreach $s (@{$subclusters}[1..$#{$subclusters}]) {
			print "\n<TR>";
			$s->printHTML($in_sp);
			print "\n</TR>";
		}
	}
}

package main;
if ($0 eq __FILE__) {
	$cltab = ClustTab->read($ARGV[0], nested=>$nested);
	$cltab->printTab(format=>HTML);
}

1;
