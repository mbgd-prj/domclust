#!/usr/bin/perl -s

$tab1 = $ARGV[0];
$tab2 = $ARGV[1];
@Colors = ('#ff0000', '#00ff00', '#0000ff', '#ff00ff', '#ffff00', '#00ffff');
@TexColors = ('red', 'green', 'blue', 'magenta', 'cyan', 'yellow');

$tabsize = 50 if (! $tabsize);

#%skipspec = (uur=>1,nme=>1,cje=>1,buc=>1,xfa=>1,pae=>1,vch=>1,dra=>1,ape=>1,pab=>1,mpn=>1,eco=>1);

&read_clust($tab1, *Clust1, *DomNum1);
&read_clust($tab2, *Clust2, *DomNum2);

sub read_clust {
	my($tab, $Clust, $DomNum) = @_;
	open(T, $tab) || die;
	while (<T>) {
		chop;
		if (/^Cluster\s*(\d+)/) {
			if ($clustnum) {
				$Clust[$clustnum] = $mem;
			}
			$clustnum = $1;
			$mem = [];
		} else {
			($d) = split;
			if ($d =~ /:/) {
				($sp,$d) = split(/:/,$d);
				$sp =~ s/hpy2/hpj/g;
				$d =~ tr/a-z/A-Z/;
				$d = "$sp:$d";
			} else {
				$d =~ tr/a-z/A-Z/;
				$d = "$sp:$d";
			}
			if ($d =~ /\(([0-9]+)\)$/) {
				$domnum = $1;
				$d =~ s/\([0-9]+\)$//;
				$DomNum{$d,$ln} = $domnum;
			}
			push(@{$mem}, $d);
			push(@{$ClNum{$d}}, $ln);
		}
	}
	close(T2);
}

&TexHeader if($out eq 'tex');

for ($i = 1; $i < @Clust1; $i++) {
	foreach $d (@{$Clust1[$i]}) {
		foreach $clnum (@{$ClNum{$d}}) {
			$Cnt{$clnum}++;
		}
	}
	$j = 1;
	foreach $k (sort {$Cnt{$b}<=>$Cnt{$a}} (keys %Cnt)) {
		$Conv{$k} = $j;
		$j++;
	}
	$OutTab[$ln]->[0] = $Data[0];
	for ($i = 1; $i < @Data; $i++) {
		$sp = $spec[$i];
		$data = $Data[$i];
		next if (!$sp);
		foreach $gn (split(/ /, $data)) {
			$gn =~ s/^([a-z0-9]+)://;
			$d = $gn;
			$d =~ tr/a-z/A-Z/;
			if ($d =~ /\(([0-9]+)\)$/) {
				$d =~ s/\(([0-9]+)\)$//;
			}
			$d = "$sp:$d";
			$min = 10000; $maxclnum = 0;
			foreach $clnum (@{$ClNum{$d}}) {
				if ($min > $Conv{$clnum}) {
					$min = $Conv{$clnum};
					$maxclnum = $clnum;
				}
			}
			$OutTab[$ln]->[$i] .= " " if ($OutTab[$ln]->[$i]);
			$OutTab[$ln]->[$i] .= "$gn";
                               if ($domnum = $DomNum{$d,$maxclnum}) {
                                       $OutTab[$ln]->[$i] .= "($domnum)";
                               }
			if ($min != 10000) {
				$FoundCl{$maxclnum} = 1;
				$OutTab[$ln]->[$i] .= "\[$min]";
			}
			$FoundGn{$d} = 1;
		}
	}
	foreach $clnum (keys %FoundCl) {
		$min = $Conv{$clnum};
		foreach $gn (@{$Clust[$clnum]}) {
			($sp,$orf) = split(/:/, $gn);
			next if (! $spn{$sp} || ! $spn1{$sp});
			$NotFound[$ln]->[$spn{$sp}] .= " "
				if ($NotFound[$ln]->[$spn{$sp}]);
			if ($FoundGn{$gn}) {
			} else {
				$NotFound[$ln]->[$spn{$sp}] .= "$orf";
                               	if ($domnum = $DomNum{$gn,$clnum}) {
                                      	    $NotFound[$ln]->[$spn{$sp}]
					.= "($domnum)";
                               	}
				$NotFound[$ln]->[$spn{$sp}] .= "[$min]";
			}
		}
	}

}

if ($ln) {
	&tableOut;
}

&TexFooter if($out eq 'tex');


sub tableOut {
	local($i, $j) = @_;
	if ($out eq 'tex') {
		&TexOut($i, $j);
	} else {
		&HtmlOut($i, $j);
	}
}
sub HtmlOut {
	local($i, $j);
	print "<HTML><BODY>\n";
	print "<h3>$Title</h3>\n" if ($Title);
	print "<h4>$SubTitle</h4>\n" if ($SubTitle);
	print "<TABLE BORDER=1>\n";
	print "<TR>\n";
	for ($j = 0; $j < @spec; $j++) {
#		next if ($statout && $skipspec{$spec[$j]});
		next if (! $spn{$spec[$j]} || ! $spn1{$spec[$j]});
		print "<TH>$spec[$j]</TH>\n";
	}
	print "</TR>\n";
	for ($i = 1; $i < @OutTab; $i++) {
		print "<TR>\n";
		$diff1 = $diff2 = $diff3 = 0;
		for ($j = 0; $j < @spec; $j++) {
#			next if ($statout && $skipspec{$spec[$j]});
			next if (! $spn{$spec[$j]} || ! $spn1{$spec[$j]});
			print "<TD>";
			$cnt0 = $cnt1 = $cnt2 = 0;
			foreach $gnam (split(/ /, $OutTab[$i]->[$j])) {
				if (($gn) = ($gnam =~ /\[([0-9]+)\]/)) {
					$col = $Colors[$gn-1];
				} else {
					$col = "#000000";
				}
				if ($gn <= 6) {
					$gnam =~ s/\[[0-9]+\]//;
				}
				if ($j > 0 && $statout) {
					if ($gn == 1) {
						$cnt0++;
					} else {
						$cnt1++;
					}
				} elsif ($j > 0 && $symbout) {
					print "<font color=$col><b>*</b></font> ";
				} else {
					print "<font color=$col><b>$gnam</b></font> ";
				}
			}
			print "\n";
			foreach $gnam (split(/ /, $NotFound[$i]->[$j])) {
				if (($gn) = ($gnam =~ /\[([0-9]+)\]/)) {
					$col = $Colors[$gn-1];
				} else {
					$col = "#000000";
				}
				if ($gn <= 6) {
					$gnam =~ s/\[[0-9]+\]//;
				}
				if ($statout) {
					$cnt2++;
				} elsif ($symbout) {
					print "<font color=$col><b>.</b></font> ";
				} else {
					print "<font color=$col><i>$gnam</i></font> ";
				}
			}
			if ($statout) {
				if ($cnt0) {
				  print "<font color=$Colors[0]><b>*</b></font> ";
				}
				if ($cnt1) {
				  print "<font color=$Colors[2]><b>*</b></font> ";
				}
				if ($cnt2) {
				  print "<font color=$Colors[0]><b>.</b></font> ";
				}
				if (! $cnt0 && $cnt1 && ! $cnt2) {
					$diff1++;
				} elsif (! $cnt0 && ! $cnt1 && $cnt2) {
					$diff2++;
				} elsif ($cnt1 || $cnt2) {
					$diff3++;
				}
			}
			print "</TD>\n";
		}
		print "</TR>\n";
		if ($statout) {
			if (! $diff1 && ! $diff2) {
				if (! $diff3) {
					$ExactMatch++;
				} else {
					$Match++;
				}
			}
			$Total++;
		}
	}
	print "</TABLE>\n";
	print "</BODY></HTML>\n";
}

sub TexHeader {
	print "\\documentclass{article}\n";
	print "\\usepackage[dvips]{color}\n";
	print "\\setlength{\\textwidth}{17cm}\n";
	print "\\setlength{\\textheight}{22cm}\n";
	print "\\begin{document}\n";
}
sub TexFooter {
	print "\\end{document}\n";
}
sub TexOut {
	local($i, $j);
	local($format);
	$format = "c|";
	for ($i = 1; $i < @spec; $i++) {
		$format .= "c";
	}
	print "\\section{$Title}\n" if ($Title);
	print "\\subsection{$SubTitle}\n" if ($SubTitle);
	print "\\tiny\n";
	print "\\begin{tabular}{$format}\n";
	print "\\hline\n";
	for ($j = 1; $j < @spec; $j++) {
		print "& {\\bf $spec[$j]}\n";
	}
	print "\\\\\n";
	print "\\hline\n";
	for ($i = 1; $i < @OutTab; $i++) {
		for ($j = 0; $j < @spec; $j++) {
			print "&" if ($j > 0);
			foreach $gnam (split(/ /, $OutTab[$i]->[$j])) {
				$gnam =~ s/_/\\_/g;
				$col = '';
				if (($gn) = ($gnam =~ /\[([0-9]+)\]/)) {
					$col = $TexColors[$gn-1];
				}
				if ($col eq '') {
					$col = "black";
				}
				if ($j > 0 && $symbout) {
					print "\\textcolor{$col}{\\bf *} ";
				} else {
					print "\\textcolor{$col}{\\bf $gnam} ";
				}
			}
			print "\n";
			foreach $gnam (split(/ /, $NotFound[$i]->[$j])) {
				$col = '';
				if (($gn) = ($gnam =~ /\[([0-9]+)\]/)) {
					$gnam =~ s/_/\\_/g;
					$col = $TexColors[$gn-1];
				}
				if ($col eq '') {
					$col = "black";
				}
				if ($symbout) {
					print "\\textcolor{$col}. ";
				} else {
					print "\\textcolor{$col}{\\it $gnam} ";
				}
			}
		}
		print "\\\\\n";
	}
	print "\\end{tabular}\n\n";
}
