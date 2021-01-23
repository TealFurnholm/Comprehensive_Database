#NOTE: THIS IS ONLY FOR USE REGARDING THE PUBLICATION - TESTING IT VS THE HUMAN ORAL MICROBIOME DATABASE - NOT FOR ANY SUBSEQUENT DATABASE CREATION OR UPDATING

$|=1;

#use warnings;

$inlin = 'TAXONOMY_DB_2018.txt';
open(INLIN, $inlin)||die;
while(<INLIN>){
		if($_ !~ /\w/){next;}
		$_ = uc($_);
		@stuff = split("\t",$_);
		$stuff[$#stuff] =~ s/[\r\n]+//;
		$tid = shift(@stuff);
		$PHY{$tid} = join(";", @stuff);
}

$outmusc = 'HOMDB_SEQS.txt';
open(OUTMUSC,">", $outmusc)||die;
$outfast = 'HOMDB_SEQS.fasta';
open(OUTFAST, ">", $outfast)||die;
$outdb = 'HOMDB_FROM_HUGEB.txt';
open(OUTDB, ">", $outdb)||die;
$outread = 'HOMDB_READS.fastq';
open(OUTREAD, ">", $outread)||die;

$done =0;
$I = "I" x 100;
$ingen = 'HOMDB_FROM_HUGE.txt';
open(INGEN, $ingen)||die;
while(<INGEN>){
	if($_ !~ /\w/){next;}
	$_ = uc($_);
	@stuff = split("\t",$_);
	$stuff[$#stuff] =~ s/[\r\n]+//;
	
	#FIX SOURCE DUPLICATION
	@GENES = split(";", $stuff[2]);
	@LOCI = ();
	%GLIST =();
	foreach(@GENES){ $GLIST{$_}=1; $_ =~ /\%([^\%]+)\%/; $locus = $1;  $locus =~ s/\&.*//; push(@LOCI, $locus); }
	$kc = keys %GLIST;
	#REMOVE DUPLICATE LOCI
	$p = $#LOCI-1;
	foreach my $x (0..$p){
		$q = $x+1;
		foreach my $y ($q..$#LOCI){
			if($LOCI[$x] eq $LOCI[$y]){ 
				if($GENES[$x] !~ /^[ABEM]\|GENE/){ delete($GLIST{$GENES[$y]}); $kc = keys %GLIST;}
				elsif($GENES[$y] !~ /^[ABEM]\|GENE/){ delete($GLIST{$GENES[$x]}); $kc = keys %GLIST;}
				else{  delete($GLIST{$GENES[$y]}); $kc = keys %GLIST;}
			}
		}
	}					
	@GENES = ();
	@TIDS = ();
	foreach my $gene (keys %GLIST){ $gene =~ s/^I/M/; $gene =~ /\|(\d+)\%/; $tid = $1; push(@GENES, $gene); push(@TIDS, $tid);}
	$stuff[2]=join(";", @GENES);
	my %HGS;
	@TIDS = grep { !$HGS{$_}++ } @TIDS;
	$stuff[3]=join(";", @TIDS);
	
	
	#ADD "OTHER" ANNOTATION
	$len = length($stuff[22]);
		$stuff[22] =~ s/(^X+|X+$)//g;
	if($stuff[22] =~ /^[MV]/){$PROTTYPE="MREG";}
	elsif($stuff[22] =~ /X/){$PROTTYPE="MIDSTOP";}
	elsif($len <=10 ){ $PROTTYPE="SHORT";}
	else{$PROTTYPE="OTHER"; if($stuff[8] eq "NORMAL"){$stuff[8]="OTHER";}}
	
	
	#PRINT SEQS
	$output1 = join("\t", @stuff);
	print OUTMUSC "$stuff[0]\t$stuff[23]\n";
	print OUTFAST ">$stuff[0]\n$stuff[23]\n";
	print OUTDB "$output1\n";
	$doing++;
			
	#SKIP BAD GENES OR IF MAX READS
	if($done == 400){print "done $done skipping doing $doing\n"; next;}
	if($stuff[3] =~ /\;/){ next;} #NOT SINGLE SOURCE GENE
	if($PHY{$stuff[3]} =~ /(YERSINIA|MYCOBACTERIUM|AGROBACTERIUM|RALSTONIA|RHODOBACTER|BACILLUS_ANTHRACIS)/){next;}
	if($stuff[6] !~ /^([^\;]+\;[^\;]+\;[^\;]+\;[^\;]+\;[^\;]+)/){next;}
	if($TIDGC{$stuff[3]} == 250){next; }  #HAVE MAXIMUM NUMBER OF GENES/READS FOR THE TID

	#make reads for the gene
	my %keep;
	$count = 0;
	$kc = keys %keep;
	while( $kc < 20 ){
		
		@reads = ($stuff[23] =~ /[ATGC]{100}/g);
		foreach(@reads){ if(!exists($track{$_})){$keep{$_}=1;} }
		$kc = keys %keep;
		
		$stuff[23]=~ s/^...//;
		$count++;
		if($count > 10){last;}
	}
	$kc = keys %keep;
	if($kc<20){next;}
	$count=0;
	foreach my $read (keys %keep){ 
		$track{$read}="@".$stuff[0]."-".$stuff[3]."\n$read\n+\n$I\n"; 
		$count++; if($count == 20){last;} }
	$tkc = keys %track;
	
	
	#increment the tid since the gene produced 10 good reads
	$TIDGC{$stuff[3]}++;
	if($TIDGC{$stuff[3]} == 250){ $DONE{$stuff[3]}=1; $done = keys %DONE; next; } #the tid has 500 genes/5000 reads so add to done list
	$done = keys %DONE;
	
	if($on%1000==0){$time = localtime; $tkc = keys %track; 
	print "on $on time $time doing $doing tkc $tkc done $done tidgc $TIDGC{$stuff[3]}\n"; }
	$on++;
}


foreach my $read (keys %track){ print OUTREAD "$track{$read}"; }

