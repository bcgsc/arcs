#!/home/martink/perl/5.10.0/bin/perl

## EDITED SCRIPT TO INTRODUCE Z GAPS INSTEAD OF N GAPS

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
#LINKS: Long Interval Nucleotide K-mer Scaffolder 

#SYNOPSIS

#DOCUMENTATION
#   LINKS-readme.txt distributed with this software @ www.bcgsc.ca
#   http://www.bcgsc.ca/platform/bioinfo/software/links
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use LINKS, the LINKS code or ideas, please cite our work

#LICENSE
#LINKS Copyright (c) 2014-2017 British Columbia Cancer Agency Branch.  All rights reserved.
#LINKS is released under the GNU General Public License v3

use strict;
use POSIX;
use FindBin;
# use lib "$FindBin::Bin/./lib/bloomfilter/swig";
# use BloomFilter;
use Time::HiRes;
use Getopt::Std;
use Net::SMTP;
use vars qw($opt_f $opt_s $opt_d $opt_k $opt_e $opt_l $opt_a $opt_v $opt_b $opt_t $opt_p $opt_o $opt_r $opt_x $opt_m $opt_z);
getopts('f:s:d:k:e:l:a:v:b:t:p:o:r:x:m:z:');
my ($bf_file,$base_name,$distances,$k,$insert_stdev,$min_links,$max_link_ratio,$verbose,$step,$offset,$fpr,$bfoff,$readlength,$min_size)=("","",4000,15,0.1,5,0.3,0,2,0,0.001,0,0,500);
my $last_step  = $step; ### default for last_step set as step
my $version = "[v1.8.5]";
my $dev = "rwarren\@bcgsc.ca";

#-------------------------------------------------

if(! $opt_f || ! $opt_s){
   print "Usage: $0 $version\n";
   print "-f  sequences to scaffold (Multi-FASTA format, required)\n"; 
   print "-s  file-of-filenames, full path to long sequence reads or MPET pairs [see below] (Multi-FASTA/fastq format, required)\n";
   print "-m  MPET reads (default -m 1 = yes, default = no, optional)\n";
   print "\t! DO NOT SET IF NOT USING MPET. WHEN SET, LINKS WILL EXPECT A SPECIAL FORMAT UNDER -s\n";
   print "\t! Paired MPET reads in their original outward orientation <- -> must be separated by \":\"\n"; 
   print "\t  >template_name\n\t  ACGACACTATGCATAAGCAGACGAGCAGCGACGCAGCACG:ATATATAGCGCACGACGCAGCACAGCAGCAGACGAC\n";	
   print "-d  distance between k-mer pairs (ie. target distances to re-scaffold on. default -d $distances, optional)\n";
   print "\tMultiple distances are separated by comma. eg. -d 500,1000,2000,3000\n";
   print "-k  k-mer value (default -k $k, optional)\n";
   print "-t  step of sliding window when extracting k-mer pairs from long reads (default -t $step, optional)\n";
   print "\tMultiple steps are separated by comma. eg. -t 10,5\n";
   print "-o  offset position for extracting k-mer pairs (default -o $offset, optional)\n";
   print "-e  error (%) allowed on -d distance   e.g. -e 0.1  == distance +/- 10% (default -e $insert_stdev, optional)\n";
   print "-l  minimum number of links (k-mer pairs) to compute scaffold (default -l $min_links, optional)\n";
   print "-a  maximum link ratio between two best contig pairs (default -a $max_link_ratio, optional)\n";
   print "\t *higher values lead to least accurate scaffolding*\n";
   print "-z  minimum contig length to consider for scaffolding (default -z $min_size, optional)\n";
   print "-b  base name for your output files (optional)\n";
   print "-r  Bloom filter input file for sequences supplied in -s (optional, if none provided will output to .bloom)\n";
   print "\t NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME FILE SUPPLIED IN -f WITH SAME -k VALUE\n";
   print "\t IF YOU DO NOT SUPPLY A BLOOM FILTER, ONE WILL BE CREATED (.bloom)\n";
   print "-p  Bloom filter false positive rate (default -p $fpr, optional; increase to prevent memory allocation errors)\n";
   print "-x  Turn off Bloom filter functionality (-x 1 = yes, default = no, optional)\n";
   print "-v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n"; 
   die "\nError: Missing mandatory options -f and -s.\n\n";
}

my $assemblyfile = $opt_f;
my $longfile = $opt_s;
$distances = $opt_d if($opt_d);
$k = $opt_k if($opt_k);
$verbose = $opt_v if($opt_v);
$min_links = $opt_l if($opt_l);
$min_size = $opt_z if($opt_z);
$max_link_ratio = $opt_a if($opt_a);
$step = $opt_t if($opt_t);
###ADDED FOR MPET
my $readlength = $opt_m if($opt_m);###MPET
$insert_stdev = 0.5 if($opt_m);    ###MPET (Need to adjust to a wider-range of distances when dealing with MPET)
$insert_stdev = $opt_e if($opt_e); ###When set, this will override the MPET-induced changes on -e
$base_name = $opt_b if($opt_b);
$offset = $opt_o if($opt_o);
$bf_file = $opt_r if($opt_r);
$fpr = $opt_p if($opt_p);
$bfoff = $opt_x if($opt_x);

my $assemblyruninfo="";


if(! -e $assemblyfile){
   die "Invalid file: $assemblyfile -- fatal\n";
}


### Naming output files
if ($base_name eq ""){

   $base_name = $assemblyfile . ".scaff_s-" . $longfile . "_d" . $distances . "_k" . $k . "_e" . $insert_stdev . "_l" . $min_links . "_a" . $max_link_ratio . "_z" . $min_size . "_t" . $step . "_o" . $offset . "_r-" . $bf_file . "_p" . $fpr . "_x" . $bfoff . "_m" . $readlength;

   my $pid_num = getpgrp(0);
   $base_name .= "_pid" . $pid_num;
}

my $log = $base_name . ".log";
my $scaffold = $base_name . ".scaffolds";
my $issues = $base_name . ".pairing_issues";
my $distribution = $base_name . ".pairing_distribution.csv";
my $bfout = $base_name . ".bloom";
my $graph = $base_name . ".gv";
my $numnamecorr = $base_name . ".assembly_correspondence.tsv";
my $tigpair_checkpoint = $base_name . ".tigpair_checkpoint.tsv";### add a checkpoint file, prevent re-running LINKS from scratch if crash
my $simplepair_checkpoint = $base_name . ".simplepair_checkpoint.tsv";### add a checkpoint file, prevent re-running LINKS from scratch if crash

open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";

#-------------------------------------------------
#
my $init_message = "\nRunning: $0 $version\n-f $assemblyfile\n-s $longfile\n-m $readlength\n";
$init_message .= "-d $distances\n-k $k\n-e $insert_stdev\n-l $min_links\n-a $max_link_ratio\n-t $step\n-o $offset\n-z $min_size\n-b $base_name\n-r $bf_file\n-p $fpr\n-x $bfoff\n\n----------------- Verifying files -----------------\n\n";

print $init_message;
print LOG $init_message;
$assemblyruninfo=$init_message;

#-------------------------------------------------
my $file_message = "";
my $ct_fof_line = 0;

if(! -e $longfile){
   $file_message = "Invalid file of filenames: $longfile -- fatal\n";
   print $file_message;
   print LOG $file_message;
   exit;
}else{
   open(FOF,$longfile) || die "Can't open file of filenames $longfile for reading -- fatal.\n";
   while(<FOF>){
      chomp;
      $ct_fof_line++;
      $file_message = "Checking $_...";
      print $file_message;
      print LOG $file_message;
      if(! -e $_){
         $file_message = "na\n*** File does not exist -- fatal (check the path/file and try again)\n";
         if(/^\>/){
           $file_message .= "It appears that the file you supplied in -s is in fasta format where it should be a file of filenames, listing the fullpath/fasta,fastq files to consider (1 per line).\n";
         }elsif(/^\@/){
           $file_message .= "It appears that the file you supplied in -s is in fastq format where it should be a file of filenames, listing the fullpath/fasta, fastq files to consider (1 per line).\n";
         }
         print $file_message;
         print LOG $file_message;
         exit;
      }else{
         $file_message = "ok\n";
         print $file_message;
         print LOG $file_message;
      }
   }
   close FOF;
}

if(! -e $assemblyfile){
   $file_message = "\nInvalid file: $assemblyfile -- fatal\n";
   print $file_message;
   print LOG $file_message;
   exit;
}else{
   $file_message = "Checking sequence target file $assemblyfile...ok\n";
   print $file_message;
   print LOG $file_message;
}

#-------------------------------------------------
#LINKS STARTS
#-------------------------------------------------

my $date = `date`;
chomp($date);
my $bloom;

eval{

my ($contigpairs,$simplepair,$tig_length)=("","","");

my ($tighash, $tignames, $tig_length);
if(! -e $tigpair_checkpoint){###MAR2016 no assembly checkpoint file detected this is the most crucial file

   my $reading_tigbloom_message = "\n\n=>Reading contig/sequence assembly file : $date\n";
   print $reading_tigbloom_message;
   print LOG $reading_tigbloom_message;
   $assemblyruninfo.=$reading_tigbloom_message;
   if($bfoff==0){

    my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($assemblyfile);
    my $bfelements = int($size);
    my $m = ceil((-1 * $bfelements * log($fpr)) / (log(2) * log(2))); #the number of bits
    ### ensures $m is a multiple of 8 to avoid error in BloomFilter
    my $rem = 64 - ($m % 64);
    $m = $m + $rem;
    #$m = (int($m/8) + 1) * 8;### ensures $m is a multiple of 8 to avoid error in BloomFilter
    my $hashfct = floor(($m / $bfelements) * log(2));# the number of hash functions

    if($bf_file eq ""){

      my $bfmessage = "Building a Bloom filter using $k-mers derived from sequences in -f $assemblyfile...\n";
      print LOG $bfmessage;
      print $bfmessage;

      $date = `date`;
      chomp($date);

      my $bfspecs = "*****\nBloom filter specs\nelements=$bfelements\nFPR=$fpr\nsize (bits)=$m\nhash functions=$hashfct\n*****\n";
      print $bfspecs;
      print LOG $bfspecs;
      $bloom = new BloomFilter::BloomFilter($m, $hashfct, $k);
      $bloom = &contigsToBloom($assemblyfile,$k,$bloom,$hashfct,$min_size);
      $date = `date`;
      chomp($date);
      my $writing_tigbloom_message = "\n\n=>Writing Bloom filter to disk ($bfout) : $date\n";
      print $writing_tigbloom_message;
      print LOG $writing_tigbloom_message;
      $assemblyruninfo.=$writing_tigbloom_message;
      $bloom->storeFilter($bfout);
    }else{
      my $bffile_message="";
      my $bfreusemessage = "A Bloom filter was supplied ($bf_file) and will be used instead of building a new one from -f $assemblyfile\n";
      print LOG $bfreusemessage;
      print $bfreusemessage;

      if(! -e $bf_file){
         $bffile_message = "\nInvalid file: $bf_file -- fatal\n";
         print $bffile_message;
         print LOG $bffile_message;
         exit;
      }else{
         $bffile_message = "Checking Bloom filter file $bf_file...ok\n";
         print $bffile_message;
         print LOG $bffile_message;
      }

      $bffile_message="Loading bloom filter of size $m from $bf_file\n";
      print $bffile_message;
      print LOG $bffile_message;

      $bloom = new BloomFilter::BloomFilter($bf_file);
    }
   }else{###no BF functionality
    my $nobfmessage = "The Bloom filter-building functionality is turned off (-x $bfoff), skipping...\n";
    print LOG $nobfmessage;
    print $nobfmessage;
   }

   #-------------------------------------------------
   $date = `date`;
   chomp($date);

   my $reading_reads_message = "\n\n=>Reading long reads, building hash table : $date\n";
   print $reading_reads_message;
   print LOG $reading_reads_message;
   $assemblyruninfo.=$reading_reads_message;

   my $matepair;
   my $pairct;
   my @distance_array = split(/,/,$distances);
   my @step_array = split(/,/,$step);
   my $totalpairs = 0;
   my $initpos = $offset;
   my $array_element = 0;

   foreach my $frag_dist (@distance_array){### v1.8 April 2016, Iterate through distances supplied. eg. -d 500,1000,2000 

      my $ctrd = 0;
      my $delta = $frag_dist - ( 2 * $k);
      my $file_ct = 0;
      my $distpairs = 0;

      if(defined $step_array[$array_element]){
         $last_step = $step_array[$array_element];### records the last usable t params, even if array doesn't match with d
      }

      my  $processed_reads_message = "Reads processed k=$k, dist=$frag_dist, offset=$initpos nt, sliding step=$last_step nt:\n";
      print $processed_reads_message;
      print LOG $processed_reads_message;

      open(FOF,$longfile) || die "Can't open file of filenames $longfile for reading -- fatal\n";
      while(<FOF>){
         chomp;
         $file_ct++;
         #if(/\.bam/i){
         #   if (! -e $SAMPATH){#only care if there are bam files to read
         #      my $samtool_err_mess = "The executable $SAMPATH does not exist, please revise -- fatal.\n";
         #      print LOG $samtool_err_mess;
         #      print $samtool_err_mess;
         #      close LOG;
         #      exit;
         #   }
         #   ($set,$bin,$sumall,$ctall) = &readBAM($sumall,$ctall,$set,$bin,$_,$encoded,$seedsplit,$r_clip,$q_clip,$c_clip,$e_ascii,$file_ct,$ct_fof_line,$targetwordlen);
         #}else{
         #($set,$bin,$sumall,$ctall) = &readFastaFastq($sumall,$ctall,$set,$bin,$_,$encoded,$seedsplit,$r_clip,$q_clip,$c_clip,$e_ascii,$file_ct,$ct_fof_line,$targetwordlen);
         if(! -e $_){ die "WARNING: Your file $_ does not exist -- fatal.\n";    }
         if($bfoff){
            ($matepair,$pairct,$ctrd) = &readFastaFastqBFoff($file_ct,$ct_fof_line,$ctrd,$_,$frag_dist,$k,$last_step,$matepair,$delta,$initpos,$readlength);
         }else{
            ($matepair,$pairct,$ctrd) = &readFastaFastq($file_ct,$ct_fof_line,$ctrd,$_,$frag_dist,$k,$last_step,$matepair,$delta,$initpos,$bloom,$readlength);
         }
	 $distpairs+=$pairct;
         #}
         $initpos += $offset;###ensures a larger kmer spectrum is explored, conflicting kmer pairs are minimized, when $offset <> 0
      }### FOF ends
      close FOF;

      my $kmerpairmessage = "\nExtracted $distpairs $k-mer pairs at -d $frag_dist, from all $ctrd sequences provided in $longfile\n\n";
      print $kmerpairmessage;
      print LOG $kmerpairmessage;
      $totalpairs+=$distpairs;
      $array_element++;
   }### iterate through distances

   my $totalkmerpairmessage = "\nExtracted $totalpairs $k-mer pairs overall. This is the set that will be used for scaffolding\n";
   print $totalkmerpairmessage;
   print LOG $totalkmerpairmessage;

   #-------------------------------------------------
   $date = `date`;
   chomp($date);

   my $reading_seqs_message = "\n\n=>Reading sequence contigs (to scaffold), tracking k-mer positions : $date\n";
   print $reading_seqs_message;
   print LOG $reading_seqs_message;
   $assemblyruninfo.=$reading_seqs_message;
   my ($track_all);
   ($track_all,$tig_length) = &readContigs($assemblyfile,$matepair,$k,$min_size);

   #-------------------------------------------------
   $date = `date`;
   chomp($date);

   my $sc_start_message = "\n\n=>Scaffolding initiated : $date\n";
   print $sc_start_message;
   print LOG $sc_start_message;
   $assemblyruninfo.= $sc_start_message . "\n";

   ($contigpairs,$simplepair) = &pairContigs($longfile,$matepair,$track_all,$tig_length,$issues,$distribution,$totalpairs,$tigpair_checkpoint,$simplepair_checkpoint,$verbose);
   #-------------------------------------------------

   ### Read contig names and sequences from the FASTA file.
   ($tighash, $tignames, $tig_length) = &readContigsMemory($assemblyfile);
}else{######MAR2016 check point file exists, skip above
   $date = `date`;
   chomp($date);
   my $sc_start_message = "\n\nScaffolding mid-point files:\n$tigpair_checkpoint\n$simplepair_checkpoint\ndetected; LINKS will skip contig pairing and use info from files instead : $date\n";
   print $sc_start_message;
   print LOG $sc_start_message;
   $assemblyruninfo.= $sc_start_message . "\n";

   ### Read contig names and sequences from the FASTA file.
   my $message = "\nReading sequences from $assemblyfile : $date\n";
   print $message;
   print LOG $message;
   ($tighash, $tignames, $tig_length) = &readContigsMemory($assemblyfile);

   $date = `date`;
   chomp($date);
   $message = "\ndone.\nReading sequence pairs from $tigpair_checkpoint : $date\n";
   print $message;
   print LOG $message;

   ($contigpairs,$simplepair) = &readContigPairs($tigpair_checkpoint,$simplepair_checkpoint,$tig_length);
}

$date = `date`;
chomp($date);
my $message = "\ndone.\nBuilding scaffolds : $date\n";
print $message;
print LOG $message;

###Building chains
$simplepair = &buildScaffolds($contigpairs,$tig_length,$simplepair,$scaffold,$longfile,$verbose);

$date = `date`;
chomp($date);
my $gv_start_message = "\ndone.\nBuilding .gv graph : $date\n";
print $gv_start_message;
print LOG $gv_start_message;

&buildGraph($simplepair,$tig_length,$graph,$version,$base_name,$date); 

$assemblyruninfo.= $gv_start_message . "\n";

$date = `date`;
chomp($date);

my $sc_end_message = "\n\n=>Scaffolding ended : $date\nScaffolds layout in: $scaffold\nScaffold graph in: $graph\n";
print $sc_end_message;
print LOG $sc_end_message;
$assemblyruninfo.= $sc_end_message . "\n";

#-------------------------------------------------
$date = `date`;
chomp($date);

my $sc_fasta_message = "\n\n=>Making FASTA file : $date\n";
print $sc_fasta_message;
print LOG $sc_fasta_message;
$assemblyruninfo.= $sc_fasta_message . "\n";

my $scaffold_fasta = &buildScaffoldFasta($scaffold,$tighash,$assemblyfile,$numnamecorr,$tignames);

my $sc_fasta_done_message = "Scaffolds FASTA in : $scaffold_fasta\n";
print $sc_fasta_done_message;
print LOG $sc_fasta_done_message;
$assemblyruninfo.= $sc_fasta_done_message . "\n";

my $corr_message = "\n\n=>Wrote correspondence file tracking LINKScontigID <=> OriginalContigNames : $numnamecorr\n";
print $corr_message;
print LOG $corr_message;
$assemblyruninfo.= $corr_message . "\n";

};###end eval block

#-------------------------------------------------
$date = `date`;
chomp($date);

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n";
   print $failure;
   print LOG $failure;
   $assemblyruninfo.=$failure . "\n";
}else{
   my $success = "\nScaffolding executed normally $date\n";
   print $success;
   print LOG $success;
   $assemblyruninfo.=$success . "\n";
}
close LOG;

print "$0 $version terminated successfully on $date\n";
exit;


###for dev. test purposes
eval{
   my $wdir = `pwd`;
   chomp($wdir);
   my $smtp = Net::SMTP->new('mailhost');
   $smtp->mail("LINKS\@bcgsc.ca");
   $smtp->to($dev);
   $smtp->data();
   $smtp->datasend("Subject: Your LINKS run\n");
   $smtp->datasend("At: $wdir\n");
   $smtp->datasend($assemblyruninfo);
   $smtp->dataend();
   $smtp->quit;
};

exit;


#----------------
sub contigsToBloom{
   my ($file,$k,$bloom,$hashfct,$min_size) = @_;

   my $prevhead = "";
   my $seq = "";
   my $cttig=0;

   ###Support for compressed files MAR2016
   if($file=~/zip$/i){
      open(IN,"unzip -p $file|") || die "Error reading $file -- fatal\n";
   }elsif($file=~/gz$/i || $file=~/gzip$/i){
      open(IN,"gunzip -c $file|") || die "Error reading $file -- fatal\n";
   }else{
      open(IN,$file) || die "Error reading $file -- fatal\n";
   }

   my $contigs_processed_message = "Contigs (>= $min_size bp) processed k=$k:\n";
   print $contigs_processed_message;
   print LOG $contigs_processed_message;
   ###
   while(<IN>){
      chomp;
      if(/^\>(\S+)/){
         my $head=$1;

         if ($head ne $prevhead && length($seq) >= $min_size && $prevhead ne ''){
            $cttig++;
            print "\r$cttig";
            $|++;
            $bloom = &kmerizeContigBloom_newloop(uc($seq),$bloom,$hashfct,$k);
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= $_;
      }
   }
   $cttig++;
   print "\r$cttig";
   print LOG "$cttig\n\n";
   $|++;
   $bloom = &kmerizeContigBloom_newloop(uc($seq),$bloom,$hashfct,$k) if(length($seq) >= $min_size);
   close IN;

   return $bloom;
}


#----------------
sub readContigs{
   my ($file,$matepair,$k,$min_size) = @_;

   my ($track_all,$tig_length);
   my $prevhead = "";
   my $seq = "";
   my $cttig=0;

   ###Support for compressed files MAR2016
   if($file=~/zip$/i){
      open(IN,"unzip -p $file|") || die "Error reading $file -- fatal\n";
   }elsif($file=~/gz$/i || $file=~/gzip$/i){
      open(IN,"gunzip -c $file|") || die "Error reading $file -- fatal\n";
   }else{
      open(IN,$file) || die "Error reading $file -- fatal\n";
   }

   my $contigs_processed_message = "Contigs (>= $min_size bp) processed k=$k:\n";
   print $contigs_processed_message;
   print LOG $contigs_processed_message;
   ###
   while(<IN>){
      chomp;
      if(/^\>(.*)/){
         my $head=$1;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $cttig++;
            print "\r$cttig";
            $|++;

            my $tiglen = length($seq);
            if($tiglen >= $min_size){
               $track_all = &kmerizeContig(uc($seq),$track_all,$matepair,$k,$cttig,0);
               my $revcomp = &reverseComplement($seq);
               $track_all = &kmerizeContig($revcomp,$track_all,$matepair,$k,$cttig,1);
               $tig_length->{$cttig} = $tiglen;
            }
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= $_;
      }
   }
   $cttig++;
   print "\r$cttig";
   print LOG "$cttig\n\n";
   $|++;

   if(length($seq) >= $min_size){
      $track_all = &kmerizeContig(uc($seq),$track_all,$matepair,$k,$cttig,0);
      my $revcomp = &reverseComplement($seq);
      $track_all = &kmerizeContig($revcomp,$track_all,$matepair,$k,$cttig,1);
      $tig_length->{$cttig} = length($seq);
   }
   ###
   close IN;

   return $track_all,$tig_length;
}

#----------------
sub readFastaFastq{

   my ($file_ct,$ct_fof_line,$ctrd,$file,$frag_dist,$k,$step,$matepair,$delta,$initpos,$bloom,$readlength) = @_;
   ###
   my $readinput_message = "\nReads processed from file $file_ct/$ct_fof_line, $file:\n";
   print $readinput_message;
   print LOG $readinput_message;


   ###Support for compressed files MAR2016
   if($file=~/zip$/i){
      open(IN,"unzip -p $file|") || die "Error reading $file -- fatal\n";
   }elsif($file=~/gz$/i || $file=~/gzip$/i){
      open(IN,"gunzip -c $file|") || die "Error reading $file -- fatal\n";
   }else{
      open(IN,$file) || die "Error reading $file -- fatal\n";
   }

   my $prevhead = "";
   my $seq = "";
   my $pairct = 0;
   my $quad = 0;
   my $endposition = $readlength-$k;###MPET
   my $buffer=$readlength+1;###MPET

   LINE:
   while(<IN>){
      chomp;
      $quad++;
      #if(/^([ATGC]+\:[ATGC]+)$/i && $quad<4){$quad=2;}###MPET SPECIFIC IGNORE Ns in SEQUENCE (don't have to)
      if(/^([ACGTMRWSYKVHDBN\:]+)$/i && $quad<4){$quad=2;}###MPET OR FASTA/FASTQ
      if($quad==1 || $quad==3){
         #print "1||3 $_\n";
         next LINE if($_ eq "+");
         my $head = $1 if(/^\S{1}(\S+)/);
#print "$quad **$_** $prevhead ne  && $head ne $prevhead\n";
         if($prevhead ne "" && $head ne $prevhead){# && length($seq)>=50){
            #my $qclipflag = 0;
            #if($qual ne ""){$qclipflag = $r_clip;}
            #print "loadSeq: $prev with $seq and $qual\n";
#print ">$prevhead QUAD:$quad\n$seq\n";

            if(! $readlength){###NO MPET
               $endposition = (length($seq)-$frag_dist);
               $buffer=$k+$delta;
            }else{
               $readlength=length($1) if($seq=~/^([ACGTMRWSYKVHDBN]+)/i);###MPET
               $endposition = $readlength-$k;###MPET
               $buffer = $readlength+1;###MPET
            }

            ($matepair,$pairct) = &kmerize(uc($seq),$frag_dist,$k,$matepair,$step,$prevhead,$endposition,$initpos,$pairct,$bloom,$buffer,$readlength);
            $seq = "";
            $quad=1;
         }
         $ctrd++;
         print "\r$ctrd";
         $|++;
         $prevhead = $head;
      }elsif($quad==2){
         #print "SEQ $_\n";
         $seq .= $_;
      }elsif($quad==4){
         $quad=0;
      }
   }
   if($prevhead ne ""){# && length($seq)>=50){
      if(! $readlength){###NO MPET
         $endposition = (length($seq)-$frag_dist);
         $buffer=$k+$delta;
      }else{
         $readlength=length($1) if($seq=~/^([ACGTMRWSYKVHDBN]+)/i);###MPET
         $endposition = $readlength-$k;###MPET
         $buffer = $readlength+1;###MPET
      }###MPET


      ($matepair,$pairct) = &kmerize(uc($seq),$frag_dist,$k,$matepair,$step,$prevhead,$endposition,$initpos,$pairct,$bloom,$buffer,$readlength);
   }

   close IN;
   print LOG "$ctrd\n\n";

   #####

   #print "\nExtracted $pairct total $k-mer pairs.\n";  
   #print LOG "\nExtracted $pairct total $k-mer pairs.\n"; 

   return $matepair,$pairct,$ctrd;
}

#----------------
sub kmerize{

   my ($seq,$frag_dist,$k,$matepair,$step,$head,$endposition,$initpos,$pairct,$bloom,$buffer,$readlength) = @_;

   for(my $pos=$initpos;$pos<=$endposition;$pos+=$step){###MPET

      my $rd1 = substr($seq,$pos,$k);
      $rd1 = &reverseComplement($rd1) if($readlength);###MPET
      my $secondstart = $pos + $buffer;###MPET
      my $rd2ss = substr($seq,$secondstart,$k);
      my $rd2 = &reverseComplement($rd2ss);

      if(defined $matepair->{$rd2}{$rd1}){my $trd1=$rd2;my $trd2=$rd1;$rd1=$trd1;$rd2=$trd2;}

      if($bloom->contains($rd1) && $bloom->contains($rd2)){ ###Don't bother hashing k-mer pairs if assembly doesn't have these kmers
         $matepair->{$rd1}{$rd2}{'is'} = $frag_dist;
         #$matepair->{$rd1}{$rd2}{'rd'}{$head}++;  ### this will be used to track uniqueness in pairing (should expect 1 pair per [long]read and one [long]read with that specific pair. HAS LITTLE TO NO EFFECT BUT REQ MORE MEM
         $matepair->{$rd1}{$rd2}{'bt'} = 0;
         $pairct++;
      }

   }
   return $matepair,$pairct;
}

#----------------
sub readFastaFastqBFoff{

   my ($file_ct,$ct_fof_line,$ctrd,$file,$frag_dist,$k,$step,$matepair,$delta,$initpos,$readlength) = @_;
   ###
   my $readinput_message = "\nReads processed from file $file_ct/$ct_fof_line, $file:\n";
   print $readinput_message;
   print LOG $readinput_message;

   ###Support for compressed files MAR2016
   if($file=~/zip$/i){
      open(IN,"unzip -p $file|") || die "Error reading $file -- fatal\n";
   }elsif($file=~/gz$/i || $file=~/gzip$/i){
      open(IN,"gunzip -c $file|") || die "Error reading $file -- fatal\n";
   }else{
      open(IN,$file) || die "Error reading $file -- fatal\n";
   }

   my $prevhead = "";
   my $seq = "";
   my $pairct = 0;
   my $quad = 0;
   my $endposition = $readlength-$k;###MPET
   my $buffer = $readlength+1;###MPET

   LINE:
   while(<IN>){
      chomp;
      $quad++;
      #if(/^([ATGC]+\:[ATGC]+)$/i && $quad<4){$quad=2;}###MPET IGNORE Ns in SEQUENCE (don't have to)
      if(/^([ACGTMRWSYKVHDBN\:]+)$/i && $quad<4){$quad=2;} ###I am sure this will work too. MPET
      if($quad==1 || $quad==3){
         #print "1||3 $_\n";
         next LINE if($_ eq "+");
         my $head = $1 if(/^\S{1}(\S+)/);
#print "$quad **$_** $prevhead ne  && $head ne $prevhead\n";
         if($prevhead ne "" && $head ne $prevhead){# && length($seq)>=50){
            #my $qclipflag = 0;
            #if($qual ne ""){$qclipflag = $r_clip;}
            #print "loadSeq: $prev with $seq and $qual\n";
#print ">$prevhead QUAD:$quad\n$seq\n";

            if(! $readlength){###NO MPET
               $endposition = (length($seq)-$frag_dist);
               $buffer=$k+$delta;
            }else{
               $readlength=length($1) if($seq=~/^([ACGTMRWSYKVHDBN]+)/i);###MPET
               $endposition = $readlength-$k;###MPET
               $buffer = $readlength+1;###MPET
            }
            ($matepair,$pairct) = &kmerizeBFoff(uc($seq),$frag_dist,$k,$matepair,$step,$prevhead,$endposition,$initpos,$pairct,$buffer,$readlength);
            $seq = "";
            $quad=1;
         }
         $ctrd++;
         print "\r$ctrd";
         $|++;
         $prevhead = $head;
      }elsif($quad==2){
         #print "SEQ $_\n";
         $seq .= $_;
      }elsif($quad==4){
         $quad=0;
      }
   }
   if($prevhead ne ""){# && length($seq)>=50){
      if(! $readlength){###NO MPET
        $endposition = (length($seq)-$frag_dist);
        $buffer=$k+$delta;
      }else{
        $readlength=length($1) if($seq=~/^([ACGTMRWSYKVHDBN]+)/i);###MPET
        $endposition = $readlength-$k;###MPET
        $buffer = $readlength+1;###MPET
      }###MPET
      ($matepair,$pairct) = &kmerizeBFoff(uc($seq),$frag_dist,$k,$matepair,$step,$prevhead,$endposition,$initpos,$pairct,$buffer,$readlength);
   }

   close IN;
#   print LOG "$ctrd\n\n";

   #####

   #print "\nExtracted $pairct total $k-mer pairs\n";
   #print LOG "\nExtracted $pairct total $k-mer pairs\n";

   return $matepair,$pairct,$ctrd;
}

#----------------
sub kmerizeBFoff{

   my ($seq,$frag_dist,$k,$matepair,$step,$head,$endposition,$initpos,$pairct,$buffer,$readlength) = @_;

    for(my $pos=$initpos;$pos<=$endposition;$pos+=$step){###MPET

      my $rd1 = substr($seq,$pos,$k);
      $rd1 = &reverseComplement($rd1) if($readlength);###MPET
      my $secondstart = $pos + $buffer;###MPET
      my $rd2ss = substr($seq,$secondstart,$k);
      my $rd2 = &reverseComplement($rd2ss);

      if(defined $matepair->{$rd2}{$rd1}){my $trd1=$rd2;my $trd2=$rd1;$rd1=$trd1;$rd2=$trd2;}

      $matepair->{$rd1}{$rd2}{'is'} = $frag_dist;
      #$matepair->{$rd1}{$rd2}{'rd'}{$head}++;  ### this will be used to track uniqueness in pairing (should expect 1 pair per [long]read and one [long]read with that specific pair. HAS LITTLE TO NO EFFECT BUT REQ MORE MEM
      $matepair->{$rd1}{$rd2}{'bt'} = 0;
      $pairct++;
   }
   return $matepair,$pairct;
}

#----------------
sub kmerizeContigBloom_newloop{
    my ($seq,$bloom,$hashfct,$k) = @_;

    BloomFilter::insertSeq($bloom, $seq, $hashfct, $k);
    return $bloom
}

#----------------
sub kmerizeContigBloom{
   my ($seq,$bloom,$k,$head,$rc) = @_;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $bloom->add($kmer);
   }
   return $bloom;
}

#----------------
sub kmerizeContig{

   my ($seq,$track_all,$matepair,$k,$head,$rc) = @_;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $rd = substr($seq,$pos,$k);
      if(defined $matepair->{$rd}){
         $track_all->{$rd}{'tig'}   = $head;
         $track_all->{$rd}{'start'} = $pos;
         $track_all->{$rd}{'end'}   = $pos + $k;
         if($rc==1){
            $track_all->{$rd}{'start'} = length($seq) - $pos;
            $track_all->{$rd}{'end'}   = length($seq) - ($pos + $k);
         }
         $track_all->{$rd}{'multiple'}++;
      }
   }
   return $track_all;
}

#-----------------------
sub buildGraph{###function builds a basic graph of connections between sequences, highlighting those merged by LINKS

   my ($simplepair,$tl,$graph,$version,$basename,$date) = @_;

   open(GV,">$graph") || die "Can't write to $graph -- fatal\n";
   print GV "graph LINKS{\n";
   print GV "\tlabel=\"LINKS $version $basename $date;\"\n\trankdir=LR;\n\tnode [shape = circle];\n";
   foreach my $tig1(sort {$a<=>$b} keys %$simplepair){
      my $list=$simplepair->{$tig1};
      foreach my $tig2(sort {$a<=>$b} keys %$list){
         if($simplepair->{$tig1}{$tig2}{'links'} >= $min_links){ ### will only output a graph representative of the -l threshold imposed
            my $t1 = $tig1;
            my $t2 = $tig2;
            my $s1 = sprintf "%.1f", ($tl->{$tig1}/1000);
            my $s2 = sprintf "%.1f", ($tl->{$tig2}/1000);
            if($simplepair->{$tig1}{$tig2}{'linked'}){
               my $label = $tig1 . "_" . $s1 . "kb - $tig2" . "_" . $s2 . "kb l=$list->{$tig2}{'links'} g=$simplepair->{$tig1}{$tig2}{'mean'} type=$list->{$tig2}{'type'}";
               print GV "\t$t1 [style=filled, fillcolor=deepskyblue, color=deepskyblue]\n\t$t2 [style=filled, fillcolor=deepskyblue, color=deepskyblue]\n";
               print GV "\t$t1 -- $t2 [ label = \"$label\", penwidth=2.0, color=deepskyblue ]\n";
            }else{
               my $average = sprintf "%.1f", ($list->{$tig2}{'sumgaps'} / $list->{$tig2}{'links'}) if($list->{$tig2}{'links'});
               my $label = $tig1 . "_" . $s1 . "kb - $tig2" . "_" . $s2 . "kb l=$list->{$tig2}{'links'}, g=$average, type=$list->{$tig2}{'type'}";
               print GV "\t$t1 -- $t2 [ label = \"$label\" ]\n";
            }
        }
     }
   }
   print GV "}\n";
   close GV;

}

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/; 
   return (reverse());
}

#------------------------------------
#Order and orient contigs into scaffolds
sub buildScaffolds{

   my ($pair, $tig_length, $simplepair, $scaffold, $longfile, $verbose) = @_;

   open (SC, ">$scaffold") || die "\nCan't write to $scaffold -- fatal\n";
   my $seen_it;
   my $sc_ct = 0;
 
   #print SC "Scaffold Number,Scaffold Size (only contig lengths considered),Scaffold Chain: e.g. _f127z7068k12a0.58m42_r3090z62k7r0.14m76_  means: contig127(+ strand=f), size 7068 (z) has 12 links (k), link ratio of 0.58 (a) and with a mean gap/overlap of 42nt (m)  with reverse (r) of contig3090 (size 62) on the right.\n";

   SEED:
   foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){
      my $ftig = "f" . $tig;
      my $rtig = "r" . $tig;

      if(! defined $seen_it->{$tig}){##should prevent re-using a contig as seed if it's already been incorporated into a scaffold

         $sc_ct++;

         my $chainleft = "";
          
         my $ori_chainright = $ftig . "Z" . $tig_length->{$tig};
         my $chainright = $ori_chainright;
         my $total = $tig_length->{$tig};
         ($total, $chainright, $seen_it, $simplepair) = &computeLayout("R", $chainright, $ftig, $pair, $tig_length, $total, $seen_it, $tig, $simplepair, $longfile);
         ($total, $chainleft, $seen_it, $simplepair) = &computeLayout("L", $chainleft, $rtig, $pair, $tig_length, $total, $seen_it, $tig, $simplepair, $longfile);
         delete $pair->{$ftig};
         delete $pair->{$rtig};
         #delete $tig_length->{$tig};
         $seen_it->{$tig}++;  ### this code was commented out in LINKS v1.8.1 and v1.8.2, resulting in sequence overuse (bug)

         my $scaffold = $chainleft . $chainright;
         print SC "scaffold" . $sc_ct . ",$total,$scaffold\n";
      }
   }
   close SC;
   return $simplepair;
}

#------------------------------------
# links contigs together into a chain - must satisfy user-defined criterions (-k -a)
sub computeLayout{

   my ($ext, $chain, $tig, $pair, $tig_length, $total, $seen_it, $orig_tig_number, $simplepair, $longfile) = @_;

   my $orig_tig = $tig;
   my $extension = 1;

   EXTENSION:
   while($extension){

      my $tnum = $1 if($tig=~/[fr](\d+)/);
      my $tnumf = "f" . $tnum;
      my $tnumr = "r" . $tnum;
      my $priflag = 1;

      if(defined $pair->{$tig}{'layer'}){## XX
 
         #$verbose=1 if($tnum==67961);     
            
         my $szlist = $pair->{$tig}{'layer'};      

         if(defined $pair->{$tig}{'layer'}{10} && -z $longfile){### indicates an ambiguous ARC/KS gap
            my $neighborList = $pair->{$tig}{'layer'}{10};
            my $ambiguous=0;
            SEARCH:
            foreach my $neighbor(sort {$neighborList->{$b}{'links'}<=>$neighborList->{$a}{'links'}} keys %$neighborList){
               $ambiguous=1 if($neighborList->{$neighbor}{'links'} >= $min_links);
               last SEARCH if($ambiguous || $neighborList->{$neighbor}{'links'} < 2);###only need one ambiguous gap est with enough support
            }
            $szlist = $pair->{$tig}{'all'} if($ambiguous);
            $priflag=0;
         }
         #$verbose =1 if($priflag);        
         print "\nVERTEX FORK -- PRIORITYFLAG=$priflag VERBOSE=$verbose\n" if($verbose);
 
         DISTANCE:
         foreach my $distance (sort {$a<=>$b} keys %$szlist){
            print "Attempt to extend $tig, at $distance distance range\n" if ($verbose);
            my $list = $szlist->{$distance};
            my ($match1,$link1,$gaps1,$match2,$link2,$gaps2,$cntloop)=("",0,0,"",0,0,0);

            LINK:
            foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
               if($cntloop){
                  ($match2,$link2,$gaps2) = ($match,$list->{$match}{'links'},$list->{$match}{'gaps'});
                  print "$tig links second best $match2 (links:$link2 total sz:$gaps2)\n" if ($verbose);
                  last LINK;
               }else{
                  ($match1,$link1,$gaps1) = ($match,$list->{$match}{'links'},$list->{$match}{'gaps'});
                  print "$tig links best $match1 (links:$link1 total sz:$gaps1)\n" if ($verbose);
               }
               $cntloop++;
            }

            ###ratio
            my $ratio = 0.00;
            $ratio = $link2 / $link1 if ($link1);        ## relative ratio of the two most abundant contig pairs
            if ($ratio =~ /(\d+\.\d{2})/){$ratio = $1;}
            ###mean
            my $mean = 0;
            $mean = int($gaps1 / $link1) if ($link1);

            my $tempnum = $1 if($match1 =~ /[fr](\d+)/);
            my $tempnumf = "f" . $tempnum;

            #### Assessment
            if(defined $seen_it->{$tempnum} || $link1 < $min_links || $ratio > $max_link_ratio || $tempnum == $orig_tig_number){
               $extension = 0;
               print "defined seen_it->{ $tempnum } || $link1 < $min_links || $ratio > $max_link_ratio\n L1:$link1 L2:$link2  M1:$match1 M2:$match2 G1:$gaps1 G2:$gaps2 "  if ($verbose);

               next DISTANCE; ### CONDITIONS NOT SATISFIED, CHECK NEXT DISTANCE
            }{### pass filter.. does this contig 
               print "$ext extension.  mean: $mean links:$link1 linkratio:$ratio\n" if ($verbose);
               my ($order1,$order2)= $tnum < $tempnum ? ($tnum,$tempnum) : ($tempnum,$tnum);
               $simplepair->{$order1}{$order2}{'linked'}=1;
               $simplepair->{$order1}{$order2}{'mean'}=$mean;

               if($ext eq "R"){
                  $chain .= "k" . $link1 . "a" . $ratio . "m" . $mean . "_" . $match1 . "z" . $tig_length->{$tempnum};
               }else{
                  my $temp_match = "";
                  if($match1 =~ /^r(\d+)/){$temp_match = "f" . $1;}else{$temp_match = "r". $1;}            
                  $chain = $temp_match . "z" . $tig_length->{$tempnum} . "k" . $link1 . "a" . $ratio . "m" . $mean . "_" . $chain;
               }   
               $total += $tig_length->{$tempnum};

               print "NEXT TIG TO LOOK AT= $match1\n" if ($verbose);
               $tig = $match1;
               $extension = 1; 

               $seen_it->{$tnum} = 1;
               $seen_it->{$tempnum} = 1;          

               print "Will flag $tnum as seen  (only if $tnumf != $orig_tig)." if ($verbose);
   
               #if($tnumf ne $orig_tig){
               #   delete $pair->{$tnumf}{$distance};
               #   delete $pair->{$tnumr}{$distance};
               #   delete $tig_length->{$tnum};
               #}else{
               #   delete $pair->{$tnumf}{$distance};
               #}
               last DISTANCE;#### MOVE TO NEXT CONTIG
            }
         }### cycle through short to long gaps
      }else{
         print "NO MORE MATCH FOR $tig in hash: pair>>\n" if ($verbose);
         $extension = 0;
         last EXTENSION;
      }
   }### pair is defined
   return $total, $chain, $seen_it, $simplepair;
}

#------------------------------------
sub getDistance{

   my ($insert_size, $length_i, $start_i, $start_j) = @_;

   # L  ------  --------- R
   # i    ->        <-    j
   #      ....  ......    insert_span
   #      ============    insert_size

   my $insert_span = ($length_i - $start_i) + $start_j;
   my $gap_or_overlap = $insert_size - $insert_span;

   return $gap_or_overlap;
}

#-----------------
#build contig pairs based on template information  -  must satisfy user-defined criterions (-d -e)
sub pairContigs{

   my ($longfile,$matepair,$track,$tig_length,$issues,$distribution,$totalpairs,$tigpair_checkpoint,$simplepair_checkpoint,$verbose) = @_;
   my ($ct_illogical, $ct_ok_contig, $ct_ok_pairs, $ct_problem_pairs, $ct_iz_issues, $ct_single, $ct_both)= (0,0,0,0,0,0,0);
   my $ct_illogical_hash;
   my $ct_ok_contig_hash;
   my $ct_ok_pairs_hash;
   my $ct_problem_pairs_hash;
   my $ct_iz_issues_hash;
   my $ct_single_hash;
   my $ct_both_hash;
   my $ct_multiple;

   my ($simplepair,$pair,$err,$track_insert);###simplepair is a simple pair to track node to node in a graph

   print "Pairing contigs...\n" if ($verbose);

   open(PET, ">$issues") || die "Can't open $issues for writing -- fatal\n";

   foreach my $read_a (keys %$matepair){ 

      my $mateslist = $matepair->{$read_a};

      foreach my $read_b (keys %$mateslist){

          if($matepair->{$read_a}{$read_b}{'bt'}==0 && $track->{$read_a}{'multiple'}==1 && $track->{$read_b}{'multiple'}==1){ ###This has little if no effect, but negative for some odd reason

            ##below indicates this specific pair has been seen
            $matepair->{$read_a}{$read_b}{'bt'}=1;

            my $insert_size = $mateslist->{$read_b}{'is'};
            my $min_allowed = -1 * ($insert_stdev * $insert_size);
            my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);

            print "Pair read1=$read_a read2=$read_b\n" if ($verbose);

            if(defined $track->{$read_a}{'tig'} && defined $track->{$read_b}{'tig'}){### both pairs assembled

               $ct_both++;
               $ct_both_hash->{$insert_size}++;

               my $tig_a = $track->{$read_a}{'tig'};
               my $tig_b = $track->{$read_b}{'tig'};

               my $ftig_a = "f" . $tig_a;
               my $ftig_b = "f" . $tig_b;

               my $rtig_a = "r" . $tig_a;
               my $rtig_b = "r" . $tig_b;

               my $A_length = $tig_length->{$tig_a};
               my $A_start = $track->{$read_a}{'start'};
               my $A_end = $track->{$read_a}{'end'};
 
               my $B_length = $tig_length->{$tig_b};
               my $B_start = $track->{$read_b}{'start'} ;
               my $B_end = $track->{$read_b}{'end'};

               if ($tig_a != $tig_b){####paired reads located on <> contigs

                  ####Determine most likely possibility
                  if ($track->{$read_a}{'start'} < $track->{$read_a}{'end'}){

                     if ($track->{$read_b}{'end'} < $track->{$read_b}{'start'}){####-> <- :::  A-> <-B  /  rB -> <- rA
                         my $d = &getDistance($insert_size, $A_length, $A_start, $B_start);
                         print "A-> <-B  WITH $tig_a -> <- $tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen, Astart,Bstart\n" if($verbose);
                         if($d >= $min_allowed){

                            my $isz = $d < 0 ? -1 : $d == 10 ? 10 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 5000 ? 5000 : 10000;###distance categories
                            #my $isz = $d < 0 ? -1 : $d < 200 ? 200 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 2500 ? 2500 : $d < 5000 ? 5000 : $d < 10000 ? 10000 : 20000;###distance categories
                            $pair->{$ftig_a}{'layer'}{$isz}{$ftig_b}{'links'}++;
                            $pair->{$ftig_a}{'layer'}{$isz}{$ftig_b}{'gaps'} += $d;                  
                            $pair->{$rtig_b}{'layer'}{$isz}{$rtig_a}{'links'}++;
                            $pair->{$rtig_b}{'layer'}{$isz}{$rtig_a}{'gaps'} += $d;
                            $pair->{$ftig_a}{'all'}{0}{$ftig_b}{'links'}++;
                            $pair->{$ftig_a}{'all'}{0}{$ftig_b}{'gaps'} += $d;
                            $pair->{$rtig_b}{'all'}{0}{$rtig_a}{'links'}++;
                            $pair->{$rtig_b}{'all'}{0}{$rtig_a}{'gaps'} += $d;

                            my ($order1,$order2)=($tig_a<$tig_b) ? ($tig_a,$tig_b) : ($tig_b,$tig_a);
                            $simplepair->{$order1}{$order2}{'links'}++;
                            $simplepair->{$order1}{$order2}{'sumgaps'}+=$d;                            
                            $simplepair->{$order1}{$order2}{'type'}='11';

                            $ct_ok_pairs++;
                            $ct_ok_pairs_hash->{$insert_size}++;
                         }else{
                            my $err_pair = $ftig_a . "-". $ftig_b;
                            $err->{$err_pair}{'links'}++;
                            $err->{$err_pair}{'gaps'} += $d;
                            $ct_problem_pairs++;
                            $ct_problem_pairs_hash->{$insert_size}++;
                            print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#$tig_a -> $d <- tig#$tig_b, A=$A_length nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                         }
                      }else{#### -> -> ::: A-> <-rB  / B-> <-rA 
                         my $rB_start = $B_length - $B_start;
                         my $d = &getDistance($insert_size, $A_length, $A_start, $rB_start);
                         print "A-> <-rB  WITH $tig_a -> <- r.$tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen,Astart,rBstart\n" if($verbose);
                         if($d >= $min_allowed){

                            my $isz = $d < 0 ? -1 : $d == 10 ? 10 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 5000 ? 5000 : 10000;###distance categories
                            #my $isz = $d < 0 ? -1 : $d < 200 ? 200 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 2500 ? 2500 : $d < 5000 ? 5000 : $d < 10000 ? 10000 : 20000;###distance categories
                            $pair->{$ftig_a}{'layer'}{$isz}{$rtig_b}{'links'}++;
                            $pair->{$ftig_a}{'layer'}{$isz}{$rtig_b}{'gaps'} += $d;
                            $pair->{$ftig_b}{'layer'}{$isz}{$rtig_a}{'links'}++;
                            $pair->{$ftig_b}{'layer'}{$isz}{$rtig_a}{'gaps'} += $d;
                            $pair->{$ftig_a}{'all'}{0}{$rtig_b}{'links'}++;
                            $pair->{$ftig_a}{'all'}{0}{$rtig_b}{'gaps'} += $d;
                            $pair->{$ftig_b}{'all'}{0}{$rtig_a}{'links'}++;
                            $pair->{$ftig_b}{'all'}{0}{$rtig_a}{'gaps'} += $d;
               
                            my ($order1,$order2)=($tig_a<$tig_b) ? ($tig_a,$tig_b) : ($tig_b,$tig_a);
                            $simplepair->{$order1}{$order2}{'links'}++;
                            $simplepair->{$order1}{$order2}{'sumgaps'}+=$d;
                            $simplepair->{$order1}{$order2}{'type'}='10';

                            $ct_ok_pairs++;
                            $ct_ok_pairs_hash->{$insert_size}++;
                         }else{
                            my $err_pair = $ftig_a . "-". $rtig_b;
                            $err->{$err_pair}{'links'}++;
                            $err->{$err_pair}{'gaps'} += $d;
                            $ct_problem_pairs++;
                            $ct_problem_pairs_hash->{$insert_size}++;
                            print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-rB  WITH tig#$tig_a -> $d <- tig#r.$tig_b, A=$A_length  nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                         }
                      }
                  }else{

                     if ($track->{$read_b}{'end'} > $track->{$read_b}{'start'}){####<-  -> ::: B-> <-A / rA -> <- rB
                        my $d = &getDistance($insert_size, $B_length, $B_start, $A_start);
                        print "B-> <-A  WITH $tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,Bstart,Astart\n" if($verbose);
                        if($d >= $min_allowed){

                           my $isz = $d < 0 ? -1 : $d == 10 ? 10 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 5000 ? 5000 : 10000;###distance categories
                           #my $isz = $d < 0 ? -1 : $d < 200 ? 200 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 2500 ? 2500 : $d < 5000 ? 5000 : $d < 10000 ? 10000 : 20000;###distance categories
                           $pair->{$ftig_b}{'layer'}{$isz}{$ftig_a}{'links'}++;
                           $pair->{$ftig_b}{'layer'}{$isz}{$ftig_a}{'gaps'} += $d;
                           $pair->{$rtig_a}{'layer'}{$isz}{$rtig_b}{'links'}++;
                           $pair->{$rtig_a}{'layer'}{$isz}{$rtig_b}{'gaps'} += $d;
                           $pair->{$ftig_b}{'all'}{0}{$ftig_a}{'links'}++;
                           $pair->{$ftig_b}{'all'}{0}{$ftig_a}{'gaps'} += $d;
                           $pair->{$rtig_a}{'all'}{0}{$rtig_b}{'links'}++;
                           $pair->{$rtig_a}{'all'}{0}{$rtig_b}{'gaps'} += $d;

                           my ($order1,$order2)=($tig_a<$tig_b) ? ($tig_a,$tig_b) : ($tig_b,$tig_a); 
                           $simplepair->{$order1}{$order2}{'links'}++;
                           $simplepair->{$order1}{$order2}{'sumgaps'}+=$d;
                           $simplepair->{$order1}{$order2}{'type'}='11';

                           $ct_ok_pairs++;
                           $ct_ok_pairs_hash->{$insert_size}++;
                        }else{
                           my $err_pair = $ftig_b . "-". $ftig_a;
                           $err->{$err_pair}{'links'}++;
                           $err->{$err_pair}{'gaps'} += $d;
                           $ct_problem_pairs++;
                           $ct_problem_pairs_hash->{$insert_size}++;
                           print PET "Pairs unsatisfied in distance within a contig pair.  B-> <-A  WITH tig#$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                        }
                     }else{                          ####<- <-  :::  rB-> <-A / rA-> <-B
                        my $rB_start = $B_length - $B_start;
                        my $d = &getDistance($insert_size, $B_length, $rB_start, $A_start);
                        print "rB-> <-A WITH r.$tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,rBstart,Astart\n" if($verbose);
                        if($d >= $min_allowed){

                           my $isz = $d < 0 ? -1 : $d == 10 ? 10 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 5000 ? 5000 : 10000;###distance categories
                           #my $isz = $d < 0 ? -1 : $d < 200 ? 200 : $d < 500 ? 500 : $d < 1000 ? 1000 : $d < 2500 ? 2500 : $d < 5000 ? 5000 : $d < 10000 ? 10000 : 20000;###distance categories

                           $pair->{$rtig_b}{'layer'}{$isz}{$ftig_a}{'links'}++;
                           $pair->{$rtig_b}{'layer'}{$isz}{$ftig_a}{'gaps'} += $d;
                           $pair->{$rtig_a}{'layer'}{$isz}{$ftig_b}{'links'}++;
                           $pair->{$rtig_a}{'layer'}{$isz}{$ftig_b}{'gaps'} += $d;
                           $pair->{$rtig_b}{'all'}{0}{$ftig_a}{'links'}++;
                           $pair->{$rtig_b}{'all'}{0}{$ftig_a}{'gaps'} += $d;
                           $pair->{$rtig_a}{'all'}{0}{$ftig_b}{'links'}++;
                           $pair->{$rtig_a}{'all'}{0}{$ftig_b}{'gaps'} += $d;

                           my ($order1,$order2)=($tig_a<$tig_b) ? ($tig_a,$tig_b) : ($tig_b,$tig_a);
                           $simplepair->{$order1}{$order2}{'links'}++;
                           $simplepair->{$order1}{$order2}{'sumgaps'}+=$d;
                           $simplepair->{$order1}{$order2}{'type'}='01';

                           $ct_ok_pairs++;
                           $ct_ok_pairs_hash->{$insert_size}++;
                        }else{
                           my $err_pair = $rtig_b . "-". $ftig_a;
                           $err->{$err_pair}{'links'}++;
                           $err->{$err_pair}{'gaps'} += $d;
                           $ct_problem_pairs++;
                           $ct_problem_pairs_hash->{$insert_size}++;
                           print PET "Pairs unsatisfied in distance within a contig pair.  rB-> <-A WITH tig#r.$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                        }
                     }
                  }
               }else{###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
           
                  print "Pair ($read_a and $read_b) located on same contig $tig_a ($A_length nt)\n" if ($verbose);
                  my $pet_size = 0;

                  if ($A_start > $B_start && ($B_start < $B_end) && ($A_start > $A_end)){    # B --> <-- A
                     $pet_size = $A_start - $B_start;
                     $track_insert->{$pet_size}++;
                     if($pet_size >= $low_iz && $pet_size <= $up_iz){
                        $ct_ok_contig++;
                        $ct_ok_contig_hash->{$insert_size}++;
                     }else{
                        print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                        $ct_iz_issues++;
                        $ct_iz_issues_hash->{$insert_size}++;
                     }
                  }elsif($B_start > $A_start && ($B_start > $B_end) && ($A_start < $A_end)){ # A --> <-- B
                     $pet_size = $B_start - $A_start;
                     $track_insert->{$pet_size}++;
                     if($pet_size >= $low_iz && $pet_size <= $up_iz){
                        $ct_ok_contig++;
                        $ct_ok_contig_hash->{$insert_size}++;
                     }else{
                        print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                        $ct_iz_issues++;
                        $ct_iz_issues_hash->{$insert_size}++;
                     }
                  }else{
                     $ct_illogical++;
                     $ct_illogical_hash->{$insert_size}++;
                     print PET "Pairs unsatisfied in pairing logic within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end\n";
                  }
               }
            }else{###^both pairs assembled
               $ct_single++;
               $ct_single_hash->{$insert_size}++;
            }
         }else{#if unseen
            $ct_multiple++ if( $matepair->{$read_a}{$read_b}{'bt'}==0 );
	 }
     }#pairing read b
   }#read a

   ### summary of contig pair issues
   print PET "------------- Putative issues with contig pairing - Summary  ----------------\n";
   foreach my $err_pair (sort {$err->{$b}{'links'}<=>$err->{$a}{'links'}} keys %$err){
      my $mean_iz = 0;
      $mean_iz = $err->{$err_pair}{'gaps'} / $err->{$err_pair}{'links'} if ($err->{$err_pair}{'links'});
      print PET "Pair $err_pair has $err->{$err_pair}{'links'} links and mean distance = $mean_iz\n";
   }
   close PET;
 
   my $satisfied = $ct_ok_pairs + $ct_ok_contig;
   my $unsatisfied = $ct_problem_pairs + $ct_iz_issues + $ct_illogical;
   my $ct_both_reads = $ct_both * 2;

   print LOG "\n===========PAIRED K-MER STATS===========\n";
   print LOG "Total number of pairs extracted from -s $longfile: $totalpairs\n";
   print LOG "At least one sequence/pair missing from contigs: $ct_single\n";
   print LOG "Ambiguous kmer pairs (both kmers are ambiguous): $ct_multiple\n";
   print LOG "Assembled pairs: $ct_both ($ct_both_reads sequences)\n";
   print LOG "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $ct_ok_contig\n";
   print LOG "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues\n";
   print LOG "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical\n";
   print LOG "\t---\n";
   print LOG "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs\n";
   print LOG "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs\n";
   print LOG "\t---\n";
   print LOG "Total satisfied: $satisfied\tunsatisfied: $unsatisfied\n\nBreakdown by distances (-d):\n";

   foreach my $izopt(sort {$a<=>$b} keys %$ct_both_hash){
      print LOG "--------k-mers separated by $izopt bp (outer distance)--------\n";
      my $maopt = -1 * ($insert_stdev * $izopt);
      my ($low_izopt, $up_izopt) = ($izopt + $maopt, $izopt - $maopt);
      print LOG "MIN:$low_izopt MAX:$up_izopt as defined by $izopt * $insert_stdev\n";
      print LOG "At least one sequence/pair missing: $ct_single_hash->{$izopt}\n";
      print LOG "Assembled pairs: $ct_both_hash->{$izopt}\n";
      print LOG "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $ct_ok_contig_hash->{$izopt}\n";
      print LOG "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues_hash->{$izopt}\n";
      print LOG "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical_hash->{$izopt}\n";
      print LOG "\t---\n";
      print LOG "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs_hash->{$izopt}\n";
      print LOG "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs_hash->{$izopt}\n";
   }
   print LOG "============================================\n";

   open (CSV, ">$distribution") || die "Can't open $distribution for writing -- fatal";

   foreach my $is (sort {$a<=>$b} keys %$track_insert){
      print CSV "$is,$track_insert->{$is}\n";
   }
   close CSV;

   ####mid-scaffolding checkpoint
   ####Save datastructure pair to disk to recover quickly in case of crash
   open(CP,">$tigpair_checkpoint") || die "Can't write to $tigpair_checkpoint -- fatal.\n";
   foreach my $tig1 (keys %$pair){
      my $distlist = $pair->{$tig1}{'layer'};
      foreach my $distance (keys %$distlist){
         my $secpair=$distlist->{$distance};
         foreach my $tig2 (keys %$secpair){
           print CP "$distance\t$tig1\t$tig2\t$secpair->{$tig2}{'links'}\t$secpair->{$tig2}{'gaps'}\n";
         }
      }
   }
   print CP "~~~end~~~\n";
   close CP;

   open(SP,">$simplepair_checkpoint") || die "Can't write to $simplepair_checkpoint -- fatal.\n";
   foreach my $tig1(keys %$simplepair){
      my $secpair=$simplepair->{$tig1};
      foreach my $tig2(keys %$secpair){
        print SP "$tig1\t$tig2\t$secpair->{$tig2}{'links'}\t$secpair->{$tig2}{'sumgaps'}\t$secpair->{$tig2}{'type'}\n";
      }
   }
   print SP "~~~end~~~\n";
   close SP;

   return $pair,$simplepair;
}

#----------------
sub readContigPairs{

   my ($contigpairfile,$simplepairfile,$tig_length) = @_;

   my $pair;
   my $simplepair;
   my $simplepair_flag = 1;
   if(-e $simplepairfile){
      ###Check if file exists. If it does, don't want to override its content
      $simplepair_flag = 0;
   }

   open(IN,$contigpairfile) || die "Can't read $contigpairfile -- fatal.\n";
   while(<IN>){
      chomp;
      my @info=split(/\t/);

      # Skip contigs shorter than -z min_size.
      my $uid = substr $info[1], 1;
      my $vid = substr $info[2], 1;
      next if $tig_length->{$uid} < $min_size || $tig_length->{$vid} < $min_size;

      ### information organized as such:
      ### distance	tig1    tig2    numLINKS    gapSize 
      ### 500     r1000005        r1856837        1       10
      if($info[2]){
         $pair->{$info[1]}{'layer'}{$info[0]}{$info[2]}{'links'}+=$info[3];
         $pair->{$info[1]}{'layer'}{$info[0]}{$info[2]}{'gaps'}+=$info[4];
         $pair->{$info[1]}{'all'}{0}{$info[2]}{'links'}+=$info[3];
         $pair->{$info[1]}{'all'}{0}{$info[2]}{'gaps'}+=$info[4];

	 ###added FEB2017
         if($simplepair_flag){### will only proceed to create this data structure if the simplepair file is missing

            my ($orDir1,$orTig1) = ($1,$2) if($info[1] =~ /^([rf])(\d+)$/);
	    my ($orDir2,$orTig2) = ($1,$2) if($info[2] =~ /^([rf])(\d+)$/);
            my ($stripTig1,$stripTig2,$dir1,$dir2) = ($orTig1,$orTig2,$orDir1,$orDir2);

	    if($orTig1<$orTig2){
		    #$stripTig1 = $orTig2;
		    #   $stripTig2 = $orTig1;
		    #   $dir1 = $orDir2;
		    #   $dir2 = $orDir1;
		    #}
	       $simplepair->{$stripTig1}{$stripTig2}{'links'} += $info[3];
               $simplepair->{$stripTig1}{$stripTig2}{'sumgaps'} += $info[4];

	       my $oritype = "";
               if($dir1 eq "f"){
	   	 $oritype .= "1";
                 if($dir2 eq "f"){
                    $oritype .= "1";
                 }else{
                    $oritype .= "0";
                 }
               }else{
                  $oritype .= "0";
                    if($dir2 eq "f"){
                       $oritype .= "1";
                    }else{
                       $oritype .= "0";
                    }
               }
	       if ($info[3] >= $simplepair->{$stripTig1}{$stripTig2}{'max'}){
                  $simplepair->{$stripTig1}{$stripTig2}{'type'} = $oritype;###direction with most support
	          $simplepair->{$stripTig1}{$stripTig2}{'max'} = $info[3];
               }
            }
         }
      }
   }
   close IN;

   open(SPR,$simplepairfile);###don't die if doesn't exists, not crucial file
   while(<SPR>){
      chomp;
      my @info=split(/\t/);
      ### information organized as such:
      ### tig1    tig2    numLINKS    sumgap   type
      ### 32      51      1       8036    10
      if($info[2]){
         $simplepair->{$info[0]}{$info[1]}{'links'}=$info[2];
         $simplepair->{$info[0]}{$info[1]}{'sumgaps'}=$info[3];
         $simplepair->{$info[0]}{$info[1]}{'type'}=$info[4];
      }
   }
   close SPR;

   return $pair, $simplepair;
}

#-------------------
sub readContigsMemory{
 
   my $file = shift;

   my $tig_length;
   my $tignames;
   my $fh;
   my $prevhead="";
   my $seq="";
   my $cttig=0;

   ###Support for compressed files MAR2016
   if($file=~/zip$/i){
      open(IN,"unzip -p $file|") || die "Error reading $file -- fatal\n";
   }elsif($file=~/gz$/i || $file=~/gzip$/i){
      open(IN,"gunzip -c $file|") || die "Error reading $file -- fatal\n";
   }else{
      open(IN,$file) || die "Error reading $file -- fatal\n";
   }

###TO REVIEW BELOW

#   while(<IN>){
#      chomp;
#      if (/\>(\S+)/){
#         $cttig++;
#         my $head=$cttig;
#         #$seq =~ s/[BDEFHIJKLMOPQRSUVWXYZ]/N/g;
#         if($prev ne $head && $prev ne "NA"){
#             $fh->{$prev} = $seq;
#         }
#         $prev = $head;
#         $seq='';
#      }elsif(/^(\S+)$/){
#         $seq .= uc($1);
#      }
#   }
#   $fh->{$prev} = $seq;
#   close IN;

###
   while(<IN>){
      chomp;
      if(/^\>(.*)/){
         my $head=$1;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $cttig++;
            $tignames->{$cttig} = $prevhead;
            $fh->{$cttig} = $seq;
            $tig_length->{$cttig} = length($seq);
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= uc($_);
      }
   }
   $cttig++;
   $tignames->{$cttig} = $prevhead;
   $fh->{$cttig} = $seq;
   $tig_length->{$cttig} = length($seq);
###

   return $fh,$tignames,$tig_length;
}

#-------------------
sub buildScaffoldFasta{

   my ($dotscaffold,$fh,$assemblyfile,$numnamecorr,$tignames) = @_;

   open(IN,$dotscaffold) || die "Cannot open $dotscaffold for reading -- fatal.\n";

   my $scaffold_fasta = $dotscaffold . ".fa";

   open(OUT,">$scaffold_fasta") || die "can't write to $scaffold_fasta -- fatal\n";
   open(CORROUT,">$numnamecorr") || die "Can't write to $numnamecorr -- fatal\n";
   print CORROUT "#This correspondence $numnamecorr file lists the scaffold ID and contig ID used by LINKS (reported in the .scaffolds file), with the original contig name in the provided assembly file (-f $assemblyfile)\n";
   print CORROUT "LINKS_scaffold_ID\tLINKS_contig_ID\toriginal_name\torientation(f=forward/r=reverse)\tnumber_of_links\tlinks_ratio\tgap_or_overlap(-)\n";

   my $tot=0;
   my $ct=0;
   my $sct=0;
   while(<IN>){
      chomp;   
      my $sc="";
      my @a = split(/\,/);
      my @tig;
   
      if($a[2]=~/\_/){
         @tig = split(/\_/,$a[2]);
      }else{
         push @tig, $a[2];
      }

      $sct++;
      my $tigsum=0;
      print OUT ">$_\n";
      foreach my $t (@tig){
         $ct++;
    
         if($t=~/([fr])(\d+)z(\d+)(\S+)?/i){

            my $orient = $1;
            my $tnum=$2;
            my $head = $orient . $tnum;
            my $search = $tnum;
            my $other = $4;
            $tot+= $3;
            $tigsum +=$3;

            my $gap = "NA"; 
            my $numlinks = "NA";
            my $linksratio = "NA";

            my $gapseq = "";

            $numlinks = $1 if($other=~/k(\d+)/);
            $linksratio = $1 if($other=~/a(\d+.*)m/);
	    $gap = $1 if($other=~/m(\-?\d+)/);

            print CORROUT "$a[0]\t$tnum\t$tignames->{$tnum}\t$orient\t$numlinks\t$linksratio\t$gap\n"; 

            my $seq = $fh->{$search};
            $seq = reverseComplement($seq) if($orient eq "r");

            $gapseq = "Z" x $gap if($gap > 0); #!!changed to Z for testing gap sizes March 28
            $gapseq = "n" if($gap ne "NA" && $gap <= 0 );
            $seq .= $gapseq;

            print OUT "$seq";         
                
         }#tig regex
   
      }#each tig
      print OUT "\n";
   }

   close CORROUT;
   close IN;
   close OUT;

   return $scaffold_fasta;
}

## We hope this code is useful to you -- Please send comments & suggestions to rwarren at bcgsc.ca
