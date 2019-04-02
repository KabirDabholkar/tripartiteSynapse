#!/usr/bin/perl -w

# Note: stimuli are 5 ms wide.
# $interval is time between end of one stim to beginning of next.

$delay=0;
$interval=10e-3 - 5e-3;  #note: stim is 5 ms long
$n_repeats=30;

@infile=();
@outfile=();

$infile[1]="VDCC_PQ_C01.dat";
$outfile[1]="VDCC_PQ_C01_100hz_300ms.dat";
$infile[2]="VDCC_PQ_C10.dat";
$outfile[2]="VDCC_PQ_C10_100hz_300ms.dat";
$infile[3]="VDCC_PQ_C12.dat";
$outfile[3]="VDCC_PQ_C12_100hz_300ms.dat";
$infile[4]="VDCC_PQ_C21.dat";
$outfile[4]="VDCC_PQ_C21_100hz_300ms.dat";
$infile[5]="VDCC_PQ_C23.dat";
$outfile[5]="VDCC_PQ_C23_100hz_300ms.dat";
$infile[6]="VDCC_PQ_C32.dat";
$outfile[6]="VDCC_PQ_C32_100hz_300ms.dat";
$infile[7]="VDCC_PQ_C34.dat";
$outfile[7]="VDCC_PQ_C34_100hz_300ms.dat";
$infile[8]="VDCC_PQ_C43.dat";
$outfile[8]="VDCC_PQ_C43_100hz_300ms.dat";
$infile[9]="VDCC_PQ_Ca.dat";
$outfile[9]="VDCC_PQ_Ca_100hz_300ms.dat";

for ($i=1;$i<=$#infile;$i++) {

#################
# Read the stim #
#################
  @x=();
  @y=();
  $filename=$infile[$i];
  open(FILE,"< $filename");
  while(defined($line=<FILE>)) {
    chomp($line);
    @tokens=split(' ',$line);
    $x[$#x+1]=$tokens[0];
    $y[$#y+1]=$tokens[1];
  }
  close(FILE);


    #################
##### Make the stim #####
    #################
  $filename=$outfile[$i];
  open(FILE,"> $filename");
  $val=$y[1];
  printf FILE ("%.15g %.15g\n",0,$val);

  $duration=$x[$#x];
  for ($j=1;$j<=$n_repeats;$j++) {
    $t_offset=$delay+($j-1)*($interval+$duration);
    for ($k=1;$k<=$#x;$k++) {
      printf FILE ("%.15g %.15g\n",$x[$k]+$t_offset,$y[$k]);
    }
  }
  close(FILE);

}

exit(0);

