# Written  by Dan Gusfield,  copyright 2009
#  Permission to use this program comes with no guarantees. You may use it if you don't laugh at it and
#  #  don't expect the author to understand or explain the code in the future.
#  #  Permission is not granted for redistribution or distribution of modified versions. Please contact
#  #  the author concerning distributing any modifications.
#
# July 18, 2008 modify multextract5.pl to work for a variable number of states, k, input in the call.
#multextract.pl
# Dec. 30, 2007
#
#
#
$infile = "$ARGV[0]"; # a binary haplotype file produced by ms 
$prob = "$ARGV[1]"; # this is the probability that a value should be deleted.
$k = "$ARGV[2]"; # this is the maximum number of states allowed in a character.

open (IN, "$infile");
open OUTB, ">$ARGV[0]b";  # where the cleaned up binary data is output 
open (OUTD, '>datalist');  # list of the files for the extracted data, with possible deletions
open (OOUTD, '>odatalist');  # list of the files for the original extracted data, with no deletions
#open (STAT, '>>stats');
$s = 0;

$line = <IN>;
if ($line =~ /ms (\d+) (\d+) -s (\d+) -r (\d+)/){  #extract the header information
                                                     #from data in the header line in the ms file
#   print STAT "$line";
   $numi = $1;
   $numreps = $2;
   $numj = $3;
   $recrate = $4;
   $kthnumj = $numj/($k-1);
   print "$numi, $kthnumj\n";
#   print STAT "$numreps reps of $numi by $kthnumj with r = $recrate:\n";
}

@lines = <IN>;
close(IN);

$i = 0;
foreach $line (@lines) {
  if ($line =~ /segsites/) {
     close OUT;
     $s++;
     print OUTB "\n set $s\n";

     $outprefix = "$k" . 'state';
     $outfile = $outprefix . $numi . '.' . $kthnumj . '.' . $s;
     $originaloutfile = "o$outfile";  # this will hold the state data for one problem instance
     print OUTD "$outfile\n";
     print OOUTD "$originaloutfile\n";
     open (OUT, ">$outfile");
     open (OOUT, ">$originaloutfile");
     print OUT "$s\n";
     print OUT "$numi $kthnumj\n";
     print OOUT "$s\n";
     print OOUT "$numi $kthnumj\n";

     for ($j = 0; $j <= $numj - $k + 1; $j = $j + $k - 1) {
         $next[$j] = 0;
     }

     %bits = (); # clear out the hash bits each time we see a new dataset.
  }

  if ($line =~ /^[01]+$/) {
  $i++;
  $line =~ tr/ //d;
  print OUTB "$line";
  chomp $line;
  @row = split (//,$line);

#  print "processing lines\n";

   $data = "";
   $odata = "";

   for ($j = 0; $j <= $numj - $k + 1; $j = $j + $k - 1) {

      $key = $j;
        for ($i = 0; $i < $k - 1; $i++) {
	  $key = "$key . $row[$j + $i]";
	}

       if (defined $bits{$key}) {
          $value = $bits{$key};
       }
       else {
          $value = $next[$j];  # next[j] is initialized to 0, so the first assigned value for any j is 0.
            if ($value < $k) { 
	    $bits{$key} = $value;
	    $next[$j]++;
	    }
	    else {
	    	$value = int(rand($k));
            }
       }

       $randvalue = rand(1);
 #       print "randomvalue $randvalue\n";
       if ($randvalue <= $prob) {
         $data .= '? ';
       }
       else {
             $data .= "$value ";
       }
      $odata .= "$value ";
   }

    print OUT "$data\n";
    print OOUT "$odata\n";

  }  # end of rowwise processing for a single row of binary data

}
close (OUTB);
close (OUTD);
close (OUT);
close (OOUT);
