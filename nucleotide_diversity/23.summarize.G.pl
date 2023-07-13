#!/usr/bin/perl
use strict;
use warnings;

# summarize the pi values from the pixy runs

# configuration variables
my ( $pfx ) = $0 =~ m/\/(\d+)[^\/]*$/;



# input files
my $tmpdir = "/dauc3/tmp/$pfx.kevinJ";



# output files
my $notebook = '/nas-simonlab/simonlabweb/bioinformatics/up/notebook/0082';



# global variables
my %data;



my $outfilename = "$notebook/$pfx.pixy.G.summary.tsv";
unlink( $outfilename );
for my $windowsize (100000)
  {
    for my $suptbl (33, 39)
      {
        my $label = "$windowsize.$suptbl";
        my $pixydir = "$tmpdir/Gpixy/$label";
        my $title = 'Low Admixture Results, Window Size='.$windowsize;
        if ($suptbl == 33) { $title = 'Full Sample Set Results, Window Size='.$windowsize; }
        %data = ();
        for my $chr (1..9)
          {
            my $pifile = "$pixydir/Chr$chr/pixy_pi.txt";
            loadinput( $pifile );
          }
        saveoutput( $outfilename, $title );
      }
  }



system( "cp -puv --no-preserve=owner $0 $notebook/" );
tlog( "$0 Done" );
exit 0;



###############################################################
sub loadinput { my ( $infilename ) = @_;
###############################################################
  tlog( "Loading \"$infilename\"" );
  my @indata;
  my $nlines = 0;
  my $INF = stdopen ( '<', $infilename );
  while ( my $aline = <$INF> )
    {
      $nlines++;
      next if ($nlines <= 1);
      $aline =~ s/[\r\n]//g;
      my @cols = split ( /\t/, $aline );
#pop     chromosome      window_pos_1    window_pos_2    avg_pi  no_sites        count_diffs     count_comparisons       count_missing
#Improved_Cultivar       DCARv3_Chr1     1       10000   0.0     1154    0       166315475       27452665
#Early_Cultivar  DCARv3_Chr1     1       10000   0.0     1154    0       15268576        2100278
#Landrace-B      DCARv3_Chr1     1       10000   0.0     1154    0       19606608        4288116
      my $pop = $cols[0];

      # pixy recommended method to average is diffs/comps
      if ( ( $cols[6] ne 'NA' ) and ( $cols[7] ne 'NA' ) )
        {
          $data{$pop}{n} ++;
          $data{$pop}{countdiffs} += $cols[6];
          $data{$pop}{countcomps} += $cols[7];
        }
    } # while <$INF>
  stdclose ( $INF );
  tlog( commify($nlines) . " lines processed" );
  return ( @indata );
} # sub loadinput



###############################################################
sub saveoutput { my ( $outfilename, $title ) = @_;
###############################################################
  tlog( "Saving \"$outfilename\"" );
  my $OUTF = stdopen ( '>>', $outfilename );
  print $OUTF "$title\n";
  print $OUTF join( "\t", 'pop', 'diffs/comps', 'diffs/comps*1000', 'n_regions', 'total_count_diffs', 'total_count_comparisons' ), "\n";
  for my $pop (sort keys %data)
    {
      my $n = $data{$pop}{n};
      my $cd = $data{$pop}{countdiffs};
      my $cc = $data{$pop}{countcomps};
      my $pr = $cd / $cc;
      print $OUTF join( "\t", $pop, $pr, $pr*1000, $n, $cd, $cc ), "\n";
    }
  print $OUTF "\n";
  stdclose ( $OUTF );
} # sub saveoutput



###############################################################
sub commify
###############################################################
# http://perldoc.perl.org/perlfaq5.html#How-can-I-output-my-numbers-with-commas
  {
    local $_ = shift;
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
  } # commify



###############################################################
sub tlog { my ( $prefix, $suffix, $noreturn ) = @_;
###############################################################
  $prefix //= '';
  $suffix //= '';
  @_ = localtime(time);
  my $t = $prefix . '  '
    . sprintf("%04d/%02d/%02d %02d:%02d:%02d", $_[5]+1900, $_[4]+1, $_[3], @_[2,1,0])
    . '  ' . $suffix;
  $t =~ s/^ +//;
  $t =~ s/ +$//;
  print $t;
  unless ( $noreturn ) { print "\n"; }
} # sub tlog



###############################################################
sub run { my ( $command, $errorokay ) = @_;
###############################################################
# Run a system command, but use bash. Perl by default uses sh.
# Command may be piped multiple commands or redirected from or to files.
  my $result = system( 'bash', '-o', 'pipefail', '-c', $command );
  if ( ( $result ) and ( ! $errorokay ) )
    {
      # exitvalue 141 can occur if piping to head or grep -q
      my $exitvalue = $result >> 8;
      my $signal = $result & 255;
      die( "Error $exitvalue:$signal running command \"$command\"\n" );
    }
  return( $result );
} # sub run



###############################################################
sub stdopen { my ( $mode, $filename, $extratext ) = @_;
###############################################################
# a replacement for the three-parameter open which also allows
# the use of "-" as the file name to mean STDIN or STDOUT
  my $fh;  # the file handle
  if ( $filename eq "-" )  # only exact match to "-" has special meaning
    {
      if ( $mode =~ m/>/ )
        { $fh = *STDOUT }
      else
        { $fh = *STDIN }
    }
  else
    {
      # supplemental passed text for error messages, need one more space
      if ( defined $extratext )
        { $extratext .= " " }
      else
        { $extratext = "" }

      my $text;  # this is only used for error message
      if ( $mode =~ m/^\+?>>/ )  # ">>" or "+>>"
        { $text = "append" }
      elsif ( $mode =~ m/^\+?>/ )  # ">" or "+>"
        { $text = "output" }
      elsif ( $mode =~ m/^\+?</ )  # "<" or "+<"
        { $text = "input" }
      elsif ( $mode eq "-|" )
        { $text = "piped input" }
      elsif ( $mode eq "|-" )
        { $text = "piped output" }
      else
        { die "Error, unsupported file mode \"$mode\" specified to stdopen( $mode, $filename, $extratext )\n"; }

      # if file name ends in ".gz", gzip compression is assumed, and handle it transparently
      if ( $filename =~ m/\.gz$/ )
        {
          if ( $mode =~ m/^>$/ ) # output mode
            { $mode = "|-"; $filename = "gzip -c > \"$filename\""; }
          elsif ( $mode =~ m/^<$/ ) # input mode
            { $mode = "-|"; $filename = "gunzip -c \"$filename\""; }
          elsif ( $mode =~ m/^>>$/ ) # append mode
            { $mode = "|-"; $filename = "gzip -c >> \"$filename\""; }
          else
            { die "Error, can't handle gzip compression with mode \"$mode\" for file \"filename\"\n"; }
        } # if gzip compressed file
      elsif ( $filename	=~ m/\.bz2$/ )
       	{
          if ( $mode =~ m/^>$/ ) # output mode
            { $mode = "|-"; $filename = "bzip2 -c > \"$filename\""; }
          elsif ( $mode =~ m/^<$/ ) # input mode
            { $mode = "-|"; $filename = "bunzip2 -c \"$filename\""; }
          else
            { die "Error, can't handle bzip2 compression with mode \"$mode\" for file \"filename\"\n"; }
       	}
      open ( $fh, $mode, $filename ) or die ( "Error opening ${extratext}file \"$filename\" for $text: $!\n" );
    }
  # return the opened file handle to the caller
  return $fh;
} # sub stdopen



###############################################################
sub stdclose { my ( $fh ) = @_;
###############################################################
# same as built-in close, except in case of STDIN or STDOUT,
# and in this case the file handle is not closed

  unless ( fileno($fh) <= 2 )  # if file number is this low, is stdin or stdout or stderr
    { close ( $fh ) or die ( "Error closing file handle: $!\n" ); }

} # sub stdclose



#eof
