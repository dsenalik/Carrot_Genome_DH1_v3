#!/usr/bin/perl
use strict;
use warnings;

# report on gene changes

#if ( $#ARGV < 3 )
#  { die "Syntax: gff3 oldfasta newfasta report\n"; }

# configuration variables
my ( $pfx ) = $0 =~ m/\/(\d+)[^\/]*$/;



# input files
my $gfffilename = $ARGV[0] // 'DCv3.all.20220211.gff3.gz';
my $oldfasta = $ARGV[1] // '/vcru_share_s/carrot/dcarv3d/ann/DCv3.pep.fa.gz';
my $newfasta = $ARGV[2] // '20220211.unwrap.faa';



# output files
my $reportfilename = $ARGV[3] // "notebook/$pfx.report.20220211.tsv";



# global variables
my %ispseudo;
my %oldseq;
my %newseq;



loadgff3( $gfffilename );
loadfasta( $oldfasta, \%oldseq );
loadfasta( $newfasta, \%newseq );
makereport( $reportfilename, \%oldseq, \%newseq );



system( "cp -puv --no-preserve=owner $0 notebook/" );
tlog( "$0 Done" );
exit 0;



###############################################################
sub loadgff3 { my ( $infilename ) = @_;
###############################################################
  tlog( "Loading \"$infilename\"" );
  my $nlines = 0;
  my $INF = stdopen ( '<', $infilename );
  while ( my $aline = <$INF> )
    {
      $nlines++;
      $aline =~ s/[\r\n]//g;
      my @cols = split ( /\t/, $aline );
      if ( ( $cols[2] ) and ( $cols[2] eq 'gene' ) and ( $cols[8] =~ m|pseudo=| ) )
        {
          unless ( $cols[8] =~ m|ID=([^;]+)| )
            { die "No ID= line $nlines \"$aline\"\n"; }
          my $id = $1;
          $ispseudo{$id} = 1;
        }
    } # while <$INF>
  stdclose ( $INF );
  tlog( commify($nlines) . " lines processed" );
} # sub loadgff3



###############################################################
sub loadfasta { my ( $infilename, $dataref ) = @_;
###############################################################
  tlog( "Loading \"$infilename\"" );
  my $nlines = 0;
  my $INF = stdopen ( '<', $infilename );
  my $id;
  while ( my $aline = <$INF> )
    {
      $nlines++;
      $aline =~ s/[\r\n]//g;
      if ( $aline =~ m|>(\S+)| )
        {
          $id = $1;
          $id =~ s|\.mRNA.*$||;
        }
      else
        { $dataref->{$id} .= $aline; }
    } # while <$INF>
  stdclose ( $INF );
  tlog( commify($nlines) . " lines processed" );
} # sub loadfasta



###############################################################
sub makereport { my ( $reportfilename, $olddataref, $newdataref ) = @_;
###############################################################
  tlog( "Generating report to \"$reportfilename\"" );
  my %counts;
  my %seen;
  my $REPT = stdopen ( '>', $reportfilename );
  print $REPT join( "\t", 'geneid', 'status' ), "\n";
  for my $id ( sort keys %$olddataref )
    {
      $seen{$id} = 1;
      my $oldseq = $olddataref->{$id};
      my $newseq = $newdataref->{$id};
      my $status = 'same';
      if ( $ispseudo{$id} )
        { $status = 'pseudo'; }
      elsif ( !$newseq )
        { $status = 'deleted'; }
      elsif ( $oldseq ne $newseq )
        { $status = 'changed'; }
      $counts{$status}++;
      print $REPT join ( "\t", $id, $status ), "\n";
    } # while <$INF>

  # genes in new not in old
  for my $id ( sort keys %$newdataref )
    {
      unless ( ( $seen{$id} ) or ( $id =~ m|T| ) )  # excludes organellar
        {
          my $status = 'new';
          $counts{$status}++;
          print $REPT join ( "\t", $id, $status ), "\n";
        }
    }

  stdclose ( $REPT );
  foreach my $key ( sort keys %counts )
    {
      print commify($counts{$key}), "\t", $key, "\n";
    }
} # sub makereport



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
