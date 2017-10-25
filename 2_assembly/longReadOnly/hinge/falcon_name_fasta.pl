#!/usr/bin/env perl
#
## falcon_name_fasta.pl -- Erich Schwarz <ems394@cornell.edu>, 11/13/2015.
## Purpose: rename FASTA reads with format acceptable to FALCON; as side product, generate old-to-new name table.
#
use strict;
use warnings;
use Getopt::Long;

my $i            = 0;
my $prefix       = q{};
my $table        = q{};
my $infile       = q{};
my $outfile      = q{};
my $header       = q{};
my $orig_seqname = q{};
my @orig_seqs    = ();

my $data_ref;
my $help;

GetOptions ( 'infile=s'  => \$infile,
             'outfile:s' => \$outfile,
                          'prefix:s'  => \$prefix,
                                       'table:s'   => \$table,
                                                    'help'      => \$help,  );

                                                    if ( $prefix =~ /\W/xms ) {
                                                        die "Cannot accept prefix with space or non-word characters (such as \".\"): $prefix\n";
                                                        }

                                                        if ( $outfile =~ /[^\w\.]/xms ) {
                                                            die "Cannot accept outfile with space or non-standard characters: $outfile\n";
                                                            }

                                                            # Have human-readable default value:
                                                            $prefix ||= 'falcon_read';
                                                            $prefix = safename($prefix);

                                                            $outfile ||= "$infile.outfile";
                                                            $outfile = safename($outfile);

                                                            $table ||= 'old2new_names.txt';
                                                            $table = safename($table);

                                                            if ( $help or (! $infile ) ) { 
                                                                die "Format: falcon_name_fasta.pl\n",
                                                                        "    --infile|-i    <input stream/files>\n",
                                                                                "    --outfile|-o   [output file name; default is (infile).outfile]\n",
                                                                                        "    --prefix|-p    [optional prefix; default is \"falcon_read\"]\n",
                                                                                                "    --table|-t     [table of old to new read names; default is \"old2new_names.txt\"]\n",
                                                                                                        "    --help|-h      [print this message]\n",
                                                                                                        }

                                                                                                       # Accept either a stream from '-' or a standard file.
                                                                                                        my $INFILE;
                                                                                                        if ($infile eq '-') {
                                                                                                            # Special case: get the stdin handle
                                                                                                                $INFILE = *STDIN{IO};
                                                                                                                }
                                                                                                                else {
                                                                                                                    # Standard case: open the file
                                                                                                                        open $INFILE, '<', $infile or die "Can't open input file $infile. $!\n";
                                                                                                                        }

                                                                                                                        open my $OUTFILE, '>', $outfile;
                                                                                                                        open my $TABLE, '>', $table;

                                                                                                                        while (my $input = <$INFILE>) { 
                                                                                                                            chomp $input;
                                                                                                                                if ( $input =~ /\A > (\S+) /xms ) { 
                                                                                                                                        $orig_seqname = $1;
                                                                                                                                                if ( exists $data_ref->{'orig_seqname'}->{$orig_seqname} ) { 
                                                                                                                                                            die "Redundant sequence name: $orig_seqname\n";
                                                                                                                                                                    }
                                                                                                                                                                            $data_ref->{'orig_seqname'}->{$orig_seqname}->{'seen'} = 1;
                                                                                                                                                                                    push @orig_seqs, $orig_seqname;
                                                                                                                                                                                        }
                                                                                                                                                                                            elsif ( $input =~ / \A \s* [A-Za-z] /xms ) { 
                                                                                                                                                                                                    $input =~ s/\s//g;
                                                                                                                                                                                                            if ( $input =~ / [^ACGTNacgtn] /xms ) { 
                                                                                                                                                                                                                        die "Can't parse: $input\n";
                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                        $data_ref->{'orig_seqname'}->{$orig_seqname}->{'sequence'} .= $input;
                                                                                                                                                                                                                                            }
                                                                                                                                                                                                                                                else { 
                                                                                                                                                                                                                                                        if ( $input !~ /\A \s* \z/xms ) { 
                                                                                                                                                                                                                                                                    die "Can't parse: $input\n";
                                                                                                                                                                                                                                                                            }
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                }
                                                                                                                                                                                                                                                                                close $INFILE;

                                                                                                                                                                                                                                                                                foreach my $orig_seq1 ( @orig_seqs ) {
                                                                                                                                                                                                                                                                                    $data_ref->{'orig_seqname'}->{$orig_seq1}->{'length'} 
                                                                                                                                                                                                                                                                                            = length( $data_ref->{'orig_seqname'}->{$orig_seq1}->{'sequence'} );
                                                                                                                                                                                                                                                                                            }

                                                                                                                                                                                                                                                                                            my $seq_count = @orig_seqs;
                                                                                                                                                                                                                                                                                            my $DIGITS = length($seq_count);

                                                                                                                                                                                                                                                                                            foreach my $orig_seq2 (@orig_seqs) { 
                                                                                                                                                                                                                                                                                                $i++;
                                                                                                                                                                                                                                                                                                    my $serial_no = $i;
                                                                                                                                                                                                                                                                                                        my $sf_format = '%0' . $DIGITS . 'u';
                                                                                                                                                                                                                                                                                                            $serial_no = sprintf($sf_format, $serial_no) or die "Can't zero-pad serial number $serial_no\n";

                                                                                                                                                                                                                                                                                                                my $new_name = $prefix . q{/} . $serial_no . q{/0_} . $data_ref->{'orig_seqname'}->{$orig_seq2}->{'length'};

                                                                                                                                                                                                                                                                                                                    print $OUTFILE '>' . "$new_name\n";
                                                                                                                                                                                                                                                                                                                        print $TABLE "$orig_seq2\t$new_name\n";
                                                                                                                                                                                                                                                                                                                            my @output_lines 
                                                                                                                                                                                                                                                                                                                                    = unpack("a60" 
                                                                                                                                                                                                                                                                                                                                                     x ($data_ref->{'orig_seqname'}->{$orig_seq2}->{'length'}/60 + 1), 
                                                                                                                                                                                                                                                                                                                                                                      $data_ref->{'orig_seqname'}->{$orig_seq2}->{'sequence'});
                                                                                                                                                                                                                                                                                                                                                                          foreach my $output_line (@output_lines) { 
                                                                                                                                                                                                                                                                                                                                                                                  if ($output_line =~ /\S/) { 
                                                                                                                                                                                                                                                                                                                                                                                              print $OUTFILE "$output_line\n";
                                                                                                                                                                                                                                                                                                                                                                                                      }
                                                                                                                                                                                                                                                                                                                                                                                                          }
                                                                                                                                                                                                                                                                                                                                                                                                          }

                                                                                                                                                                                                                                                                                                                                                                                                          close $OUTFILE;
                                                                                                                                                                                                                                                                                                                                                                                                          close $TABLE;

                                                                                                                                                                                                                                                                                                                                                                                                          sub safename {
                                                                                                                                                                                                                                                                                                                                                                                                              my $_filename = $_[0];
                                                                                                                                                                                                                                                                                                                                                                                                                  my $_orig_filename = $_filename;
                                                                                                                                                                                                                                                                                                                                                                                                                      if (-e $_orig_filename) {
                                                                                                                                                                                                                                                                                                                                                                                                                              my $_suffix1 = 1;
                                                                                                                                                                                                                                                                                                                                                                                                                                      $_filename = $_filename . ".$_suffix1";
                                                                                                                                                                                                                                                                                                                                                                                                                                              while (-e $_filename) {
                                                                                                                                                                                                                                                                                                                                                                                                                                                          $_suffix1++;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      $_filename =~ s/\.\d+\z//xms;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  $_filename = $_filename . ".$_suffix1";
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  return $_filename;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  }


