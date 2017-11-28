#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Spec;
use File::stat;

my ($run_summary, $sample, $sequence_type, $help, $debug, %run, %lane, %read_type, $command, $result);


#-- prints usage if no command line parameters are passed or there is an unknown.
usage() if ( @ARGV < 4 or
          ! GetOptions('h|?' => \$help, 'r:s' => \$run_summary, 's:s' => \$sample, 't:s' => \$sequence_type)
           );


#-- checks if the Run_Summary file exists.
$run_summary = File::Spec->rel2abs($run_summary);
die "Error : $run_summary not found.\n" unless ( -e $run_summary);

#-- check if genoem or exome
die "Error : invalid sequece type. Use -t [genome|exome].\n" unless($sequence_type =~ /(exome|genome)/i);

#-- pre-defined variables
$debug = 0;
%run = ("1" => "Run1", "2" => "Run2", "3" => "Run3", "4" => "Run4", "5" => "Run5", "6" => "Run6", "7" => "Run7");
%lane = (
	"L001" => 1, "L002" => 2, "L003" => 3, "L004" => 4, 
	"L005" => 5, "L006" => 6, "L007" => 7, "L008" => 8
	); 
%read_type = ("R1" => 1, "R2" => 2);


#-- check if sample name exist in Run-Summary.
$command = "grep -w '$sample' $run_summary";
$result = `$command`;
die "Error : failed to grep $sample at $run_summary\n" unless ($? == 0);
print ">The sample name, $sample,  has been checked ...\n" if($debug == 1);


my $num_error = 0;
my @runs = split(/\n/,$result);
foreach my $line (@runs)
{
  $line =~ s/\r//g;
  print "$line\n" if($debug == 1);
  my @cols  = split(/\t/,$line);
  my $myRun = $cols[0]; #$1;
  my $myLane = $cols[1]; #$2;
  my $myPath = $cols[2]; #$3;
  my $myType = $cols[3]; #$4;
  my $mySampleNickname = $cols[9];
  my $myTopstrand = $cols[10]; #$5;
  my $myExomeKit = $cols[11];
  my $myDir = $myPath."/".$myRun;


  #-- check if run-directory exist.
  $myDir = File::Spec->rel2abs($myDir);
  die "Error : the run-directory,$myDir, not found.\n" unless ( -e $myDir);   
  print ">Ther run directory, $myDir, has been checked ...\n" if($debug == 1);


  #-- for each run, create a hash like {lane => {read_type => [fastq1, fastq2, ... and so forth]}}.
  my $hash_ref = {};
  $command = "ls $myDir/*.fastq* | xargs -n1 basename";
  print "command : $command\n" if($debug == 1);
  $result = `$command`;
  chomp($result);
  my @fastq_names = split(/\n/,$result);
  foreach my $fastq_name (@fastq_names)
  {
	print "fastq_name : $fastq_name\n" if($debug ==1);

        my @tmp_arr = split(/[_.]/,$fastq_name);
  	print "$tmp_arr[0],$tmp_arr[1],$tmp_arr[2],$tmp_arr[3]\n" if($debug == 1);
        my ($tmp_lane, $tmp_read)  = ( $lane{$tmp_arr[2]}, $read_type{$tmp_arr[3]} );
        $tmp_arr[4] =~ s/^0+//;
        my $tmp_num = $tmp_arr[4];
        ${$hash_ref}{$tmp_lane}{$tmp_read}[$tmp_num-1] = $fastq_name;
        print "$tmp_lane, $tmp_read, $tmp_num\n" if($debug == 1);
  }


  #-- check if Run_Summary information is valid
  $num_error++ unless( isLaneValid($hash_ref, $myDir, $myLane, $myType) );
  $num_error++ unless( isLaneNumberValid($hash_ref, $myDir, $myLane, $myType) );
  $num_error++ unless( isReadTypeValid($hash_ref, $myDir, $myLane, $myType) );
  if($myType eq 'PE')
  {
  	$num_error++ unless( isEqualArraySize($hash_ref, $myDir, $myLane, $myType) );
  }
  unless($myTopstrand =~ /^na$/i || $myTopstrand eq '')
  {
  	$num_error++ unless( isTopstrandFileValid($sample,$mySampleNickname,$myTopstrand) ); 
  }

  die "Error : invalid the exome_kit value. Please double check.\n" if($sequence_type =~ /exome/i && $myExomeKit !~ /(37|50|65|nimble|[Rr]oche|VCRome2_1|AgilentV4UTR|AgilentV4|AgilentCRE|AgilentV5|[Aa]lopecia|IDTERPv1)/);
}
die "Invalid information between Run-Summary1.8 and sequence directory!\n" if ($num_error > 0);


sub isTopstrandFileValid
{
  my ($tmpSample, $tmpSampleID, $tmpTopstrandFile) = @_;

  $command = "ssh 10.73.50.41 'ls $tmpTopstrandFile'";
  $result = `$command`;

  unless($? == 0)
  {
  	print "Error : the topstand file, $tmpTopstrandFile,  is not found.\n";
  	print ERROR "Error : the topstand file, $tmpTopstrandFile, is not found.\n";
  	return 0;
  }  	
  else
  {
  	$tmpTopstrandFile =~ /(BUILD37)/i;
  	print "$tmpSample\'s tostrand version : $1\n"; 
  }

  #$command = "ssh 10.73.50.41 'grep -m 1 $tmpSampleID $tmpTopstrandFile | wc -l'";
  $command = "ssh 10.73.50.41 \"awk -F\",\" \'\{ print \$2 \}\' $tmpTopstrandFile | grep -w -m 1 $tmpSampleID | wc -l\"";
  $result = `$command`;
  unless($result == 1){
	print "Error : the genotype information of $tmpSample do not exist in $tmpTopstrandFile\n";
  	print ERROR "Error : the genotype information of $tmpSample do not exist in $tmpTopstrandFile\n";
	return 0;
  }
  return 1;
}


sub isLaneValid
{
  my ($tmpHash, $tmpDir, $tmpLane, $tmpType) = @_;
  my $not_found = 0;

  $tmpLane =~ s/ //g;
  if($tmpLane =~ /(\d)-(\d)/)		#-- consecutive multi lanes with delimited dash
  {
  	foreach($1 .. $2)
	{
		unless(exists $tmpHash->{$_})
		{
			print "Error : The lane $_ is not found in $tmpDir\n";
			print ERROR "Error : The lane $_ is not found in $tmpDir\n";
			$not_found++;
		}	
		
	}	
	return 0 if($not_found > 0);
  }
  elsif($tmpLane =~ /,/)		#-- multi lanes with delimited comma
  {
	foreach(split(/,/,$tmpLane))
	{
                unless(exists $tmpHash->{$_})
                {
                        print "Error : The lane $_ is not found in $tmpDir\n";
                        print ERROR "Error : The lane $_ is not found in $tmpDir\n";
                        $not_found++;
                }
		
	}
	return 0 if($not_found > 0);
  }
  else					#-- single lane
  {
  	foreach(split(//,$tmpLane))
	{
  		unless(exists $tmpHash->{$_})
		{
        		print "Error : The lane $_ is not found in $tmpDir\n";
        		print ERROR "Error : The lane $_ is not found in $tmpDir\n";
			$not_found++;
		}
  	}
	return 0 if($not_found > 0);
  }

  return 1;
}


sub isLaneNumberValid
{
  my ($tmpHash, $tmpDir, $tmpLane, $tmpType) = @_;
  my $laneNum = keys %{$tmpHash};
  my $lane_count = 0;

  $tmpLane =~ s/ //g;
  if($tmpLane =~ /(\d)-(\d)/)           #-- consecutive multi lanes with delimited dash
  {
 
        foreach($1 .. $2)
        {
        	$lane_count++;
        }
        unless($laneNum == $lane_count)
  	{
  		print "Error : the number of lanes is different at $tmpLane in $tmpDir\n";
  		print ERROR "Error : the number of lanes is different at $tmpLane in $tmpDir\n";
		return 0;
	}
  }
  elsif($tmpLane =~ /,/)                #-- multi lanes with delimited comma
  {
        foreach(split(/,/,$tmpLane))
        {
  		$lane_count++;
        }
        unless($laneNum == $lane_count)
        {
                print "Error : the number of lanes is different at $tmpLane in $tmpDir\n";
                print ERROR "Error : the number of lanes is different at $tmpLane in $tmpDir\n";
                return 0;
        }

  }
  else                                  #-- single lane
  {
	foreach(split(//,$tmpLane))
        {
                $lane_count++;
        }
        unless($laneNum == $lane_count)
        {
                print "Error : the number of lanes is different at $tmpLane in $tmpDir\n";
                print ERROR "Error : the number of lanes is different at $tmpLane in $tmpDir\n";
                return 0;
        }
  }

  return 1;
}



sub isReadTypeValid
{
  my ($tmpHash, $tmpDir, $tmpLane, $tmpType) = @_;
  my $invalid = 0;

  if($tmpType eq 'PE')
  {
	#my $laneNum = keys %{$tmpHash};
	#print "laneNum : $laneNum\n";
  	foreach my $lane (keys %{$tmpHash})
	{
		my $typeNum = keys %{$tmpHash->{$lane}};
		unless($typeNum == 2)
		{
			print "Error : the read-type is not paied reads (PE) at $tmpDir.\n";
			print ERROR "Error : the read-type is not paied reads (PE) at $tmpDir.\n";
			$invalid++;
		}
	}
  	return 0 if($invalid > 0);
  }
  elsif($tmpType eq 'SE')
  {
        foreach my $lane (keys %{$tmpHash})
        {
                my $typeNum = keys %{$tmpHash->{$lane}};
                unless($typeNum == 1)
                {
                        print "Error : the read-type is not single read (SE) $tmpDir.\n";
                        print ERROR "Error : the read-type is not single read (SE) $tmpDir.\n";
                        $invalid++;
                }
        }
  	return 0 if($invalid > 0);

  }
  else
  {
  	print "Error : cannot determin the read type (PE, SE) at $tmpDir.\n";
	print ERROR "Error : cannot determin the read type (PE, SE) at $tmpDir.\n";
	return 0;
  }
  
  return 1;
}


sub isEqualArraySize
{
  my ($tmpHash, $tmpDir, $tmpLane, $tmpType) = @_;
  my $not_equal = 0;

  foreach my $lane (keys %{$tmpHash})
  {
  	my $arr1_ref_num = @{${$tmpHash->{$lane}}{"1"}};
  	my $arr2_ref_num = @{${$tmpHash->{$lane}}{"2"}};
  	#print "numbers : $arr1_ref_num, $arr2_ref_num\n";
  	unless($arr1_ref_num == $arr2_ref_num)
  	{
  		print "Error : the number of split-reads is different at $tmpDir.\n";
  		print ERROR "Error : the number of split-reads is different at $tmpDir.\n";
  		$not_equal++;
  	}
  }
  ($not_equal == 0) ? (return 1) : (return 0);
}


sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "#-- fastq_check_cassava18.pl :\n";
  print "#   is written by Hee Shin Kim.\n";
  print "#   reads Run_Summary_v1.8.txt for a sample and checks validation such as sequence_location, lanes, read_type, and read_number.\n";
  print "#   returns error message if Run_Summary_v1.8 is invalid.\n";
  print "#\n";
  print "#   usage for single sample: fastq_check_cassava18.pl -r Run_Summary_v1.8.txt -s your_sample_name -t [exome|genome][-h]\n";
  print "#\n";
  print "# ex) perl fastq_check_cassava18.pl -r /nfs/chgv/seqpipe01/WGseq/Run_Summary_v1.8.txt -s your_sample -t genome\n";
  print "#\n";
  print "#   usage for mutilple samples : for f in `your_sample_lis.txt`\; do perl fastq_check_cassava18.pl -r /Run_Summary_v1.8.txt -s \$f -t genome\; done\n";
  print "#\n";
  print "# ex)  for f in `your_sample_lis.txt`\; do perl fastq_check_cassava18.pl -r /nfs/chgv/seqpipe01/WGseq/Run_Summary_v1.8.txt -s \$f -t genome\; done\n";
  exit;
}

