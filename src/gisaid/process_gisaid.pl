#! /usr/bin/perl -w

use strict;
use Getopt::Std;

my $MAX_PER_MONTH = 50000;
my $MIN_LENGTH = 25000;

my $USAGE = "$0 [options] input_metadata.tar.xz input_sequences.tar.xz output_metadata.tsv output_sequences.fasta.gz

Options:
-M <integer>   Maximum number of sequences per month, default $MAX_PER_MONTH
-L <integer>   Minimum length of a sequence, excluding unknwown bases, default $MIN_LENGTH


The script performs several preprocessing steps on files downloaded from GISAID
- records with incomplete date or incomplete sequence (see -L option) are skipped
- duplicate records with the same virus name are also skipped
- spaces in the virus name are converted to underscores, as tools working with fasta format do not handle sequence names with spaces
- from each month, a random sample of record is selected, to achive approximately the size given by -M option. 
  If there are fewer record in that month than a given threshold, all records are selected

Input metadata should be a .tar.xz file containing metadata.tsv file, 
input sequences should be a .tar.xz file containing sequence.fasta file.

Selected records are written to a tsv file, which is not compressed,
selected sequences are written to a gzipped fasta file.

Tools tar and gz are used by te script and should be available on the system.
";

my %Options;
getopts('M:L:', \%Options);
if(exists $Options{'M'}) {
    $MAX_PER_MONTH = $Options{'M'};
}
if(exists $Options{'L'}) {
    $MIN_LENGTH = $Options{'L'};
}

die $USAGE unless @ARGV==4;

my ($meta_in, $fasta_in, $meta_out, $fasta_out) = @ARGV;

die "$meta_in does not exists" unless -r $meta_in;
die "$fasta_in does not exists" unless -r $fasta_in;


my $counts = {};         # dictionary with count of smaples for each month
my $selected = {};       # dictionary with a list of selected samples
my $duplicates = {};     # dictionary with a list of duplicates

# in phase 1, read the file, find duplicates and compute counts
read_metadata($meta_in, $meta_out, $counts, $selected, $duplicates, 1);
# in phase 2, sample according to counts, skip duplicates, print output metadata and add names to $selected
read_metadata($meta_in, $meta_out, $counts, $selected, $duplicates, 2);

# read fasta, print selected record to output
read_fasta($fasta_in, $fasta_out, $selected);

##############
sub read_fasta {
    my ($fasta_in, $fasta_out, $selected) = @_;

    # open tar file for reading (reading file sequences.fasta from the tar file)
    my $in;
    open $in, "tar -xOf $fasta_in sequences.fasta | "  or die "Cannot open $fasta_in";

    # open file for writing, output is sent to gzip
    my $out;
    open $out, "| gzip -9 -c > $fasta_out"; 

    # read input line by line
    my $print = 0;   # we are currently not printing input lines
    while(my $line = <$in>) {

	if($line=~/>([^\|]+)/) {   # line with name (everything after | is ignored)
	    my $name = $1;
	    $name =~ s/\s+$//;    # remove whitespace at the end, if any
	    $name =~ s/ /_/g;     # replace other spaces by underscores
	    if(exists $selected->{$name}) {		# check if sequence was selected, if yes, start printing
		$print = 1;
		$line = ">$name\n";
	    } else {   # sequence was not selected, set printing to 0
		$print = 0;
	    }
	}
	
	if($print) {   # print the line if we are currently in printing mode
	    print $out $line;
	}
    }
    # close both files
    close $in;
    close $out or die;
}

##############
sub read_metadata {
    my ($meta_in, $meta_out, $counts, $selected, $duplicates, $phase) = @_;

    # open tar file for reading (reading file metadata.tsv from the tar file)
    my $in;
    open $in, "tar -xOf $meta_in metadata.tsv | "  or die "Cannot open $meta_in";

    # open tsv file for writing
    my $out;
    if($phase == 2) {
	open $out, ">", $meta_out or die "Cannot open $meta_out";
    }
    
    my $line_num = 0;
    my $ids = {};  # list of all encountered IDs, for each store its month
    while(my $line = <$in>) {
	$line_num++;
	chomp $line;  # rmeove end of line
	my @parts = split "\t", $line;  # split to columns
	next if @parts==0;  # skip empty lines
	
        # >>> added by Askar
        next if @parts < 24;
        die "bad line $line (expected 24 columns)" unless @parts == 24;  # check that we have 22 columns
        # <<< added by Askar
	
        my ($name, $date, $ncount, $seq_len) = @parts[0, 5, 22, 8];  # store needed columns to variables
	
	if($line_num == 1) {    
	    # check header line, if selected columns have the correct names
            
            my $separated_line = $line;
            $separated_line =~ s/\t/|/g;
            print $separated_line;
            
	    die "bad header line $line" unless
		$name eq "Virus name"
		and $date eq "Collection date"
		and $ncount eq "N-Content"
		and $seq_len eq "Sequence length";

	    if($phase == 2) {  # print header in phase 2
		print $out $line, "\n";
	    }
	    next;  # go to the next line
	} 

	# length of good sequence is the total length minus Ns
	if($ncount eq "") {
	    $ncount = 0;
	}
	my $good_seq = $seq_len * (1 - $ncount);	


	$name =~ s/ /_/g;  # replaces spaces in names by underscores
	my $month = substr($date, 0, 7);  # get month from the date

	my $is_good;
	# mark incomplete dates and short sequences as "bad"
	if(length($date) < 10 || $good_seq < $MIN_LENGTH) {
	    $is_good = 0;
	    $month = "xx";  # fake month for bad sequences
	} else {
	    $is_good = 1;
	}
		    
	if($phase==1) {
	    # in phase 1, check duplicates and count records for each month
	    if(! exists $ids->{$name} && ! exists $duplicates->{$name}) {    # new id - count in stats
		$ids->{$name} = $month;
		$counts->{$month}++;
	    } elsif(! exists $duplicates->{$name}) {
		# duplicated id - mark as such when seen the first time
		$duplicates->{$name} = 1;
		my $prev_month = $ids->{$name};
		$counts->{$prev_month}--;  # remove from counts
	    }
	} else {
	    # phase 2

	    # skip duplicates and bad records
	    next if exists $duplicates->{$name} || $is_good == 0;  
    
	    my $n = $counts->{$month};  # real size of the month group
	    if(rand() < $MAX_PER_MONTH / $n) {   # subsample based on group size
		$parts[0] = $name;                    # rewrite virus name with our version (no spaces)
		print $out join("\t", @parts), "\n";  # write current line to output, with our virus name
		$selected->{$name} = 1;               # add to the dictionary of selected records
	    }
	}
    }
    
    $ids = undef;   # free list of all ids
    close $in;
    if($phase == 2) {
	close $out or die;
	$duplicates = undef;   # free the list of duplicates
    }
}
    
  
