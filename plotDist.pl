#!/usr/bin/perl
use strict;
use warnings;

# Check that a file name has been provided
if (@ARGV < 1) {
    die "Usage: $0 <csv_file>\n";
}

my $csv_file = $ARGV[0];

# Open the CSV file
open my $fh, '<', $csv_file or die "Could not open file '$csv_file': $!";

# Prepare an array of [ED, TSS] pairs
my @data_points;

my $ed_idx = -1;
my $tss_idx = -1;
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split /,/, $line;
    if ( $line =~ /^s/ ) {
      for ( my $i = 0; $i < @fields; $i++ ) {
        if ( $fields[$i] eq "ED" ) { 
           $ed_idx = $i;
        }elsif ( $fields[$i] eq "TSS" ) {
           $tss_idx = $i;
        }
      }
      #print STDERR "ED index: $ed_idx, TSS index: $tss_idx\n";
      next;
    }

    # fields: s1, s2, ED, TSS
    my ($ed, $tss) = ($fields[$ed_idx], $fields[$tss_idx]);

    # Make sure ED and TSS are numeric
    next unless defined $ed and defined $tss;
    next if $ed  !~ /^[0-9.+-eE]+$/;
    next if $tss !~ /^[0-9.+-eE]+$/;

    push @data_points, [$ed, $tss];
}
close $fh;

#
# Now, we print out the HTML for Google Charts.
#

# Print the HTML header
print <<"END_HTML";
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>TSS vs ED</title>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['corechart']});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {
        // Create a data table
        var data = new google.visualization.DataTable();
        data.addColumn('number', 'ED');
        data.addColumn('number', 'TSS');

        data.addRows([
END_HTML

# Print the data rows
foreach my $pair (@data_points) {
    my ($ed, $tss) = @$pair;
    print "          [$ed, $tss],\n";
}

# Print the final part of the HTML / JS
print <<"END_HTML";
        ]);

        var options = {
          title: 'TSS vs ED',
          hAxis: {title: 'ED'},
          vAxis: {title: 'TSS'},
          width:  900,
          height: 500
        };

        var chart = new google.visualization.ScatterChart(document.getElementById('chart_div'));
        chart.draw(data, options);
      }
    </script>
</head>
<body>
  <div id="chart_div" style="width: 900px; height: 500px;"></div>
</body>
</html>
END_HTML

