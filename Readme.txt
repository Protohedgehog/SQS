How to execute SQS in both R and Perl.

First step is to download a list of fossil occurrences from the Paleobiology Database: https://paleobiodb.org/cgi-bin/bridge.pl?a=displayDownloadGenerator

Handily, this generates a .csv file, which makes it dead easy to use in R and Perl.

PERL
Make sure Perl is installed. Have the time bins.txt, Perl script, and occurrences file in the same working directory. Make sure the Perl script is pointing to each of these. There are many options you can play with here, and it's worth consulting Alroy's papers so you know what these do in advance.
Simply double click the script then to run it. It should cycle through and spit out several files. Then simply go through the results and voila!

R
I think you have to do this one time bin by time bin, which can get a bit tedious. Simply read in a list of taxonomic names, and convert this to an occurrence list using the table() function.
You can play with singletons and dominant taxa, but these shouldn't have a major impact on your resulting curves.
The bootstrap code will generate the median and 95% confidence intervals, which are essential when looking at the shape of diversity.
