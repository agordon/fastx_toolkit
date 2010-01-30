#!/usr/bin/env bash

#    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
#    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    This script was modified by Oliver Tam (tam@cshl.edu)

usage()
{
	echo "FASTA/Q Nucleotide Distribution Plotter"
	echo
	echo "Usage: $0 [-i INPUT.TXT] [-t TITLE] [-p] [-o OUTPUT]"
	echo
	echo "  [-p]           - Generate PostScript (.PS) file. Default is PNG image."
	echo "  [-i INPUT.TXT] - Input file. Should be the output of \"fastx_quality_statistics\" program."
	echo "  [-o OUTPUT]    - Output file name. default is STDOUT."
	echo "  [-t TITLE]     - Title - will be plotted on the graph."
	echo
	exit 
}

#
# Input Data columns: #pos	cnt	min	max	sum       	mean	Q1	med	Q3	IQR	lW	rW A_Count	C_Count	G_Count	T_Count	N_Count
#  As produced by "fastq_quality_statistics" program

TITLE=""					# default title is empty
FILENAME=""
OUTPUTTERM="set term png size 1048,768"		# default output terminal is "PDF"
OUTPUTFILE="/dev/stdout"   			# Default output file is simply "stdout"
while getopts ":t:i:o:ph" Option
	do
	case $Option in
		t ) TITLE="for $OPTARG" ;;
		i ) FILENAME=$OPTARG ;;
		o ) OUTPUTFILE="$OPTARG" ;;
		p ) OUTPUTTERM="set term postscript enhanced color \"Helvetica\" 8" ;;
		h ) usage ;;
		* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
	esac
done
shift $(($OPTIND - 1)) 


if [ -z "$FILENAME" ]; then
	usage
fi

if [ ! -r "$FILENAME" ]; then
	echo "Error: can't open input file ($1)." >&2
	exit 1
fi
##
## Input validation
## Too many users (in galaxy) try to plot a FASTQ file
## (without using the 'fastq statistics' tool first).
##
## gnuplot's error in that case is crypt, and support emails are annoying.
##
## try to detect FASTA/FASTQ input, and give a detailed, easy-to-understand warning.
##
##
AWK_FASTX_DETECTION='
NR==1 && $0 ~ /^>/ { fasta_id = 1 }
NR==1 && $0 ~ /^@/ { fastq_id = 1 }
NR==2 && $0 ~ /^[ACGT][ACGT]*$/ { nucleotides = 1 }
NR>3 { exit }
END { if ( fasta_id && nucleotides ) { print "FASTA" }
      if ( fastq_id && nucleotides ) { print "FASTQ" }
}'

INPUT_TYPE=$(awk "$AWK_FASTX_DETECTION" "$FILENAME")

if [ "x$INPUT_TYPE" = "xFASTA" ] ; then
#this doesn't even make sense: FASTA files don't contain any quality scores
cat>&2<<EOF
Error: It looks like your input file is a FASTA file.

This tool (fastq-quality-plot) can't use FASTA files directly - it requires a tabular text file conaining summary statistic about your FASTA file.

In Galaxy,
Please use the "Compute Quality Statistics (improved)" tool (in the "NGS: QC and Manipulation" category) to compute the quality statistics report, and then use this tool with the new statistics report.

On the command line,
Please use the "fastx_quality_stats" program to create the statistics report.
EOF
exit 1
fi
if [ "x$INPUT_TYPE" = "xFASTQ" ] ; then
cat>&2<<EOF
Error: It looks like your input file is a FASTQ file.

This tool (fastq-quality-plot) can't use FASTQ files directly - it requires a tabular text file conaining summary statistic about your FASTQ file.

In Galaxy,
Please use the "Compute Quality Statistics (improved)" tool (in the "NGS: QC and Manipulation" category) to compute the quality statistics report, and then use this tool with the new statistics report.

On the command line,
Please use the "fastx_quality_stats" program to create the statistics report.
EOF
exit 1
fi

##
## Even if this is not a FASTA/FASTQ file,
## users can still use incompatible input files.
## Try to detect it and abort with a warning.
AWK_VALID_STAT='NR==1 && $1=="cycle" && $2=="count" && $3=="min" { exit 2 } NR>1 { exit }'

awk "$AWK_VALID_STAT" "$FILENAME"
if [ $? -ne 2 ] ; then
cat>&2<<EOF
Error: Input file is not a valid statistics report.

This tool (fastq-quality-plot) requires a tabular text file conaining summary statistic about your FASTQ file.

In Galaxy,
Please use the "Compute Quality Statistics" tool (in the "NGS: QC and Manipulation" category) to compute the quality statistics report, and then use this tool with the new statistics report.

On the command line,
Please use the "fastx_quality_stats" program to create the statistics report.
EOF
exit 1
fi

GNUPLOTCMD="
$OUTPUTTERM
#set boxwidth 0.75 absolute
set size 1,1
#set style fill solid 1.00 border -1
set xlabel \"read position\"
set title \"Nucleotides distribution $TITLE\" 
set ylabel \"% of total (per read position)\" 
#set grid noxtics nomxtics ytics nomytics noztics nomztics \
# nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
#set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
set key outside right top vertical Left reverse enhanced autotitles columnhead nobox
set key invert samplen 4 spacing 1 width 0 height 0 
#set style histogram rowstacked 
#set style data histograms 
set ytics 10
set xtics 1
set yrange [ 0.00000 : 100.000 ] noreverse nowriteback

plot '$FILENAME' using (100.*column(13)/column(18)):xtic(1) title \"A\" with linespoints lt rgb \"#5050ff\", \
       '' using (100.*column(14)/column(18)) title \"C\" with linespoints lt rgb \"#e00000\", \
       '' using (100.*column(15)/column(18)) title \"G\" with linespoints lt rgb \"#00c000\", \
       '' using (100.*column(16)/column(18)) title \"T\" with linespoints lt rgb \"#e6e600\", \
       '' using (100.*column(17)/column(18)) title \"N\" with linespoints lt rgb \"pink\"
"

echo "$GNUPLOTCMD" | gnuplot > "$OUTPUTFILE"
