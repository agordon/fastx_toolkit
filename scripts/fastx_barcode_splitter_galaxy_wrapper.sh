#!/bin/sh

#
#This is a shell script wrapper for 'fastx_barcode_splitter.pl'
#
# 1. Output files are saved at a predefined location
#    (Which was made publicly accessible using apache)
#
# 2. 'fastx_barcode_splitter.pl' outputs a textual table.
#    This script turns it into pretty HTML with working URL
#    (so lazy users can just click on the URLs and get thier files)

BASEPATH="/media/sdb1/galaxy/barcode_splits/"
PUBLICURL="http://tango.cshl.edu/barcode_splits/"

BARCODE_FILE="$1"
FASTQ_FILE="$2"
LIBNAME="$3"
shift 3
# The rest of the parameters are passed to the split program

if [ "$LIBNAME" == "" ]; then
	echo "Usage: $0 [BARCODE FILE] [FASTQ FILE] [LIBRARY_NAME]" >&2
	exit 1
fi

#Sanitize library name, make sure we can create a file with this name
LIBNAME=${LIBNAME//\.gz/}
LIBNAME=${LIBNAME//\.txt/}
LIBNAME=${LIBNAME//[^[:alnum:]]/_}

if [ ! -r "$FASTQ_FILE" ]; then
	echo "Error: Input file ($FASTQ_FILE) not found!" >&2
	exit 1
fi
if [ ! -r "$BARCODE_FILE" ]; then
	echo "Error: barcode file ($BARCODE_FILE) not found!" >&2
	exit 1
fi

PREFIX="$BASEPATH"`date "+%Y-%m-%d_%H%M__"`"${LIBNAME}__"
SUFFIX=".txt"

RESULTS=`zcat -f "$FASTQ_FILE" | fastx_barcode_splitter.pl --bcfile "$BARCODE_FILE" --prefix "$PREFIX" --suffix "$SUFFIX" "$@"`
if [ $? != 0 ]; then
	echo "error"
fi

#
# Convert the textual tab-separated table into simple HTML table,
# with the local path replaces with a valid URL
echo "<html><body><table border=1>"
echo "$RESULTS" | sed "s|$BASEPATH|$PUBLICURL|" | sed '
i<tr><td>
s|\t|</td><td>|g
s|http.*|<a href="&">&<\/a>|
a<\/td><\/tr>
'
echo "<p><b>Copy these files to your local computer, as they will be soon deleted.</b>"
echo "</table></body></html>"
