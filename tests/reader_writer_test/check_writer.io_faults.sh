#!/usr/bin/env bash

##
## Check handling of POSIX faults
## requires fiu-run and libfiu (http://blitiri.com.ar/p/libfiu/)

FIURUN=$(which fiu-run)

if [ -z "$FIURUN" ]; then
	echo "Running these test requires libfiu and fiu-run." >&2
	echo "Download & install them from http://blitiri.com.ar/p/libfiu/" >&2
	exit 1;
fi

PROG=../../src/fastx_reader_writer_test/fastx_reader_writer_test

if [ ! -x "$PROG" ]; then
	echo "Error: can't find executable '$PROG'" >&2
	exit 1;
fi

##
## Uncompress input files
##
zcat long.fa.gz > long.fa || { echo "Failed to uncompress test file long.fa.gz" >&2 ; exit 1 ; }
zcat long.fq.gz > long.fq || { echo "Failed to uncompress test file long.fq.gz" >&2 ; exit 1 ; }

# Check read failures
#


# Fail 10% of the times - should happen in multiple stages of the FASTA reader
for i in $(seq 1 10) ; do
	echo "$i"
	fiu-run -x -e posix/io/rw/write -p 10 $PROG -i long.fa > /dev/null
	fiu-run -x -e posix/io/rw/write -p 10 $PROG -i long.fq > /dev/null
done
