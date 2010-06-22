#!/usr/bin/env bash

PROG=../../src/fastx_reader_writer_test/fastx_reader_writer_test

if [ ! -x "$PROG" ]; then
	echo "Error: can't find executable '$PROG'" >&2
	exit 1;
fi

for i in good* ; do
	## Check with STDIN
	cat "$i" | $PROG > /dev/null
	if (( $? )) ; then
		#Program failed, report and stop
		echo "Test failed to for '$i'" >&2
		exit 1
	else
		#Program succeeded, all OK
		DUMMY=
	fi

	## Check with Inputfile
	$PROG -i "$i" > /dev/null
	if (( $? )) ; then
		#Program failed, report and stop
		echo "Test failed to for '$i'" >&2
		exit 1
	else
		#Program succeeded, all OK
		DUMMY=
	fi
done

for i in bad* ; do
	## Check with STDIN
	cat "$i" | $PROG > /dev/null
	if (( $? )) ; then
		#Program failed, which is expected
		DUMMY=
	else
		#Program succeeded, which is unexpected
		#rerun, show STDERR to the user
		cat "$i" | $PROG > /dev/null
		echo "Test failed for '$i' - program should have failed." >&2
		exit 1
	fi

	## Check with input file
	$PROG -i "$i" > /dev/null
	if (( $? )) ; then
		#Program failed, which is expected
		DUMMY=
	else
		#Program succeeded, which is unexpected
		#rerun, show STDERR to the user
		cat "$i" | $PROG > /dev/null
		echo "Test failed for '$i' - program should have failed." >&2
		exit 1
	fi
done
