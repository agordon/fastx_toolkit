#!/bin/sh

#
# Arguments check and suage information
#
SRC="."
DEST="$1"
if [ -z "$DEST" ]; then
cat<<EOF

FASTX-toolkit Galaxy Installation script.
	
This script copies the FASTX-Toolkit files into the specified Galaxy directory.
	
Usage: $0 [GALAXY-DIRECTORY]
	
	GALAXY-DIRECTORY - root directory of the Galaxy server.
	
EOF
	exit
fi


echo 
echo "FASTX-toolkit Galaxy Installation script."
echo

#
# Sanity checks for the specified galaxy directory
#
echo -n "Checking Galaxy destination directory..."
[ -d "$DEST" ] ||
    { echo "Error: directory '$DEST' does not exist!" ; exit 1 ; }
    
[ -r "$DEST/tool_conf.xml" ] ||
    { echo "Error: file '$DEST/tool_conf.xml' does not exist! (is '$DEST' the root of the Galaxy server?)" ; exit 1 ; }
    
for subdir in tools tool-data test-data static; do
	[ -d "$DEST/$subdir" ] ||
		{ echo "Error: sub-directory '$DEST/$subdir' does not exist! (is '$DEST' the root of the Galaxy server?)" ; exit 1 ; }
done
echo "ok"

#
# Sanity checks for the FASTX-toolkit files
#
echo -n "Checking FASTX-toolkit source directory..."
[ -r "$SRC/galaxy/fastx_toolkit_conf.xml" ] ||
    { echo "Error: file '$SRC/galaxy/fastx_toolkit_conf.xml' does not exist! (is '$SRC' the root of FASTX-toolkit ?)" ; exit 1 ; }
    
for subdir in tools tools/fastx_toolkit tool-data test-data static static/fastx_icons; do
	[ -d "$SRC/galaxy/$subdir" ] ||
		{ echo "Error: sub-directory '$SRC/galaxy/$subdir' does not exist! (is '$SRC' the root of FASTX-toolkit?)" ; exit 1 ; }
done
echo "ok"


#
# Copy FASTX-Toolkit files into Galaxy server
#
echo -n "Creating static/fastx_icons directory..."
mkdir -p "$DEST/static/fastx_icons" || exit 1 ;
echo "OK"

echo -n "Copying static/fastx_icons..."
cp $SRC/galaxy/static/fastx_icons/*.png "$DEST/static/fastx_icons" || exit 1 ;
echo "OK"

echo -n "Copying test-data files..."
cp $SRC/galaxy/test-data/fast* "$DEST/test-data" || exit 1 ;
echo "OK"

echo -n "Copying tool-data files..."
cp $SRC/galaxy/tool-data/fastx_clipper_sequences.txt "$DEST/tool-data/" || exit 1;
echo "OK"

echo -n "Creaing tools/fastx_toolkit directory..."
mkdir -p "$DEST/tools/fastx_toolkit" || exit 1;
echo "OK"

#
# Be extra careful when copying the XML files - 
# Ask the user for confirmation if the XML files already exists
# (so that if they were changed, they will not be blindly overwriten)
echo "==="
echo "=== NOTE:"
echo "==="
echo "If the FASTX-toolkit XML files already exist on your galaxy server,"
echo "You will be prompted to confirm overwriting them."
echo "If you have made any changes to the XML files, DO NOT overwrite your files."
echo 
echo -n "Copying FASTX-toolkit XML tool configuration..."
cp -i $SRC/galaxy/tools/fastx_toolkit/*.{xml,sh} "$DEST/tools/fastx_toolkit"
echo "ok"



#
# Instruct the user what to do next
#
cat<<EOF
FASTX-toolkit files copied to your galaxy server directory.

Additionally, you'll need to make the following manual configurations:

1. Add the content of 
	$SRC/galaxy/fastx_toolkit_conf.xml 
   to
	$DEST/tool_conf.xml
	
2. Update the adapters file:

	$DEST/tool-data/fastx_clipper_sequences.txt
	
   And add valid adapters/linkers.
	
3. Edit "fastx_barcode_splitter_galaxy_wrapper.sh", change 
   The two variables BASEPATH and PUBLICURL to valid path/URL.
   See README for detailed explanation (under the
   "Special configuration for Barcode-Splitter" section).
   
EOF

