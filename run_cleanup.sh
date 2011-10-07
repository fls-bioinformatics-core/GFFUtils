#!/bin/sh -e
#
# cleanup.sh: clean up a GFF file
#
# Runs GFFcleaner multiple times to perform each stage of clean up
# automatically
#
if [ "$#" -lt 1 ] ; then
    echo "Usage: `basename $0` <gff_file> [<mapping_file>]"
    exit
fi
# Initialise
SCRIPT_DIR=`dirname $0`
gff=$1
gff_base=${gff%.*}
mapping=$2
#
# Do prepend
echo "============= PREPEND STEP ================"
gff_in=$gff
gff_out=${gff_in%.*}_clean.gff
next_gff=${gff_base}_prepend.gff
cmd="${SCRIPT_DIR}/GFFcleaner.py --prepend=Sbay_ $gff_in"
echo Running $cmd
$cmd
if [ -f $gff_out ] ; then
    echo Output moved to $next_gff
    /bin/mv $gff_out $next_gff
fi
#
# Do clean
echo "============== CLEAN STEP ================="
gff_in=$next_gff
gff_out=${gff_in%.*}_clean.gff
next_gff=${gff_base}_cleaned.gff
cmd="${SCRIPT_DIR}/GFFcleaner.py --clean $gff_in"
echo Running $cmd
$cmd
if [ -f $gff_out ] ; then
    echo Output moved to $next_gff
    /bin/mv $gff_out $next_gff
fi
#
# Do steps involving mapping file
if [ ! -z "$mapping" ] ; then
    # Remove duplicates
    echo "=========== RESOLVE DUPLICATES ============="
    gff_in=$next_gff
    gff_out=${gff_in%.*}_clean.gff
    next_gff=${gff_base}_resolved_duplicates.gff
    cmd="${SCRIPT_DIR}/GFFcleaner.py --report-duplicates --resolve-duplicates=${mapping} $gff_in"
    echo Running $cmd
    $cmd
    if [ -f $gff_out ] ; then
	echo Output moved to $next_gff
	/bin/mv $gff_out $next_gff
    fi
    #
    # Insert missing
    echo "========= INSERT MISSING GENES ============"
    gff_in=$next_gff
    gff_out=${gff_in%.*}_clean.gff
    next_gff=${gff_base}_inserted_missing.gff
    cmd="${SCRIPT_DIR}/GFFcleaner.py --insert-missing=${mapping} $gff_in"
    echo Running $cmd
    $cmd
    if [ -f $gff_out ] ; then
	echo Output moved to $next_gff
	/bin/mv $gff_out $next_gff
    fi
fi
#
# Finish
echo "================ FINISHED =================="
final_gff=${gff_base}_final.gff
/bin/mv $next_gff $final_gff
echo Final output file: $final_gff
##
#