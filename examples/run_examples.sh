#!/bin/bash
#
# Run examples for GFFUtils
#
# Assertion functions for tests
_PASSED=0
_FAILED=0
function report_tests {
    # Report summary of tests passed and failed
    local n_tests=$((_PASSED+_FAILED))
    echo "---------------------------------------------------------"
    echo "Ran $n_tests tests: $_PASSED passed, $_FAILED failed"
    if [ $_FAILED -ne 0 ] ; then
	return 1
    else
	return 0
    fi
}
function run_test {
    # Run a command and check outputs
    # Takes following arguments:
    # --command CMD: specify the command to execute
    # --expected FILES: list of file names which should
    #                have reference versions to compare
    #                against
    # --must_exist FILES: list of file names which should
    #                exist after the command has run
    # --status INT: exit code to check (if command was
    #                run externally)
    # --strip-paths: strip leading paths found inside
    #                reference and output files when
    #                checking content matches
    local test_name=$1
    local command=
    local expected_outputs=
    local check_exists=
    local exit_status=
    local working_dir=
    local test_status=
    local strip_paths=
    # Collect arguments
    shift
    while [ $# -gt 0 ] ; do
	case $1 in
	    --command)
		command=$2
		shift
		;;
	    --expected)
		expected_outputs=$2
		shift
		;;
	    --must_exist)
		check_exists=$2
		shift
		;;
	    --status)
		exit_status=$2
		shift
		;;
	    --strip-paths)
		strip_paths=$1
		;;
	    *)
		echo "$test_name: SKIPPED (unrecognised test argument '$1')"
		return
		;;
	esac
	shift
    done
    echo "---------------------------------------------------------"
    echo test_name: $test_name
    echo command: $command
    echo expected_outputs: $expected_outputs
    echo check_exists: $check_exists
    echo exit_status: $exit_status
    echo strip_paths: $strip_paths
    echo PWD: $(pwd)
    # If command supplied then run it
    if [ ! -z "$command" ] ; then
	working_dir=$(mktemp -d --tmpdir=$(pwd))
	echo working_dir: $working_dir
	cd $working_dir
	echo "Running command"
	$command 1>STDOUT 2>STDERR
	exit_status=$?
	echo "Exit status $exit_status"
    fi
    # Check exit status
    if [ ! -z "$exit_status" ] ; then
	if [ $exit_status -ne 0 ] ; then
	    echo Failed exit status check
	    test_status=FAILED
	fi
    fi
    # Compare expected outputs
    for f in $expected_outputs ; do
	assert_equal $REF_DATA/ref_$f $f $strip_paths
	if [ $? -ne 0 ] ; then
	    echo Failed output comparison check
	    test_status=FAILED
	fi
    done
    # Check existence
    for f in $check_exists ; do
	if [ ! -e $f ] ; then
	    echo "$f: missing"
	    echo Failed output existence check
	    test_status=FAILED
	fi
    done
    # Set test status if no failures
    if [ -z "$test_status" ] ; then
	test_status=OK
    fi
    echo test_status: $test_status
    # Report logs from failed job
    if [ $test_status == FAILED ] ; then
	for f in STDOUT STDERR ; do
	    if [ -e $f ] ; then
		echo "===== $test_name: $f ====="
		cat $f
	    fi
	done
    fi
    # Clean up any working area
    if [ ! -z "$working_dir" ] ; then
	cd ..
	#rm -rf $working_dir
    fi
    # Test counts
    case $test_status in
	OK)
	    _PASSED=$((_PASSED+1))
	    ;;
	FAILED)
	    _FAILED=$((_FAILED+1))
	    ;;
    esac
    # Finish
    echo "---------------------------------------------------------"
    echo "TEST: $test_name: $test_status"
}
function assert_equal {
    # Check two files are the same
    local strip_paths=
    if [ "$3" == "--strip-paths" ] ; then
	strip_paths=yes
    fi
    if [ ! -e $1 ] ; then
	echo "$1: missing reference data"
	return 1
    elif [ ! -e $2 ] ; then
	echo "$2: missing"
	return 1
    fi
    if [ -z "$strip_paths" ] ; then
	old=$1
	new=$2
    else
	tmpdir=$(mktemp -d)
	old=$tmpdir/old
	sed 's,/.*/,,g' $1 >$old
	new=$tmpdir/new
	sed 's,/.*/,,g' $2 >$new
    fi
    diff -q $old $new
    if [ $? -ne 0 ] ; then
	echo "$2: doesn't match reference data:"
	diff $1 $2
	return 1
    else
	return 0
    fi
}
#
# Initialise and set up dir for test outputs
TEST_DIR=$(dirname $0)
if [ "$TEST_DIR" == "." ] ; then
    TEST_DIR=$(pwd)
elif [ -z "$(echo $TEST_DIR | grep ^/)" ] ; then
    TEST_DIR=$(pwd)/$TEST_DIR
fi
DATA_DIR=$TEST_DIR/data
REF_DATA=$TEST_DIR/ref-data
if [ ! -d test-output ] ; then
    mkdir test-output
else
    rm -rf test-output/*
fi
cd test-output
#
# gtf2bed.py
run_test "gtf2bed" \
    --expected "mm10.bed" \
    --command "gtf2bed -o mm10.bed ${DATA_DIR}/mm10.gtf"
#
# gff_annotation_extractor with GFF input
run_test "gff_annotation_extractor with GFF input" \
    --expected "test_dicty.out test_dicty_stats.out" \
    --command "gff_annotation_extractor -o test_dicty.out --htseq-count ${DATA_DIR}/dicty.gff ${DATA_DIR}/dicty_htseq_counts.dimA.txt ${DATA_DIR}/dicty_htseq_counts.AXA4.txt" \
    --strip-paths
#
# gff_annotation_extractor with GTF input
run_test "gff_annotation_extractor with GTF input" \
    --expected "test_mm10.out test_mm10_stats.out" \
    --command "gff_annotation_extractor -o test_mm10.out --htseq-count ${DATA_DIR}/mm10.gtf ${DATA_DIR}/mm10_htseq_counts.txt" \
    --strip-paths
#
# gff_annotation_extractor with GTF input with -i gene_name
run_test "gff_annotation_extractor with GTF input using 'gene_name' for ID" \
    --expected "test_mm10_gencode_vM5.out test_mm10_gencode_vM5_stats.out" \
    --command "gff_annotation_extractor -o test_mm10_gencode_vM5.out --htseq-count -i gene_name ${DATA_DIR}/mm10_gencode_vM5.gtf ${DATA_DIR}/mm10_day6_s1_htseq_counts.txt" \
    --strip-paths
#
# GFFcleaner
run_test "gff_cleaner" \
    --expected "dicty_cleaned.gff" \
    --command "gff_cleaner -o dicty_cleaned.gff ${DATA_DIR}/dicty.gff"
#
# gtf_extract
run_test "gtf_extract" \
    --expected "mm10_exons.gtf" \
    --command "gtf_extract -o mm10_exons.gtf -f exon ${DATA_DIR}/mm10.gtf"
#
# Finished
report_tests
exit $?
##
#
