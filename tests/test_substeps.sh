#! /bin/sh

. ./compat.sh

$PORTCULLIS prep -o temp/prep_test ${data}/spombe.III.fa ${data}/spombe.gsnap.III.25K.bam
$PORTCULLIS junc -o temp/junc_test/portcullis temp/prep_test
$PORTCULLIS filt -o temp/filt_test/portcullis temp/prep_test temp/junc_test/portcullis.junctions.tab
