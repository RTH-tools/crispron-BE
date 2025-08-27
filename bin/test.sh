#!/bin/bash


path=`pwd`

rm -rf test/
mkdir -p test || exit 1
cp test_org/seq.fasta test

./bin/CRISPRonBE.sh ABE test/seq.fasta test/ABE_test
if [[ -s "test/ABE_test/crispronABE_prediction.tsv" ]]; then
	echo ABE Test ok
else
	echo ABE Test failed
	exit 1
fi

./bin/CRISPRonBE.sh CBE test/seq.fasta test/CBE_test
if [[ -s "test/CBE_test/crispronCBE_prediction.tsv" ]]; then
	echo CBE Test ok
else
	echo CBE Test failed
	exit 1
fi

echo Test ok
echo "Please check the output in the following files and compare the contents"
echo "to that found in test_org"
echo "test/ABE_test/crispronABE_prediction.tsv"
echo "test/CBE_test/crispronCBE_prediction.tsv"
