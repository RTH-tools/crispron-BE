#!/bin/bash
BIN=`dirname -- "$( readlink -f -- "$0"; )";`
DIR=`readlink -f $BIN/..`

echo "testing python dependencies"
python3 $BIN/test.py || exit 1

echo "downloading crispron"
[[ -e $DIR/crispron-1.0.tar.gz ]] || wget https://rth.dk/resources/crispr/crispron/downloads/crispron-1.0.tar.gz -O $DIR/crispron-1.0.tar.gz || exit 1
[[ -d $DIR/crispron-1.0 ]] || tar xf $DIR/crispron-1.0.tar.gz -C $DIR || exit 1
cp -p $DIR/crispron-1.0/bin/DeepCRISPRon_eval.py $BIN


echo "downloading crisproff"
[[ -e $DIR/crisproff-1.1.2.tar.gz ]] || wget https://rth.dk/resources/crispr/crisproff/downloads/crisproff-1.1.2.tar.gz -O $DIR/crisproff-1.1.2.tar.gz || exit 1
[[ -d $DIR/crisproff-1.1.2 ]] || tar xf $DIR/crisproff-1.1.2.tar.gz  -C $DIR || exit 1
mkdir -p data/model
cp -p $DIR/crisproff-1.1.2/CRISPRspec_CRISPRoff_pipeline.py $BIN
cp -p $DIR/crisproff-1.1.2/energy_dics.pkl $DIR/data/model

echo "downloading crispron V0 models"
[[ -e $DIR/CRISPRon_models_V0-1.0.tgz ]] || wget https://rth.dk/resources/crispr/crispron/downloads/CRISPRon_models_V0-1.0.tgz -O $DIR/CRISPRon_models_V0-1.0.tgz || exit 1
[[ -d $DIR/data/CRISPRon_models_V0 ]] || tar xf $DIR/CRISPRon_models_V0-1.0.tgz -C $DIR || exit 1
[[ -d $DIR/data/CRISPRon_models_V0 ]] || mv $DIR/CRISPRon_models_V0-1.0 $DIR/data/CRISPRon_models_V0/

echo "downloading crispron-BE models"
[[ -e $DIR/CRISPRonBE_models-1.0.tgz ]] || wget https://rth.dk/resources/crispr/crispron-be/downloads/CRISPRonBE_models-1.0.tgz -O $DIR/CRISPRonBE_models-1.0.tgz || exit 1
[[ -d $DIR/data/CRISPRonBE_models ]] || tar xf $DIR/CRISPRonBE_models-1.0.tgz -C $DIR || exit 1
[[ -d $DIR/data/CRISPRonBE_models ]] || mv $DIR/CRISPRonBE_models-1.0 $DIR/data/CRISPRonBE_models  


echo "deleting $DIR/test"
[[ -d $DIR/test ]] && rm -rf $DIR/test

echo "testing CRISPRon-BE in $DIR/test"
$BIN/test.sh
