#!/bin/bash
BINDIR=`dirname $0`
DATADIR=$BINDIR/../data
EDITOR=$1
CRISPRON_FASTA=$2
OUTDIR=$3
weight=$4

if [[ -z $EDITOR ]]; then
	echo "Needs three arguments to run e.g. $0 ABEorCBE test/seq.fasta test/outdir" 1>&2
	exit 1
fi

if [[ -z $OUTDIR ]]; then
	echo "Needs three arguments to run e.g. $0 ABEorCBE test/seq.fasta test/outdir" 1>&2
	exit 1
fi

mkdir -p $OUTDIR || exit 1

if [[ ! -s $CRISPRON_FASTA ]]; then
	echo "Needs a fasta file as second input parameter to run e.g. $0 ABEorCBE indir/test.fasta outdir" 1>&2
	exit 1
fi

which python3 || exit 1
RNAfold --version || exit 1

$BINDIR/get_30mers_from_fa.py  -f ${CRISPRON_FASTA} -m $OUTDIR/30mers.fa -g $OUTDIR/23mers.fa || exit 1

if [[ ! -s $OUTDIR/30mers.fa ]]; then
	echo failed to get target sequences 1>&2
	exit 1
fi

echo "#Running CRISPROff pipeline"
$BINDIR/CRISPRspec_CRISPRoff_pipeline.py \
	--guides $OUTDIR/23mers.fa \
	--specificity_report $OUTDIR/CRISPRspec.tsv \
	--guide_params_out $OUTDIR/CRISPRparams.tsv \
	--duplex_energy_params $DATADIR/model/energy_dics.pkl \
	--no_azimuth || exit 1

echo "#Running CRISPRon prediction"
$BINDIR/DeepCRISPRon_eval.py $OUTDIR $OUTDIR/30mers.fa $OUTDIR/CRISPRparams.tsv $DATADIR/CRISPRon_models_V0/best/*/  2>&1 | grep -v tracing

if [ -z "$weight" ]; then
    if [[ $EDITOR == 'ABE' ]]; then
        weight='0.5-0.5-0-0-0'
    else
        weight='0.5-0.5-0'
    fi
fi

echo "#Running CRISPRon-BE evaluation"
$BINDIR/DeepCRISPRonBE_eval.py $EDITOR $OUTDIR/30mers.fa $OUTDIR/CRISPRparams.tsv $OUTDIR/crispron.csv $DATADIR/CRISPRonBE_models $OUTDIR $weight 2>&1 | grep -v tracing

rm -f temp_grna_id_ss.ps

if [[ -e "$OUTDIR/crispron${EDITOR}_prediction.tsv" ]]; then
    echo "#All done" 1>&2
else
    echo CRISPRon-BE failed, please provide the suitable input sequence.
fi
