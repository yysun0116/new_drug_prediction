#!/bin/bash
set -e

# Required

if [ -z $IN_FILE ]; then
    echo "Err: missing input parameters"
    exit 1
fi

if [ -z $IN_SMILES_DIR ] ; then
    echo "Err: missing input parameters"
    exit 1
fi

if [ -z $OUT_DIR ] ; then
    echo "Err: missing output parameters"
    exit 2
fi


# Make script

mkdir -p $OUT_DIR
echo -n "python3 /scDrug/script/drug_response_prediction_new_drug.py -o ${OUT_DIR} -i ${IN_FILE} -i_smiles ${IN_SMILES_DIR}" >> tmp.sh


# Execute
/usr/bin/time -f "scDrug sensitivity prediction on new drugs mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" \
    bash tmp.sh
rm tmp.sh
