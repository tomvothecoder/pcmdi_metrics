#!/bin/bash

source $CONDA_SOURCE
conda activate $CONDA_ENV_ROOT/_CMEC_pcmdi_metrics

cd $CMEC_WK_DIR

tmp_param=$CMEC_WK_DIR/mean_climate_param.py

python $CMEC_CODE_DIR/../scripts/pmp_param_generator.py $CMEC_CONFIG_DIR/cmec.json $tmp_param "mean_climate"

if [[ $? = 0 ]]; then
    mean_climate_driver.py -p $tmp_param

    # write output.json
    python $CMEC_CODE_DIR/mean_climate_output.py
else
    echo "Failure in PMP/mean_climate parameter file generation"
fi