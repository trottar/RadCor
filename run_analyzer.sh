#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-04-30 08:01:59 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#

# Flag definitions (flags: h, c)
while getopts 'hca' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./run_analyzer.sh -{flags} {variable arguments, see help}"
	echo
        echo "Description: "
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the RadCorr analysis..."
        echo "    -h, help"
        echo "    -c, compile CAnalyzer"
	echo "    -a, radiate function for all a1n d2n kinematics (11, 18, 30deg) and polarizations (unpol, long, trans)"
	echo
        exit 0
        ;;
	c) c_flag='true' ;;
	r) a_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

# Source root version
source /apps/root/6.22.06/setroot_CUE.bash

dataset_config_filename="data_example"
output_fiilename="radcor_out"

# Specific example
#dataset_config_filename="data_sets_radiate_11deg_long"
#output_fiilename="radiated_model_11deg_long"

# Compile CAnalyzer
if [[ $c_flag == "true" ]]; then
    echo
    echo
    echo "Compiling..."
    echo
    echo
    cd "CAnalyzer-master"
    make clean
    make
    cd "../"
fi

# Load ROOT macro
echo
echo
echo "Loading ROOT macro..."
echo
echo
cd "CAnalyzer-master/example/"
if [[ $a_flag == "true" ]]; then
    root -l <<EOF 
.L rad_corr.C
radiate_all()
EOF
else
    root -l <<EOF
.L rad_corr.C("configs/${dataset_config_filename}.conf","output/${output_fiilename}.dat")
rad_corr
radiate
EOF
fi
