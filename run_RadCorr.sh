#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-04-30 13:02:57 trottar"
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
	echo " Avaliable Kinematics..."
	echo "                      theta = 11, 18, 30 deg"
	echo "                      polar = long, trans, unpol"
	echo
        exit 0
        ;;
	c) c_flag='true' ;;
	a) a_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

angle=$2
polar=$3
echo "Angle must be one of - [11 - 18 - 30]"
if [[ -z "$2" || ! "$angle" =~ 11|18|30 ]]; then # Check the 2nd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid angle..."
    while true; do
	echo ""
	read -p "Angle must be one of - [11 - 18 - 30] - or press ctrl-c to exit : " angle
	case $angle in
	    '');; # If blank, prompt again
	    '11'|'18'|'30') break;; # If a valid option, break the loop and continue
	esac
    done
fi
echo "Polarization must be one of - [long - trans - unpol]"
if [[ -z "$3" || ! "$polar" =~ long|trans|unpol ]]; then # Check the 2nd argument was provided and that it's one of the valid options
    echo ""
    echo "I need a valid polarization..."
    while true; do
	echo ""
	read -p "Polarization must be one of - [long - trans - unpol] - or press ctrl-c to exit : " polar
	case $polar in
	    '');; # If blank, prompt again
	    'long'|'trans'|'unpol') break;; # If a valid option, break the loop and continue
	esac
    done
fi

# Source root version
source /apps/root/6.28.06/setROOT_CUE-gcc10.2.0.sh # Carter used 6.30.0, g++ 12.3.0, but this seems to work fine

input_dataset_filename="data_sets_radiate_${angle}deg_${polar}"
radcorr_output_fiilename="radcor_${angle}deg_${polar}"
radiate_output_fiilename="radiated_model_${angle}deg_${polar}"

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
.L rad_corr.C
rad_corr("configs/${input_dataset_filename}.conf","output/${radcorr_output_fiilename}.dat")
radiate("configs/${input_dataset_filename}.conf","output/${radiate_output_fiilename}.dat")
EOF
fi
