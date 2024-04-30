#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-04-30 12:38:54 trottar"
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
        echo "./set_Model.sh -{flags} {variable arguments, see help}"
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
	a) a_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

# Source root version
source /apps/root/6.28.06/setROOT_CUE-gcc10.2.0.sh # Carter used 6.30.0, g++ 12.3.0, but this seems to work fine

# Compile CAnalyzer
if [[ $c_flag == "true" ]]; then
    echo
    echo
    echo "Compiling..."
    echo
    echo
    cd "YoniCode_Retry_4Dave_2/Carter/"
    make clean
    make
fi

if [[ $j_flag == "true" ]]; then
    echo
    echo
    echo "Generating cross section using JAM for theta of 11, 18, 30..."
    echo
    echo
    ./xs_gen_jam
else
    echo
    echo
    echo "Generating cross section using JAM for theta of 11, 18, 30..."
    echo
    echo    
    ./xs_gen_dis6
fi

