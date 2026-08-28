#!/bin/bash

#Run plotting macro

#First user argument (${1}) is path to the input file.
#Second user argument (${2}) is the collision energy:
#9x130 GeV: 9130
#10x130 GeV: 10130

#Example usege (run from repo base folder):
# ./sripts/run_plot_cross_section_eN.sh ./production/output_RC_eN_test.root 9130

root -b -q ./macros/Plot_cross_section_EICRecon_eN.cxx+\(1,0,\"${1}\",${2}\)

#Copy and rename figures for NIM A paper
cp ./plots/eN/cross_section/sigma_red_D_Q2bin_0_x_bin_1.png ./plots/pub/NIM_A/TDIS_eD_sigma_red_D_eN.png
cp ./plots/eN/cross_section/sigma_red_D_Q2bin_0_x_bin_1.pdf ./plots/pub/NIM_A/TDIS_eD_sigma_red_D_eN.pdf

cp ./plots/eN/cross_section/sigma_red_n_corr_Q2bin_0_x_bin_1.png ./plots/pub/NIM_A/TDIS_eD_sigma_red_n_eN.png
cp ./plots/eN/cross_section/sigma_red_n_corr_Q2bin_0_x_bin_1.pdf ./plots/pub/NIM_A/TDIS_eD_sigma_red_n_eN.pdf

cp ./plots/eN/cross_section/sigma_red_n_vs_x_Q2bin_0.png ./plots/pub/NIM_A/TDIS_eD_sigma_red_n_vs_x.png
cp ./plots/eN/cross_section/sigma_red_n_vs_x_Q2bin_0.pdf ./plots/pub/NIM_A/TDIS_eD_sigma_red_n_vs_x.pdf
