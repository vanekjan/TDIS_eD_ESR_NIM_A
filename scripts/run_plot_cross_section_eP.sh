#!/bin/bash

#Run plotting macro
root -b -q ./macros/Plot_cross_section_EICRecon_eP.cxx+\(1,0,\"${1}\"\)

#Copy and rename figures for NIM A paper
cp ./plots/eP/cross_section/sigma_red_D_Q2bin_0_x_bin_1.png ./plots/pub/NIM_A/TDIS_eD_sigma_red_D_eP.png
cp ./plots/eP/cross_section/sigma_red_D_Q2bin_0_x_bin_1.pdf ./plots/pub/NIM_A/TDIS_eD_sigma_red_D_eP.pdf

cp ./plots/eP/cross_section/sigma_red_p_corr_Q2bin_0_x_bin_1.png ./plots/pub/NIM_A/TDIS_eD_sigma_red_n_eP.png
cp ./plots/eP/cross_section/sigma_red_p_corr_Q2bin_0_x_bin_1.pdf ./plots/pub/NIM_A/TDIS_eD_sigma_red_n_eP.pdf

cp ./plots/eP/cross_section/sigma_red_p_vs_x_Q2bin_0.png ./plots/pub/NIM_A/TDIS_eD_sigma_red_p_vs_x.png
cp ./plots/eP/cross_section/sigma_red_p_vs_x_Q2bin_0.pdf ./plots/pub/NIM_A/TDIS_eD_sigma_red_p_vs_x.pdf


