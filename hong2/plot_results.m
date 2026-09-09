clear; clc; close all;

load('result_dinkelbach.mat');

residual_hist_abs = abs(residual_hist);
residual_hist_abs(isnan(residual_hist_abs)) = eps;

stats = plot_dinkelbach_convergence_ieee(residual_hist_abs, 1e-4, '');