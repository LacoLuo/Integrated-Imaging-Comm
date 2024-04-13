clearvars
clc

load("BS_results/sample12/background_subtraction.mat")

RIS_DM = figure; imagesc(DM); colorbar
background_DM = figure; imagesc(DM_background); colorbar
background_subtracted_DM = figure; imagesc(BS_DM); colorbar
filtered_background_subtracted_DM = figure; imagesc(filtered_BS_DM); colorbar
binary_background_subtracted_DM = figure; imagesc(binary_BS_DM); colorbar
UE_only_DM = figure; imagesc(final_DM); colorbar