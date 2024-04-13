clearvars

load("./results/avg_top_k_bf_gain_OSF1x1.mat")
OSF_1x1_topk_bf_gain = avg_top_k_bf_gain; clear avg_top_k_bf_gain;
load("./results/avg_top_k_bf_gain_OSF4x4.mat")
OSF_4x4_topk_bf_gain = avg_top_k_bf_gain; clear avg_top_k_bf_gain;

load("./results/opt_CB_bf_gain_OSF1x1.mat") 
OSF_1x1_opt_bf_gain = avg_opt_CB_bf_gain; clear avg_opt_CB_bf_gain;
load("./results/opt_CB_bf_gain_OSF4x4.mat") 
OSF_4x4_opt_bf_gain = avg_opt_CB_bf_gain; clear avg_opt_CB_bf_gain;

EG_bf_gain = zeros(size(OSF_4x4_topk_bf_gain));

set_default_plot;

yline(0, '--',  Color='black', LineWidth=3, DisplayName='Equal-gain beamforming');
hold on;
yline(OSF_1x1_opt_bf_gain, '--', LineWidth=3,  Color='#F2606A', DisplayName='Exhaustive search, OSF=1');
hold on;
yline(OSF_4x4_opt_bf_gain, '--', LineWidth=3,  Color='#298fc2', DisplayName='Exhaustive search, OSF=4');
hold on;
plot(OSF_1x1_topk_bf_gain, Marker="o", MarkerSize=10, MarkerFaceColor='white',  Color='#F2606A', DisplayName='Proposed beam selection, OSF=1');
hold on;
plot(OSF_4x4_topk_bf_gain, Marker="o", MarkerSize=10, MarkerFaceColor='white',  Color='#298fc2', DisplayName='Proposed beam selection, OSF=4');
grid on;
box on;

set(gca, 'LooseInset', get(gca, 'TightInset'));

legend(Location='southeast')
ylim([-9, 0.3])

xlabel('Number of refinement beams k')
ylabel('Normalized beamforming gain (dB)')

