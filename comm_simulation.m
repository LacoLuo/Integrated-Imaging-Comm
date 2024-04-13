clearvars
clc

%  Azimuth and zenith angles towards UE
UE_theta_list = [95.5196, 95.1679, 95.3338, 96.6688, 93.5047, 94.0724, 94.5626, 93.7828, ...
                           94.8281, 93.9211, 96.5415, 96.6761, 96.0453, 94.6824, 94.5220, 94.2879, ...
                           94.1574, 96.1631, 97.1122, 97.2974, 96.6688, 95.8848, 95.6095, 94.6792, ...
                           97.0680, 97.1229, 97.6079, 97.9139, 97.2768, 97.0056, 95.9923, 95.4320];
UE_phi_list = 180 - [65.5650, 72.2406, 79.4517, 85.3160, 94.6840, 101.3710, 107.7594, 114.4350, ...
                                 62.7994, 70.7060, 77.8111, 86.1649, 93.8351, 102.1889, 109.2940, 116.5213, ...
                                 60.1644, 66.9958, 76.1905, 85.3160, 94.6840, 103.8095, 112.2772, 120.4738, ...
                                 55.2894, 63.4787, 73.0180, 84.4693, 96.3750, 106.9820, 116.5213, 125.8481];

% RIS codebook
oversampling_x = 1;
if oversampling_x == 1
    load("./codebooks/UPA_codebook_40x40_OSF_1x1.mat") % F_CB = [# ant., # codes]
elseif oversampling_x == 2
    load("./codebooks/UPA_codebook_40x40_OSF_2x2.mat") % F_CB = [# ant., # codes]
elseif oversampling_x == 4
    load("./codebooks/UPA_codebook_40x40_OSF_4x4.mat") % F_CB = [# ant., # codes]
end

% Antenna parameters for a UPA RIS structure
RIS_element_spacing = 0.5;
kd_RIS = 2 * pi *  RIS_element_spacing;
M_RIS_idx = antenna_channel_map(40, 1, 40, 0);
Mx = 40;
ang_conv = pi / 180;

num_samples = 96;
start = 1;
k = 25;
top_k_acc = zeros(num_samples-(start-1), k);
top_k_bf_gain = zeros(num_samples-(start-1), k);
opt_CB_bf_gain = zeros(num_samples-(start-1), 1);

for sample=start:1:num_samples
    disp(sample)
    current_power_WI = -30; % Current power in WI in dBW (Now it is set to 0 dBm (-30 dBW))

    load(strcat('./channels/sample', int2str(sample), '/APtoRIS_mat/scene_0_TX1.mat'))
    APtoRIS_params.num_paths = numel(channels{1}.paths.phase);
    APtoRIS_params.DoA_theta = channels{1}.paths.DoA_theta;
    APtoRIS_params.DoA_phi = channels{1}.paths.DoA_phi;
    % Subtract the transmit power to get accurate paths gain (Path_gain_dB = P_rec_dBm - P_transmit_dBm)
    APtoRIS_params.power = 10.^( 0.1 * ( double (channels{1}.paths.power) - (current_power_WI + 30) ) );
    APtoRIS_params.phase = channels{1}.paths.phase;
    clear channels
    
    load(strcat('./channels/sample', int2str(sample), '/RIStoUE_mat/scene_0_TX4.mat'))
    RIStoUE_params.num_paths = numel(channels{1}.paths.phase);
    RIStoUE_params.DoD_theta = channels{1}.paths.DoD_theta;
    RIStoUE_params.DoD_phi = channels{1}.paths.DoD_phi;
    % Subtract the transmit power to get accurate paths gain (Path_gain_dB = P_rec_dBm - P_transmit_dBm)
    RIStoUE_params.power = 10.^( 0.1 * ( double (channels{1}.paths.power) - (current_power_WI + 30) ) );
    RIStoUE_params.phase = channels{1}.paths.phase;
    clear channels
    
    % Far-field array response vector from AP to RIS
    APtoRIS_array_response_phases = sqrt(-1) * kd_RIS * ... 
                                                        [sind(APtoRIS_params.DoA_theta).*cosd(APtoRIS_params.DoA_phi); ...
                                                        sind(APtoRIS_params.DoA_theta).*sind(APtoRIS_params.DoA_phi); ...
                                                        cosd(APtoRIS_params.DoA_theta)];
    APtoRIS_array_response = exp(M_RIS_idx * APtoRIS_array_response_phases);
    
    % Far-field array response vector from RIS to UE
    RIStoUE_array_response_phases = sqrt(-1) * kd_RIS * ... 
                                                        [sind(RIStoUE_params.DoD_theta).*cosd(RIStoUE_params.DoD_phi); ...
                                                        sind(RIStoUE_params.DoD_theta).*sind(RIStoUE_params.DoD_phi); ...
                                                        cosd(RIStoUE_params.DoD_theta)];
    RIStoUE_array_response = exp(M_RIS_idx * RIStoUE_array_response_phases);
    
    % AP-RIS channel generation
    APtoRIS_complex_channel_gain = sqrt(APtoRIS_params.power).*exp( sqrt(-1) * APtoRIS_params.phase * ang_conv);
    APtoRIS_complex_channel = APtoRIS_array_response * APtoRIS_complex_channel_gain.';
    
    % RIS-UE channel generation
    RIStoUE_complex_channel_gain = sqrt(RIStoUE_params.power).*exp( sqrt(-1) * RIStoUE_params.phase * ang_conv);
    RIStoUE_complex_channel = RIStoUE_array_response * RIStoUE_complex_channel_gain.';
    
    % RIS interaction vector 
    % Codebook-based beam steering 
    AP_theta = APtoRIS_params.DoA_theta(1);
    AP_phi = APtoRIS_params.DoA_phi(1);
    AP_array_response_phases = sqrt(-1) * kd_RIS * ... 
                                                        [sind(AP_theta).*cosd(AP_phi); ...
                                                        sind(AP_theta).*sind(AP_phi); ...
                                                        cosd(AP_theta)];
    AP_RIS_interaction_vector = conj( exp( M_RIS_idx * AP_array_response_phases ));
    
    UE_theta = UE_theta_list(mod(sample-1, 32)+1);
    %UE_theta = UE_theta_list(sample);
    UE_phi = UE_phi_list(mod(sample-1, 32)+1);
    %UE_phi = UE_phi_list(sample);
    UE_array_response_phases = sqrt(-1) * kd_RIS * ... 
                                                        [sind(UE_theta).*cosd(UE_phi); ...
                                                        sind(UE_theta).*sind(UE_phi); ...
                                                        cosd(UE_theta)];
    UE_RIS_interaction_vector = conj( exp( M_RIS_idx * UE_array_response_phases ));
    
    RIS_interaction_vector = AP_RIS_interaction_vector.*UE_RIS_interaction_vector;
    
    % Selected the beam from codebook
    CB_similarity = abs( conj(F_CB).' * RIS_interaction_vector );
    [~, CB_idx] = max(CB_similarity);
    CB_RIS_interaction_vector = F_CB(:, CB_idx);

    [sorted_CB_sim, sorted_CB_sim_idx] = sort(CB_similarity, 'descend');
    
    % Equal-gain beamforming
    EG_RIS_interaction_vector = exp( -sqrt(-1) * angle( APtoRIS_complex_channel.*RIStoUE_complex_channel ) );
    
    % Random beamforming
    random_RIS_interaction_vector = exp( sqrt(-1) * 2*pi*randn(size(RIS_interaction_vector)) );
    
    % Channel generation
    channel = (APtoRIS_complex_channel.*RIStoUE_complex_channel).' * RIS_interaction_vector;
    EG_channel = (APtoRIS_complex_channel.*RIStoUE_complex_channel).' * EG_RIS_interaction_vector;
    random_channel = (APtoRIS_complex_channel.*RIStoUE_complex_channel).' * random_RIS_interaction_vector;
    CB_channel = double(APtoRIS_complex_channel.*RIStoUE_complex_channel).' * CB_RIS_interaction_vector;

    % Finding the optimal beam in the codebook
    all_CB_channels = abs(double(APtoRIS_complex_channel.*RIStoUE_complex_channel).' * F_CB);
    [opt_CB_channel, opt_CB_idx] = max(all_CB_channels);
    
    % Calculate channel gain
    sensing_aided_channel_gain = 10*log10(abs(channel)^2);
    EG_channel_gain = 10*log10(abs(EG_channel)^2);
    random_channel_gain = 10*log10(abs(random_channel)^2);
    est_CB_channel_gain = 10*log10(abs(CB_channel)^2);
    opt_CB_channel_gain = 10*log10(abs(opt_CB_channel)^2);
    opt_CB_bf_gain(sample) = opt_CB_channel_gain - EG_channel_gain;
    
    fprintf('Equal-gain beamforming channel gain: %f dB \n', EG_channel_gain)
    fprintf('Sensing-aided beam steering channel gain: %f dB \n', sensing_aided_channel_gain)
    fprintf('Sensing-aided codebook-based beamforming channel gain: %f dB \n', est_CB_channel_gain)
    fprintf('Optimal codebook-based beamforming channel gain: %f dB \n', opt_CB_channel_gain)
    fprintf('Random beamforming channel gain: %f dB \n', random_channel_gain)
    fprintf('\n')
    
    %%%%%%%%%%%%%%%%%% Metrics Calculation %%%%%%%%%%%%%%%%%%%%%%
    
    % Top-k beams versus beamforming gain
    best_top_k_channel_gain = est_CB_channel_gain;
    for i = 1:1:k
        CB_idx = sorted_CB_sim_idx(i);
        if CB_idx == opt_CB_idx
            for j = i:1:k
                top_k_acc(sample, j) = top_k_acc(sample, j) + 1;
            end
            break
        end
    end
    
    for i = 1:1:k
        CB_idx = sorted_CB_sim_idx(i);
        CB_RIS_interaction_vector = F_CB(:, CB_idx);
        CB_channel = double(APtoRIS_complex_channel.*RIStoUE_complex_channel).' * CB_RIS_interaction_vector;
        CB_channel_gain = 10*log10(abs(CB_channel)^2);
        if CB_channel_gain > best_top_k_channel_gain
            best_top_k_channel_gain = CB_channel_gain;
        end
        top_k_bf_gain(sample, i) = top_k_bf_gain(sample, i) + (best_top_k_channel_gain - EG_channel_gain);
    end

end

avg_opt_CB_bf_gain = 10 .* log10(mean(10.^(opt_CB_bf_gain./10)));
avg_top_k_acc = mean(top_k_acc, 1);
avg_top_k_bf_gain = 10 .* log10(mean(10.^(top_k_bf_gain./10), 1));

%% Save results
output_dir = "./results/";
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save( strcat(output_dir, 'avg_top_k_bf_gain_OSF', num2str(oversampling_x), 'x', num2str(oversampling_x), '.mat'), 'avg_top_k_bf_gain' )
save( strcat(output_dir, 'opt_CB_bf_gain_OSF', num2str(oversampling_x), 'x', num2str(oversampling_x), '.mat'), 'avg_opt_CB_bf_gain' )

%% External Functions
function antenna_map = antenna_channel_map(x, y, z, string)

    Mx_Ind=0:1:x-1;
    My_Ind=0:1:y-1;
    Mz_Ind=0:1:z-1;
    Mxx_Ind=repmat(Mx_Ind, 1, y*z)'; %col vector
    Myy_Ind=repmat(reshape(repmat(My_Ind,x,1), 1, x*y), 1, z)'; %col vector
    Mzz_Ind=reshape(repmat(Mz_Ind,x*y,1), 1, x*y*z)'; %col vector
    if string
        antenna_map = cellstr([num2str(Mxx_Ind) num2str(Myy_Ind) num2str(Mzz_Ind)]);
    else
        antenna_map = [Mxx_Ind, Myy_Ind, Mzz_Ind];
    end
end



