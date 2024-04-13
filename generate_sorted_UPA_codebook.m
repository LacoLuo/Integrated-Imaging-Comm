clearvars
clc

%% Define output directory
output_dir = "./codebooks/";
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Generate codebook
Mx = 40;
My = 1;
Mz = 40;

oversampling_x = 4;
oversampling_y = 1;
oversampling_z = 4;

codebook_size_x = Mx * oversampling_x;
codebook_size_y = My * oversampling_y;
codebook_size_z = Mz * oversampling_z;

ant_spacing = .5;
kd = 2 * pi * ant_spacing;
antx_index = 0:1:Mx-1;
anty_index = 0:1:My-1;
antz_index = 0:1:Mz-1;
M = Mx * My * Mz;

theta_qx = 0:pi/codebook_size_x:pi-1e-6; % quantized beamsteering angles
F_CBx = zeros(Mx, codebook_size_x);
for i=1:1:length(theta_qx)
    F_CBx(:,i) = exp( +1j * kd * antx_index' * cos(theta_qx(i)) );
end

theta_x = zeros(codebook_size_x, 1);
for i = 1:1:codebook_size_x
    theta_x(i) = estimate_angles(F_CBx(:, i)) / (pi/180);
end

[sorted_theta_x, sorted_theta_x_idx] = sort(theta_x, 'descend');
sorted_F_CBx = F_CBx(:, sorted_theta_x_idx);
clf
for i = 1:1:codebook_size_x
    tmp = estimate_angles(sorted_F_CBx(:, i)) / (pi/180);
end

theta_qy=0:pi/codebook_size_y:pi-1e-6; % quantized beamsteering angles
F_CBy=zeros(My,codebook_size_y);
for i=1:1:length(theta_qy)
    F_CBy(:,i)=exp(+1j*kd*anty_index'*cos(theta_qy(i)));
end
 
theta_qz=0:pi/codebook_size_z:pi-1e-6; % quantized beamsteering angles
F_CBz=zeros(Mz,codebook_size_z);
for i=1:1:length(theta_qz)
    F_CBz(:,i)=exp(+1j*kd*antz_index'*cos(theta_qz(i)));
end

clf
theta_z = zeros(codebook_size_z, 1);
for i = 1:1:codebook_size_z
    theta_z(i) = estimate_angles(F_CBz(:, i)) / (pi/180);
end

[sorted_theta_z, sorted_theta_z_idx] = sort(theta_z, 'ascend');
sorted_F_CBz = F_CBz(:, sorted_theta_z_idx);

F_CBxy=kron(F_CBy, sorted_F_CBx);
F_CB=kron(sorted_F_CBz, F_CBxy);

% Save codebook
output_filename = strcat(output_dir, "UPA_codebook_", num2str(Mx), "x", num2str(Mz), "_OSF_", ...
                                        num2str(oversampling_x), "x", num2str(oversampling_z), ".mat");
save(output_filename, "F_CB")

%%  External Functions 
function theta = estimate_angles(vec)

    My = size(vec, 1);
    over_sampling_y = 1000;
    [F,~] = UPA_codebook_generator(1,My,1,1,over_sampling_y,1,.5); %F: (#ant, #sampled_directions)
    theta_s = 0:pi/(over_sampling_y*My):pi-1e-6;
    projection = ctranspose(F)*vec;
    proj = abs(projection).^2;
    [~, argidx] = max(proj);
    theta = theta_s(argidx);

    figure(1);
    fig_position = get(gcf, 'Position');
    set(gcf, 'Position', [10, 10, fig_position(3), fig_position(4)]);
    for n=1:1:size(vec, 2) % #beams
        polarplot(theta_s, proj(:,n).')
        hold on;
    end
    grid on
    box on
    hold on

end

function [F_CB,all_beams]=UPA_codebook_generator(Mx,My,Mz,over_sampling_x,over_sampling_y,over_sampling_z,ant_spacing)

    kd=2*pi*ant_spacing;
    antx_index=0:1:Mx-1;
    anty_index=0:1:My-1;
    antz_index=0:1:Mz-1;
    M=Mx*My*Mz;
    
    % Defining the RF beamforming codebook in the x-direction
    codebook_size_x=over_sampling_x*Mx;
    codebook_size_y=over_sampling_y*My;
    codebook_size_z=over_sampling_z*Mz;
    
    
    theta_qx=0:pi/codebook_size_x:pi-1e-6; % quantized beamsteering angles
    F_CBx=zeros(Mx,codebook_size_x);
    for i=1:1:length(theta_qx)
        F_CBx(:,i)=sqrt(1/Mx)*exp(+1j*kd*antx_index'*cos(theta_qx(i)));
    end
     
    theta_qy=0:pi/codebook_size_y:pi-1e-6; % quantized beamsteering angles
    F_CBy=zeros(My,codebook_size_y);
    for i=1:1:length(theta_qy)
       F_CBy(:,i)=sqrt(1/My)*exp(+1j*kd*anty_index'*cos(theta_qy(i)));
    end
     
    theta_qz=0:pi/codebook_size_z:pi-1e-6; % quantized beamsteering angles
    F_CBz=zeros(Mz,codebook_size_z);
    for i=1:1:length(theta_qz)
        F_CBz(:,i)=sqrt(1/Mz)*exp(+1j*kd*antz_index'*cos(theta_qz(i)));
    end
    
    F_CBxy=kron(F_CBy,F_CBx);
    F_CB=kron(F_CBz,F_CBxy);
    
    beams_x=1:1:codebook_size_x;
    beams_y=1:1:codebook_size_y;
    beams_z=1:1:codebook_size_z;
    
    Mxx_Ind=repmat(beams_x,1,codebook_size_y*codebook_size_z)';
    Myy_Ind=repmat(reshape(repmat(beams_y,codebook_size_x,1),1,codebook_size_x*codebook_size_y),1,codebook_size_z)';
    Mzz_Ind=reshape(repmat(beams_z,codebook_size_x*codebook_size_y,1),1,codebook_size_x*codebook_size_y*codebook_size_z)';
    
    Tx=cat(3,Mxx_Ind',Myy_Ind',Mzz_Ind');
    all_beams=reshape(Tx,[],3);
end