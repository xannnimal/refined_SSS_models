%% EXAMPLE: OPM sensors, multi-in, single-vsh out, iterative
%investigate iterative SSS method on multi-origin internal basis

%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
center1= [-0.00350699, 0.01138051, 0.05947857]; 
center2= [-0.00433911, 0.04081329, 0.05194245]; 
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];

%% opm geometry from Peter at SANDIA
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);

%% SSS expansions- multi origin interior
%speficy sensing direction. SQUID=R_hat or EZ, OPM=Theta or phi hat
sensing_dir=phi_hat;
other_dir=theta_hat;
%calculate multi-vsh in and single-vsh out
[SNin_tot,SNout] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',other_dir',sensing_dir',ch_types,Lin,Lout);

%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]'; %(size 3xN)
dipole_data = single_dipole_sim(opm_matrix,sensing_dir,dip_pos,dip_mom);
%pick a specific channel
phi_0= dipole_data.trial{1,1}(:,:);

%% iteratively reconstruct data
%   "S" is the full normalized SSS basis, in and out
%   "PHI" is the data (a vector or a matrix)
%   "nm" is the same as L_in -1
%   "nt" equals Lout - 1
%   "ni" equals number of iterations
ni=5;
data_rec = xi([SNin_tot,SNout],phi_0,Lin-1,Lout-1,ni);

%% plot data to check
%plot data from single channel
chan_num=1; %3 corresponds to MEG0111
data_time=dipole_data.time{1,1};
data_chan_num=dipole_data.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(data_time(:,1:100), data_chan_num(:,1:100))
plot(data_time(:,1:100),data_rec(chan_num,1:100))
title('Two-Origin SSS iterative, Channel 1, Sandia Helmet phi, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','Reconstructed'},'location','northwest')
hold off



