%% EXAMPLE: SQUID sensors, multi-origin VSH in, single VSH out

%% constant variables
import mne.*
magscale = 100; % Numerical scaling factor between magnetometer and gradiometer signals
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
thresh = 0.98;
%use optimized centers for two spheres from Eric's code sphere_opt.py
center1=[-0.00350699, 0.01138051, 0.05947857];
center2=[-0.00433911, 0.04081329, 0.05194245];
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];

%% generate SQUID magnetometers
coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";
[R_mag,EX_mag,EY_mag,EZ_mag,ch_types] = gen_squid_geometry(rawfile, coordsys);

%% SSS expansions
sensing_dir=EZ_mag;
other_dir=EY_mag;
[SNin_tot,SNout] = multiVSHin_singleVSHout(center1', center2',R_mag,EX_mag,other_dir,sensing_dir,ch_types,Lin,Lout);


%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]'; %(size 3xN)
dipole_data = single_dipole_sim(R_mag',sensing_dir',dip_pos,dip_mom);
%pick a specific channel
phi_0= dipole_data.trial{1,1}(:,:);

%% reconstrct internal data
pS=pinv([SNin_tot SNout]);   
XN=pS*phi_0;
data_rec=real(SNin_tot*XN(1:size(SNin_tot,2),:)); 
%check condition numbers
condition_in = cond(SNin_tot);
condition_out= cond(SNout);
condition_both = cond([SNin_tot SNout]);

%% plot data to check
%plot data from single channel
chan_num=1; %3 corresponds to MEG0111
data_time=dipole_data.time{1,1};
data_chan_num=dipole_data.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(data_time(:,1:100), data_chan_num(:,1:100))
plot(data_time(:,1:100),data_rec(chan_num,1:100))
title('Two-Origin SSS, Channel 1, 102chan SQUID MEG, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','Reconstructed'},'location','northwest')
hold off

