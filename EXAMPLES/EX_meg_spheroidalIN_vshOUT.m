%% EXAMPLE: SQUID sensors, spheroid-in, single-vsh out
clear
%% constant variables regardless
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";
[R_mag,EX_mag,EY_mag,EZ_mag,ch_types] = gen_squid_geometry(rawfile, coordsys);

%% SSS expansions
%calculate spheroidal in and single-vsh out
[semi_major,semi_minor]=find_ellipse_axis(R_mag');
[SNin_spm,SNout] = spheroidIN_vshOUT(R_mag,EX_mag,EY_mag,EZ_mag,semi_major,semi_minor,Lin,Lout,ch_types);

%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]'; %(size 3xN)
dipole_data = single_dipole_sim(opm_matrix,sensing_dir,dip_pos,dip_mom);
%pick a specific channel
phi_0= dipole_data.trial{1,1}(:,:);

%% reconstrct internal data
pS=pinv([SNin_spm SNout]);   
XN=pS*phi_0;
data_rec=real(SNin_spm*XN(1:size(SNin_spm,2),:)); 
%check condition numbers
condition_in = cond(SNin_spm);
condition_out= cond(SNout);
condition_both = cond([SNin_spm SNout]);

%% plot data to check
%plot data from single channel
chan_num=1; %3 corresponds to MEG0111
data_time=dipole_data.time{1,1};
data_chan_num=dipole_data.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(data_time(:,1:100), data_chan_num(:,1:100))
plot(data_time(:,1:100),data_rec(chan_num,1:100))
title('spheroidal and vsh, Channel 1, 306 Chan, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','Reconstructed'},'location','northwest')
hold off

