%% EXAMPLE: OPM sensors, spheroid-in, spheroid-out
clear
%% constant variables regardless
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% opm geometry from Peter at SANDIA
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);

%% SSS expansions- multi origin interior
%speficy sensing direction. SQUID=R_hat or EZ, OPM=Theta or phi hat
sensing_dir=phi_hat;
other_dir=theta_hat;
%find major and minor axis of spheroidal ellipse
%Y must be the longest axis of opm_matrix
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_matrix);
%calculate spheroidal in and single-vsh out
%matrix,chan_ori,semi_major,semi_minor,Lin,Lout
[Sin_spm,Sout_spm] = spheroidIN_spheroidOUT(opm_matrix,phi_hat,origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm,2)
  SNin_spm(:,j) = Sin_spm(:,j)/norm(Sin_spm(:,j));
end

for j = 1:size(Sout_spm,2)
  SNout_spm(:,j) = Sout_spm(:,j)/norm(Sout_spm(:,j));
end

%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]; %(size 3xN)
phi_0 = magneticDipole_pointMags(opm_matrix',sensing_dir',dip_pos', dip_mom')';
%dipole_data = single_dipole_sim(opm_matrix,sensing_dir,dip_pos,dip_mom);
%pick a specific channel
%phi_0= dipole_data.trial{1,1}(:,:);

%% reconstrct internal data
pS=pinv([SNin_spm SNout_spm]);   
XN=pS*phi_0;
data_rec=real(SNin_spm*XN(1:size(SNin_spm,2),:)); 
%check condition numbers
condition_in = cond(SNin_spm);
condition_out= cond(SNout_spm);
condition_both = cond([SNin_spm SNout_spm]);

angle = subspace(phi_0,SNin_spm);

%% plot data to check
%plot data from single channel
% chan_num=1; %3 corresponds to MEG0111
% data_time=dipole_data.time{1,1};
% data_chan_num=dipole_data.trial{1,1}(chan_num,:); 
% 
% figure(2);
% hold on;
% plot(data_time(:,1:100), data_chan_num(:,1:100))
% plot(data_time(:,1:100),data_rec(chan_num,1:100))
% title('spheroidal and vsh, Channel 1, Sandia Helmet phi, dipole 1cm x')
% xlabel('time')
% ylabel('MEG0121')
% %ylim([-8e-12 8e-12])
% legend({'Raw Data','Reconstructed'},'location','northwest')
% hold off

