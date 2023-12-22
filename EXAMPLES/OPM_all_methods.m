%% OPM with all recons
clear
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
%find semi major and minor
%[semi_major,semi_minor]=find_ellipse_axis(opm_matrix);
for i=(1:144)
    r(i)=sqrt(opm_matrix(i,1)^2+ opm_matrix(i,2)^2+ opm_matrix(i,3)^2);
end
semi_major=max(r);
semi_minor=min(r);
semi_major=0.11;
semi_minor=0.09;
%calculate multi-vsh in and single-vsh out
[SNin_tot,SNout] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',other_dir',sensing_dir',ch_types,Lin,Lout);
[Sin_spm,Sout_spm] = spheroidIN_spheroidOUT(opm_matrix,R_hat,other_dir,sensing_dir,semi_major,semi_minor,Lin,Lout);

%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]'; %(size 3xN)
dipole_data = single_dipole_sim(opm_matrix,sensing_dir,dip_pos,dip_mom);
%pick a specific channel
phi_0= dipole_data.trial{1,1}(:,:);

%% reconstrct internal data
%multi in, vsh out
pS_multi_vsh=pinv([SNin_tot SNout]);   
XN_multi_vsh=pS_multi_vsh*phi_0;
data_rec_multi_vsh=real(SNin_tot*XN_multi_vsh(1:size(SNin_tot,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph=pinv([Sin_spm Sout_spm]);   
XN_sph_sph=pS_sph_sph*phi_0;
data_rec_sph_sph=real(Sin_spm*XN_sph_sph(1:size(Sin_spm,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh=pinv([Sin_spm SNout]);   
XN_sph_vsh=pS_sph_vsh*phi_0;
data_rec_sph_vsh=real(Sin_spm*XN_sph_vsh(1:size(Sin_spm,2),:));

%check condition numbers
cond_SNin_tot = cond(SNin_tot);
cond_SNout= cond(SNout);
condition_multi_vsh = cond([SNin_tot SNout]);
cond_SNin_spm=cond(Sin_spm);
cond_SNout_spm= cond(Sout_spm);
condition_sph_sph = cond([Sin_spm Sout_spm]);
condition_sph_vsh = cond([Sin_spm SNout]);

%% plot data to check
%plot data from single channel
chan_num=1; %3 corresponds to MEG0111
data_time=dipole_data.time{1,1};
data_chan_num=dipole_data.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(data_time(:,1:100), data_chan_num(:,1:100))
plot(data_time(:,1:100),data_rec_multi_vsh(chan_num,1:100))
plot(data_time(:,1:100),data_rec_sph_sph(chan_num,1:100))
plot(data_time(:,1:100),data_rec_sph_vsh(chan_num,1:100))
title('Two-Origin SSS, Channel 1, Sandia Helmet phi, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off

