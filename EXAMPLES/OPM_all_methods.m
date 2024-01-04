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

%% SSS expansions- multi origin interior- phi
%speficy sensing direction. SQUID=R_hat or EZ, OPM=Theta or phi hat
%find semi major and minor
%[semi_major,semi_minor]=find_ellipse_axis(opm_matrix);
for i=(1:size(opm_matrix,1))
    r(i)=sqrt(opm_matrix(i,1)^2+ opm_matrix(i,2)^2+ opm_matrix(i,3)^2);
end
semi_major=max(r);
semi_minor=min(r);
semi_major=0.11;
semi_minor=0.095;
%calculate multi-vsh in and single-vsh out
[SNin_tot_p,SNout_p] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin,Lout);
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(opm_matrix,R_hat,theta_hat,phi_hat,semi_major,semi_minor,Lin,Lout);

%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]'; %(size 3xN)
dipole_data_p = single_dipole_sim(opm_matrix,phi_hat,dip_pos,dip_mom);
%pick a specific channel
phi_0p= dipole_data_p.trial{1,1}(:,:);
%calculate B field
%magneticField = magneticDipole(dip_pos,dip_mom);

%% reconstrct internal data
%single in, single out
[Sin_p,SNin_p] = Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin);
[Sout_p,SNout_p] = Sout_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lout);
pS_p=pinv([SNin_p SNout_p]);
XN_p=pS_p*phi_0p;
data_rec_vsh_p=real(SNin_p*XN_p(1:size(SNin_p,2),:));
%multi in, vsh out
pS_multi_vsh_p=pinv([SNin_tot_p SNout_p]);   
XN_multi_vsh_p=pS_multi_vsh_p*phi_0p;
data_rec_multi_vsh_p=real(SNin_tot_p*XN_multi_vsh_p(1:size(SNin_tot_p,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph_p=pinv([Sin_spm_p Sout_spm_p]);   
XN_sph_sph_p=pS_sph_sph_p*phi_0p;
data_rec_sph_sph_p=real(Sin_spm_p*XN_sph_sph_p(1:size(Sin_spm_p,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_p=pinv([Sin_spm_p SNout_p]);   
XN_sph_vsh_p=pS_sph_vsh_p*phi_0p;
data_rec_sph_vsh_p=real(Sin_spm_p*XN_sph_vsh_p(1:size(Sin_spm_p,2),:));

%check condition numbers
cond_vsh_vsh_p=cond([SNin_p SNout_p]);
cond_SNin_p=cond(SNin_p);
cond_SNin_tot_p = cond(SNin_tot_p);
cond_SNout_p= cond(SNout_p);
condition_multi_vsh_p = cond([SNin_tot_p SNout_p]);
cond_SNin_spm_p=cond(Sin_spm_p);
cond_SNout_spm_p= cond(Sout_spm_p);
condition_sph_sph_p = cond([Sin_spm_p Sout_spm_p]);
condition_sph_vsh_p = cond([Sin_spm_p SNout_p]);

%calculate the subspace angle between the reconstructed and noiseless original data for one time instant
time=10;
sub_multi_vsh_p=subspace(phi_0p(:,time),data_rec_multi_vsh_p(:,time))*180/pi;
sub_sph_sph_p=subspace(phi_0p(:,time),data_rec_sph_sph_p(:,time))*180/pi;
sub_sph_vsh_p=subspace(phi_0p(:,time),data_rec_sph_vsh_p(:,time))*180/pi;
sub_vsh_vsh_p=subspace(phi_0p(:,time),data_rec_vsh_p(:,time))*180/pi;



%% plot data to check
%plot data from single channel
chan_num=1; %3 corresponds to MEG0111
data_time_p=dipole_data_p.time{1,1};
data_chan_num_p=dipole_data_p.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(data_time_p(:,1:100), data_chan_num_p(:,1:100))
plot(data_time_p(:,1:100),data_rec_vsh_p(chan_num,1:100))
plot(data_time_p(:,1:100),data_rec_multi_vsh_p(chan_num,1:100))
plot(data_time_p(:,1:100),data_rec_sph_sph_p(chan_num,1:100))
plot(data_time_p(:,1:100),data_rec_sph_vsh_p(chan_num,1:100))
title('All SSS Methods, Channel 1, Sandia Helmet Phi, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off

%%%%%%%%%%%%%%%
%% compare with theta as sensing dir, mVSH
%calculate multi-vsh in and single-vsh out
[SNin_tot_t,SNout_t] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',phi_hat',theta_hat',ch_types,Lin,Lout);
[Sin_spm_t,Sout_spm_t] = spheroidIN_spheroidOUT(opm_matrix,R_hat,phi_hat,theta_hat,semi_major,semi_minor,Lin,Lout);

%% generate single dipole simulated data
dip_pos = [0.01,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,0,1]'; %(size 3xN)
dipole_data_t = single_dipole_sim(opm_matrix,theta_hat,dip_pos,dip_mom);
%pick a specific channel
phi_0t= dipole_data_t.trial{1,1}(:,:);
%calculate B field
%magneticField = magneticDipole(dip_pos,dip_mom);

%% reconstrct internal data
%single in, single out
[Sin_t,SNin_t] = Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',phi_hat',theta_hat',ch_types,Lin);
[Sout_t,SNout_t] = Sout_vsh_vv([0,0,0]',opm_matrix',R_hat',phi_hat',theta_hat',ch_types,Lout);
pS_t=pinv([SNin_t SNout_t]);
XN_t=pS_t*phi_0t;
data_rec_vsh_t=real(SNin_t*XN_t(1:size(SNin_t,2),:));
%multi in, vsh out
pS_multi_vsh_t=pinv([SNin_tot_t SNout_t]);   
XN_multi_vsh_t=pS_multi_vsh_t*phi_0t;
data_rec_multi_vsh_t=real(SNin_tot_t*XN_multi_vsh_t(1:size(SNin_tot_t,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph_t=pinv([Sin_spm_t Sout_spm_t]);   
XN_sph_sph_t=pS_sph_sph_t*phi_0t;
data_rec_sph_sph_t=real(Sin_spm_t*XN_sph_sph_t(1:size(Sin_spm_t,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_t=pinv([Sin_spm_t SNout_t]);   
XN_sph_vsh_t=pS_sph_vsh_t*phi_0t;
data_rec_sph_vsh_t=real(Sin_spm_t*XN_sph_vsh_t(1:size(Sin_spm_t,2),:));

%check condition numbers
cond_vsh_vsh_t=cond([SNin_t SNout_t]);
cond_SNin_t=cond(SNin_t);
cond_SNin_tot_t = cond(SNin_tot_t);
cond_SNout_t= cond(SNout_t);
condition_multi_vsh_t = cond([SNin_tot_t SNout_t]);
cond_SNin_spm_t=cond(Sin_spm_t);
cond_SNout_spm_t= cond(Sout_spm_t);
condition_sph_sph_t = cond([Sin_spm_t Sout_spm_t]);
condition_sph_vsh_t = cond([Sin_spm_t SNout_t]);

%calculate the subspace angle between the reconstructed and noiseless original data for one time instant
time=10;
sub_multi_vsh_t=subspace(phi_0t(:,time),data_rec_multi_vsh_t(:,time))*180/pi;
sub_sph_sph_t=subspace(phi_0t(:,time),data_rec_sph_sph_t(:,time))*180/pi;
sub_sph_vsh_t=subspace(phi_0t(:,time),data_rec_sph_vsh_t(:,time))*180/pi;
sub_vsh_vsh_t=subspace(phi_0t(:,time),data_rec_vsh_t(:,time))*180/pi;



%% plot data to check
%plot data from single channel
chan_num=1; %3 corresponds to MEG0111
data_time_t=dipole_data_t.time{1,1};
data_chan_num_t=dipole_data_t.trial{1,1}(chan_num,:); 

figure(3);
hold on;
plot(data_time_t(:,1:100), data_chan_num_t(:,1:100))
plot(data_time_t(:,1:100),data_rec_vsh_t(chan_num,1:100))
plot(data_time_t(:,1:100),data_rec_multi_vsh_t(chan_num,1:100))
plot(data_time_t(:,1:100),data_rec_sph_sph_t(chan_num,1:100))
plot(data_time_t(:,1:100),data_rec_sph_vsh_t(chan_num,1:100))
title('All SSS Methods, Channel 1, Sandia Helmet Theta, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off

%% compare sig values
[U_t,sig_t,Vt_t]=svd(SNin_tot_t,"econ");
[U_p,sig_p,Vt_p]=svd(SNin_tot_p,"econ");
sig_t=diag(sig_t);
sig_p=diag(sig_p);
figure(4);
hold on;
plot(sig_t)
plot(sig_p)
legend({'theta','phi'},'location','northwest')
hold off

