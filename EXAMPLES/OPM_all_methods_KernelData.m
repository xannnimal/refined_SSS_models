%% OPM with all recons- Kernal OPM Phantom data
%From Eric: audio_ERF_notebook_portal
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

%centers optimized using ACKLEE10 mri
% center1_new=[-0.00298921,-0.0212776,0.00999846];
% center2_new=[-0.00239802,0.0143463,0.0162393];
%after mri trans- not real trans from file, Eric made it
% center1_new=[-0.00425102,-0.0093374 ,-0.0198284 ];
% center2_new=[-0.00421421 ,0.02673303,-0.0171297];


%% UCL data from MNE-Python 
% coordsys = 'device'; 
% %%%fiff_file = "phantom_32_100nam_nobads_raw.fif";
% fname="phantom_OPM_evoked.fif";
% [opm_matrix,R_hat,theta_hat,phi_hat] = fiff_getpos(fname,coordsys);
% opm_matrix=opm_matrix';
% R_hat=R_hat';
% theta_hat=theta_hat';
% phi_hat=phi_hat';
% [evoked] = fiff_read_evoked_all(fname);
% phi=evoked.evoked.epochs;
% %remove bad channels, indicies adjusted from python
% %not 36,37,38,39..42,43..54,55]
% bad_inds=[37,38,39,40,43,44,55,56];
% k=1;
% for i=(1:size(phi,1))
%     if ismember(i,bad_inds)
%         k=k;
%     else
%         phi_0p(k,:)=phi(i,:);
%         k=k+1;
%     end
% end
% time=evoked.evoked.times;
% for i=(1:size(opm_matrix,1))
%     ch_types(i)=1; %model as magnetometers
% end

%% Kernel opm data: AKCLEE_110
coordsys = 'device'; 
filename= 'C:/Users/xanmc/RESEARCH/audio_ERF_notebook_portal/audio_ERF_portal_raw.fif';
%[R,EX,EY,EZ] = fiff_getpos(file,coord_system,calfile,z_offset)
[opm_matrix,EX,EY,EZ] = fiff_getpos(filename,coordsys);
opm_matrix=opm_matrix';
R_hat=EZ';
theta_hat=EX';
phi_hat=EY';

info = fiff_read_meas_info(filename);
[raw] = fiff_setup_read_raw(filename);
[data,times] = fiff_read_raw_segment(raw);
t_start=50001; %50sec
t_end=100001;
phi_0=data(:,50001:100001);
%check sensor layout and orientations
% figure(8)
% hold on
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),theta_hat(:,1),theta_hat(:,2),theta_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),phi_hat(:,1),phi_hat(:,2),phi_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),R_hat(:,1),R_hat(:,2),R_hat(:,3))
% grid on
% rotate3d
% view(135, 20);
% hold off
%visualize center positions
% figure(7);
% hold on
% scatter3(center1(1),center1(2),center1(3), 'r*')
% scatter3(center2(1),center2(2),center2(3), 'g*')
% scatter3(center1_new(1),center1_new(2),center1_new(3), 'b*')
% scatter3(center2_new(1),center2_new(2),center2_new(3), 'bl*')
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3))
% grid on
% rotate3d
% view(135, 20);
% hold off

for i=(1:size(opm_matrix,1))
    ch_types(i)=1; %model as magnetometers
end

%read evoked data
% [evoked] = fiff_read_evoked_all('ERF_Flux_meg_evoked.fif');
% nchan=evoked.info.nchan;
% time=evoked.evoked.times;
% phi_0p=evoked.evoked.epochs;

%read raw from matrix, fif functions not working
%t_start, t_end = 114.79969620704651, 564.6109187602997 from mne-python
% rawfile='ERF_Flux_meg_matrix.mat';
% raw = load(rawfile);
% raw_data=raw.data(:,30002:111111); %for 10sec (:,28802:30802);
% raw_times=raw.times(30002:111111); %for 10 sec (28802:30802);
% time=raw_times;
% phi_0p=raw_data;

%plot raw
figure(10);
hold on;
plot(times(:,50001:100001), phi_0)
title('Kernel OPM Auditory Raw Data')
xlabel('Time')
ylabel('(T)')
hold off

return

%% plot evoked
% figure(1);
% hold on;
% plot(time, phi_0p)
% %title('Kernel OPM Auidio Evoked Data')
% title('UCL OPM Auditory')
% xlabel('Time')
% ylabel('(T)')
% hold off


%% SSS expansions- phi
%speficy sensing direction. SQUID=R_hat or EZ, OPM=Theta or phi hat
%calculate single in single out
[Sin_p,SNin_p] = Sin_vsh_vv([0,0,0]',opm_matrix',theta_hat',phi_hat',R_hat',ch_types,Lin);
[Sout_p,SNout_p] = Sout_vsh_vv([0,0,0]',opm_matrix',theta_hat',phi_hat',R_hat',ch_types,Lout);

%calculate multi-vsh in and single-vsh out
[SNin_tot_p,SNout_p] = multiVSHin_singleVSHout(center1', center2',opm_matrix',theta_hat',phi_hat',R_hat',ch_types,Lin,Lout);

%calculate spheroidal in/out
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_matrix);
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(opm_matrix,R_hat,origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm_p,2)
  SNin_spm_p(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm_p(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end


%% reconstrct internal data
%single in, single out
pS_p=pinv([SNin_p SNout_p]);
XN_p=pS_p*phi_0p;
data_rec_vsh_p=real(SNin_p*XN_p(1:size(SNin_p,2),:));
%multi in, vsh out
pS_multi_vsh_p=pinv([SNin_tot_p SNout_p]);   
XN_multi_vsh_p=pS_multi_vsh_p*phi_0p;
data_rec_multi_vsh_p=real(SNin_tot_p*XN_multi_vsh_p(1:size(SNin_tot_p,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph_p=pinv([SNin_spm_p,SNout_spm_p]);  
XN_sph_sph_p=pS_sph_sph_p*phi_0p;
data_rec_sph_sph_p=real(SNin_spm_p*XN_sph_sph_p(1:size(SNin_spm_p,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_p=pinv([SNin_spm_p SNout_p]);   
XN_sph_vsh_p=pS_sph_vsh_p*phi_0p;
data_rec_sph_vsh_p=real(SNin_spm_p*XN_sph_vsh_p(1:size(SNin_spm_p,2),:));


%% check condition numbers
cond_vsh_vsh_p=cond([SNin_p SNout_p]);
cond_SNin_p=cond(SNin_p);
cond_SNin_tot_p = cond(SNin_tot_p);
cond_SNout_p= cond(SNout_p);
condition_multi_vsh_p = cond([SNin_tot_p SNout_p]);
cond_SNin_spm_p=cond(SNin_spm_p);
cond_SNout_spm_p= cond(SNout_spm_p);
condition_sph_sph_p = cond([SNin_spm_p SNout_spm_p]);
condition_sph_vsh_p = cond([SNin_spm_p SNout_p]);

%% compare sig values
% [U,sig,Vt]=svd(SNin_tot_p,"econ");
% sig=diag(sig);
% figure(10);
% hold on;
% plot(sig)
% title('singlular values multi-VSH in/out')
% hold off

%% plot data to check
ymin=-5e-11;
ymax=5e-11;

figure(2);
hold on;
plot(time, phi_0p(1,:))
plot(time,data_rec_vsh_p(1,:))
plot(time,data_rec_multi_vsh_p(1,:))
plot(time,data_rec_sph_sph_p(1,:))
plot(time,data_rec_sph_vsh_p(1,:))
title('All SSS Methods, Kernel OPM Auidio Evoked, R-hat Sensing')
xlabel('Time')
ylabel('Mag 1 (T)')
ylim([ymin, ymax])
%legend({'Raw Data','VSH/VSH','Spm/Spm'},'location','northwest')
legend({'Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off

% Create block plots
figure(6);
hold on
t = tiledlayout(2,2); %3
ax1 = nexttile;
plot(ax1,time,data_rec_vsh_p)
ylim(ax1,[ymin, ymax])
title(ax1,'sVSH/sVSH')
ax2=nexttile;
plot(ax2,time,data_rec_multi_vsh_p)
ylim(ax2,[ymin, ymax])
title(ax2,'mVSH/sVSH')
ax3=nexttile;
plot(ax3,time,data_rec_sph_sph_p)
ylim(ax3,[ymin, ymax])
title(ax3,'Spheroid in/out')
ax4=nexttile;
plot(ax4,time,data_rec_sph_vsh_p)
ylim(ax4,[ymin, ymax])
title(ax4, 'Spheroid/sVSH')
% ax5=nexttile;
% plot(ax5,time,phi_0p)
% ylim(ax5,[ymin, ymax])
% title(ax5, 'Raw Data')
title(t, 'Kernel OPM Auidio Evoked, R-hat Sensing')
xlabel(t,'Time')
ylabel(t,'T')
% Move plots closer together
t.TileSpacing = 'compact';
hold off

return
%% SNR calculations (signal/noise)
%noise from time -0.2 to -0.05, indicies (1,31), 30
%peak from time 0.05 to 0.15, indicies (51,71), 20
%If A is a matrix whose columns are random variables and whose rows are observations 
% then S is a row vector containing the standard deviation corresponding to each column
std_noise_raw=mean(std(phi_0p(:,1:31))); %each entry is the std of all the channels at 1 time, then average over all time
std_peak_raw=mean(std(phi_0p(:,51:71)));

SNR_raw = std_peak_raw/std_noise_raw;
SNR_vsh_p = mean(std(data_rec_vsh_p(:,51:71)))/mean(std(data_rec_vsh_p(:,1:31)));
SNR_mvsh_p = mean(std(data_rec_multi_vsh_p(:,51:71)))/mean(std(data_rec_multi_vsh_p(:,1:31)));
SNR_sphsph_p = mean(std(data_rec_sph_sph_p(:,51:71)))/mean(std(data_rec_sph_sph_p(:,1:31)));
SNR_sphvsh_p = mean(std(data_rec_sph_vsh_p(:,51:71)))/mean(std(data_rec_sph_vsh_p(:,1:31)));


%% compare data reconstructions
sVSH_sVSH_p=[SNin_p SNout_p];
mVSH_sVSH_p=[SNin_tot_p SNout_p];
oid_oid_p=[Sin_spm_p,Sout_spm_p];
oid_sVSH_p=[Sin_spm_p,SNout_p];

for i=(1:size(time,2))
    check_data_vsh_vsh_p(i) = subspace(phi_0p(:,i), sVSH_sVSH_p)*180/pi;
    check_data_mvsh_vsh_p(i) = subspace(phi_0p(:,i), mVSH_sVSH_p)*180/pi;
    check_data_oid_oid_p(i) = subspace(phi_0p(:,i), oid_oid_p)*180/pi;
    check_data_oid_vsh_p(i) = subspace(phi_0p(:,i), oid_sVSH_p)*180/pi;
end

check_data_vsh_vsh_pmin = min(check_data_vsh_vsh_p);
check_data_vsh_vsh_pmax = max(check_data_vsh_vsh_p);
check_data_vsh_vsh_pav = mean(check_data_vsh_vsh_p);

check_data_mvsh_vsh_pmin = min(check_data_mvsh_vsh_p);
check_data_mvsh_vsh_pmax = max(check_data_mvsh_vsh_p);
check_data_mvsh_vsh_pav = mean(check_data_mvsh_vsh_p);

check_data_oid_oid_pmin = min(check_data_oid_oid_p);
check_data_oid_oid_pmax = max(check_data_oid_oid_p);
check_data_oid_oid_pav = mean(check_data_oid_oid_p);

check_data_oid_vsh_pmin = min(check_data_oid_vsh_p);
check_data_oid_vsh_pmax = max(check_data_oid_vsh_p);
check_data_oid_vsh_pac = mean(check_data_oid_vsh_p);

%% Compare with Iterative SSS methods %%





