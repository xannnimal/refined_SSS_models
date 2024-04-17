%% Investigating other ways to combine two interior bases together
% edited to investigate methods of combining two spheres after ILABS
% meeting beginning of March
%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
center1= [-0.00350699, 0.01138051, 0.05947857]; 
center2= [-0.00433911, 0.04081329, 0.05194245]; 
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];

%% generate SQUID magnetometers
coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";
%for 306 channels, run this
[R,EX,EY,EZ] = fiff_getpos(rawfile,coordsys);
for i=(1:size(EX,2))
    if mod(i,3)==0 %every third is a magnetometer
        ch_types(i)=1;
    else
        ch_types(i)=0;
    end
end
k=1;
for i=(1:306)
    if ch_types(i)==1 %every third is a magnetometer
        mags(k)=i;
        k=k+1;
    else
        k=k;
    end
end

%% opm geometry from Peter at SANDIA
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types_opm] = gen_opm_geometry(filename);
nchan = size(ch_types_opm,1);
mags_opm = 1:1:nchan;

%% SSS expansions- SQUID
%single VSH
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
%Multi-VSH
[SNin_tot,SNout] = multiVSHin_singleVSHout(center1', center2',R,EX,EY,EZ,ch_types,Lin,Lout);

%Investigating other methods for combinging two interior expansions
[~,SNin_1] = Sin_vsh_vv(center1',R,EX,EY,EZ,ch_types,Lin);
[~,SNin_2] = Sin_vsh_vv(center2',R,EX,EY,EZ,ch_types,Lin);
[SNin_tot_svd,sig,~]=svd([SNin_1,SNin_2],'econ');
SNin_tot_svd=SNin_tot_svd(:,1:80); %keep first 80 to make it the same size as SNin,SNout
SNin_tot_orth=orth(SNin_tot);
SNin_tot_orth=SNin_tot_orth(:,1:80);



%% SSS expansions - SANDIA OPM - Phi hat sensing
%single in, single out
[Sin_p,SNin_p] = Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types_opm,Lin);
[Sout_p,SNout_p] = Sout_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types_opm,Lout);

%calculate multi-vsh in and single-vsh out
[SNin_tot_p,SNout_p] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',theta_hat',phi_hat',ch_types_opm,Lin,Lout);
%other methods for mVSH: use svd, first 80 singular values
[~,SNin_1p] = Sin_vsh_vv(center1',opm_matrix',R_hat',theta_hat',phi_hat',ch_types_opm,Lin);
[~,SNin_2p] = Sin_vsh_vv(center2',opm_matrix',R_hat',theta_hat',phi_hat',ch_types_opm,Lin);
[SNin_tot_svd_p,~,~]=svd([SNin_1p,SNin_2p],'econ');
SNin_tot_svd_p=SNin_tot_svd_p(:,1:80); 
%other method: orth of original SNin_tot, keep first 80
SNin_tot_orth_p=orth(SNin_tot_p);
SNin_tot_orth_p=SNin_tot_orth_p(:,1:80);

%% SSS expansions - SANDIA OPM - Theta hat sensing
%single in, single out
[Sin_t,SNin_t] = Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',phi_hat',theta_hat',ch_types_opm,Lin);
[Sout_t,SNout_t] = Sout_vsh_vv([0,0,0]',opm_matrix',R_hat',phi_hat',theta_hat',ch_types_opm,Lout);

%calculate multi-vsh in and single-vsh out
[SNin_tot_t,SNout_t] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',phi_hat',theta_hat',ch_types_opm,Lin,Lout);
%other methods for mVSH: use svd, first 80 singular values
[~,SNin_1t] = Sin_vsh_vv(center1',opm_matrix',R_hat',phi_hat',theta_hat',ch_types_opm,Lin);
[~,SNin_2t] = Sin_vsh_vv(center2',opm_matrix',R_hat',phi_hat',theta_hat',ch_types_opm,Lin);
[SNin_tot_svd_t,~,~]=svd([SNin_1p,SNin_2t],'econ');
SNin_tot_svd_t=SNin_tot_svd_t(:,1:80); 
%other method: orth of original SNin_tot, keep first 80
SNin_tot_orth_t=orth(SNin_tot_t);
SNin_tot_orth_t=SNin_tot_orth_t(:,1:80);

%% check condition numbers
condition_sVSH = cond(SNin);
condition_mVSH = cond(SNin_tot);
condition_in_orth = cond(SNin_tot_orth);
condition_in_svd = cond(SNin_tot_svd);

condition_sVSH_p = cond(SNin_p);
condition_mVSH_p = cond(SNin_tot_p);
condition_in_orth_p = cond(SNin_tot_orth_p);
condition_in_svd_p = cond(SNin_tot_svd_p);

condition_sVSH_t = cond(SNin_t);
condition_mVSH_t = cond(SNin_tot_t);
condition_in_orth_t = cond(SNin_tot_orth_t);
condition_in_svd_t = cond(SNin_tot_svd_t);


% condition_both = cond([SNin SNout]);
% condition_both_m = cond([SNin_tot SNout]);
% condition_both_orth = cond([SNin_tot_orth SNout]);
% condition_both_svd = cond([SNin_tot_svd SNout]);


%% simulate current dipole
%current dipole using Samu's implementation of Sarvas
dip_mom_out=[1,0,0];
dip_pos_out = [0,5,20];
rs=[0,0,0];
q=[0,1,0]; %y direction
r0=[0.05,0,0]; %5cm along x axis
f_start = 100; % start frequency
f_end = 50; % end frequency
f_start_out = 50; % start frequency
f_end_out = 30; % end frequency
timestep = 0.0001;
T = 0.05;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;
% 
for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2));
end
% 
% %current dipole in, magnetic dipole out
for i=(1:size(times,2))
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t_out(:,i),ch_types)'*(1e13);
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)';

    phi_out_p(:,i) = magneticDipole(opm_matrix',R_hat',theta_hat',phi_hat',dip_pos_out', dip_mom_t_out(:,i),ch_types)'*1e13;
    phi_in_p(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',opm_matrix',R_hat',theta_hat',phi_hat',mags)';

    phi_out_t(:,i) = magneticDipole(opm_matrix',R_hat',phi_hat',theta_hat',dip_pos_out', dip_mom_t_out(:,i),ch_types)'*1e13;
    phi_in_t(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',opm_matrix',R_hat',phi_hat',theta_hat',mags)';
end
phi_0=phi_in +phi_out;
phi_0p = phi_in_p+phi_out_p;
phi_0t = phi_in_t+phi_out_t;


% phi_0 = dipole_field_sarvas(rs',q',r0',R,EX,EY,EZ,mags)';
% phi_0p = dipole_field_sarvas(rs',q',r0',opm_matrix',R_hat',theta_hat',phi_hat',mags_opm)';
% phi_0t = dipole_field_sarvas(rs',q',r0',opm_matrix',R_hat',phi_hat',theta_hat',mags_opm)';


%% Kernel opm data: AKCLEE_110 updated April 2024 correct sensor positions
coordsys = 'device'; 
filename= 'C:/Users/xanmc/RESEARCH/audio_ERF_notebook_portal/audio_ERF_portal_raw.fif';
info = fiff_read_meas_info(filename);
nchan_k=info.nchan;
ch_types_k=ones(nchan_k,1);
% [raw] = fiff_setup_read_raw(filename);
% [data,times] = fiff_read_raw_segment(raw);
% t_start=206022; %200sec
% t_end=260420; %260
% phi_raw=data(:,t_start:t_end);
for i=1:nchan_k
    opm_matrix_k(:,i)=info.chs(i).loc(1:3,:);
    theta_hat_k(:,i)=info.chs(i).loc(4:6,:);
    phi_hat_k(:,i)=info.chs(i).loc(7:9,:);
    R_hat_k(:,i)=info.chs(i).loc(10:12,:);
end
%read evoked
file= 'C:/Users/xanmc/RESEARCH/audio_ERF_notebook_portal/audio_ERF_portal_evoked.fif';
[evoked] = fiff_read_evoked(file);
evoked_data=evoked.evoked.epochs;
phi_0k=evoked_data;
evoked_times = evoked.evoked.times;
time=evoked.evoked.times;


%calculate single in single out
[Sin_k,SNin_k] = Sin_vsh_vv([0,0,0]',opm_matrix_k,theta_hat_k,phi_hat_k,R_hat_k,ch_types_k,Lin);
[Sout_k,SNout_k] = Sout_vsh_vv([0,0,0]',opm_matrix_k,theta_hat_k,phi_hat_k,R_hat_k,ch_types_k,Lout);
%calculate multi-vsh in and single-vsh out
[SNin_tot_k,SNout_k] = multiVSHin_singleVSHout(center1', center2',opm_matrix_k,theta_hat_k,phi_hat_k,R_hat_k,ch_types_k,Lin,Lout);
%other methods for mVSH: use svd, first 80 singular values
[~,SNin_1k] = Sin_vsh_vv(center1',opm_matrix_k,R_hat_k,theta_hat_k,phi_hat_k,ch_types_k,Lin);
[~,SNin_2k] = Sin_vsh_vv(center2',opm_matrix_k,R_hat_k,theta_hat_k,phi_hat_k,ch_types_k,Lin);
[SNin_tot_svd_k,~,~]=svd([SNin_1k,SNin_2k],'econ');
SNin_tot_svd_k=SNin_tot_svd_k(:,1:80); 
%other method: orth of original SNin_tot, keep first 80
SNin_tot_orth_k=orth(SNin_tot_k);
SNin_tot_orth_k=SNin_tot_orth_k(:,1:80);


% subspace angles 
for i=(1:size(evoked_times,2))
    angle_sVSH_k(i) = subspace(phi_0k(:,i),[SNin_k SNout_k])*180/pi;
    angle_mVSH_k(i)=subspace(phi_0k(:,i),[SNin_tot_k SNout_k])*180/pi;
    angle_mVSH_orth_k(i) = subspace(phi_0k(:,i),[SNin_tot_orth_k SNout_k])*180/pi;
    angle_mVSH_svd_k(i) = subspace(phi_0k(:,i),[SNin_tot_svd_k SNout_k])*180/pi;
end

angle_sVSH_k_min = min(angle_sVSH_k);
angle_sVSH_k_max = max(angle_sVSH_k);
angle_sVSH_k_av = mean(angle_sVSH_k);
angle_mVSH_k_min = min(angle_mVSH_k);
angle_mVSH_k_max = max(angle_mVSH_k);
angle_mVSH_k_av = mean(angle_mVSH_k);
angle_mVSH_k_orth_min = min(angle_mVSH_orth_k);
angle_mVSH_k_orth_max = max(angle_mVSH_orth_k);
angle_mVSH_k_orth_av = mean(angle_mVSH_orth_k);
angle_mVSH_k_svd_min = min(angle_mVSH_svd_k);
angle_mVSH_k_svd_max = max(angle_mVSH_svd_k);
angle_mVSH_k_svd_av = mean(angle_mVSH_svd_k);

return
%% reconstrct internal data- SQUIDS
pS=pinv([SNin SNout]);   
XN=pS*phi_0;
data_rec_sVSH=real(SNin*XN(1:size(SNin,2),:)); 

pS_m=pinv([SNin_tot SNout]);   
XN_m=pS_m*phi_0;
data_rec_mVSH=real(SNin_tot*XN_m(1:size(SNin_tot,2),:)); 

pS_mo=pinv([SNin_tot_orth SNout]);   
XN_mo=pS_mo*phi_0;
data_rec_orth=real(SNin_tot_orth*XN_mo(1:size(SNin_tot_orth,2),:)); 

pS_svd=pinv([SNin_tot_svd SNout]);   
XN_svd=pS_svd*phi_0;
data_rec_svd=real(SNin_tot_svd*XN_svd(1:size(SNin_tot_svd,2),:)); 

%% reconstrct internal data- OPM
%phi-hat
pSp=pinv([SNin_p SNout_p]);   
XNp=pSp*phi_0p;
data_rec_sVSH_p=real(SNin_p*XNp(1:size(SNin_p,2),:)); 

pS_mp=pinv([SNin_tot_p SNout_p]);   
XN_mp=pS_mp*phi_0p;
data_rec_mVSH_p=real(SNin_tot_p*XN_mp(1:size(SNin_tot_p,2),:)); 

pS_mop=pinv([SNin_tot_orth_p SNout_p]);   
XN_mop=pS_mop*phi_0p;
data_rec_orth_p=real(SNin_tot_orth_p*XN_mop(1:size(SNin_tot_orth_p,2),:)); 

pS_svd_p=pinv([SNin_tot_svd_p SNout_p]);   
XN_svd_p=pS_svd_p*phi_0p;
data_rec_svd_p=real(SNin_tot_svd_p*XN_svd_p(1:size(SNin_tot_svd_p,2),:));

%theta-hat
pSt=pinv([SNin_t SNout_t]);   
XNt=pSt*phi_0t;
data_rec_sVSH_t=real(SNin_t*XNt(1:size(SNin_t,2),:)); 

pS_mt=pinv([SNin_tot_t SNout_t]);   
XN_mt=pS_mt*phi_0t;
data_rec_mVSH_t=real(SNin_tot_t*XN_mt(1:size(SNin_tot_t,2),:)); 

pS_mot=pinv([SNin_tot_orth_t SNout_t]);   
XN_mot=pS_mot*phi_0t;
data_rec_orth_t=real(SNin_tot_orth_t*XN_mot(1:size(SNin_tot_orth_t,2),:)); 

pS_svd_t=pinv([SNin_tot_svd_t SNout_t]);   
XN_svd_t=pS_svd_t*phi_0t;
data_rec_svd_t=real(SNin_tot_svd_t*XN_svd_t(1:size(SNin_tot_svd_t,2),:));


%% subspace angles
%SQUID 
for i=(1:size(times,2))
    angle_sVSH(i)=subspace(phi_0(:,i),[SNin SNout])*180/pi;
    angle_mVSH(i)=subspace(phi_0(:,i),[SNin_tot SNout])*180/pi;
    angle_mVSH_orth(i) = subspace(phi_0(:,i),[SNin_tot_orth SNout])*180/pi;
    angle_mVSH_svd(i) = subspace(phi_0(:,i),[SNin_tot_svd SNout])*180/pi;

    angle_sVSH_p(i)=subspace(phi_0p(:,i),[SNin_p SNout_p])*180/pi;
    angle_mVSH_p(i)=subspace(phi_0p(:,i),[SNin_tot_p SNout_p])*180/pi;
    angle_mVSH_orth_p(i) = subspace(phi_0p(:,i),[SNin_tot_orth_p SNout_p])*180/pi;
    angle_mVSH_svd_p(i) = subspace(phi_0p(:,i),[SNin_tot_svd_p SNout_p])*180/pi;
    
    angle_sVSH_t(i)=subspace(phi_0t(:,i),[SNin_t SNout_t])*180/pi;
    angle_mVSH_t(i)=subspace(phi_0t(:,i),[SNin_tot_t SNout_t])*180/pi;
    angle_mVSH_orth_t(i) = subspace(phi_0t(:,i),[SNin_tot_orth_t SNout_t])*180/pi;
    angle_mVSH_svd_t(i) = subspace(phi_0t(:,i),[SNin_tot_svd_t SNout_t])*180/pi;
end

angle_sVSH_min = min(angle_sVSH);
angle_sVSH_max = max(angle_sVSH);
angle_sVSH_av = mean(angle_sVSH);
angle_mVSH_min = min(angle_mVSH);
angle_mVSH_max = max(angle_mVSH);
angle_mVSH_av = mean(angle_mVSH);
angle_mVSH_orth_min = min(angle_mVSH_orth);
angle_mVSH_orth_max = max(angle_mVSH_orth);
angle_mVSH_orth_av = mean(angle_mVSH_orth);
angle_mVSH_svd_min = min(angle_mVSH_svd);
angle_mVSH_svd_max = max(angle_mVSH_svd);
angle_mVSH_svd_av = mean(angle_mVSH_svd);

angle_sVSH_p_min = min(angle_sVSH_p);
angle_sVSH_p_max = max(angle_sVSH_p);
angle_sVSH_p_av = mean(angle_sVSH_p);
angle_mVSH_p_min = min(angle_mVSH_p);
angle_mVSH_p_max = max(angle_mVSH_p);
angle_mVSH_p_av = mean(angle_mVSH_p);
angle_mVSH_p_orth_min = min(angle_mVSH_orth_p);
angle_mVSH_p_orth_max = max(angle_mVSH_orth_p);
angle_mVSH_p_orth_av = mean(angle_mVSH_orth_p);
angle_mVSH_p_svd_min = min(angle_mVSH_svd_p);
angle_mVSH_p_svd_max = max(angle_mVSH_svd_p);
angle_mVSH_p_svd_av = mean(angle_mVSH_svd_p);

angle_sVSH_t_min = min(angle_sVSH_t);
angle_sVSH_t_max = max(angle_sVSH_t);
angle_sVSH_t_av = mean(angle_sVSH_t);
angle_mVSH_t_min = min(angle_mVSH_t);
angle_mVSH_t_max = max(angle_mVSH_t);
angle_mVSH_t_av = mean(angle_mVSH_t);
angle_mVSH_t_orth_min = min(angle_mVSH_orth_t);
angle_mVSH_t_orth_max = max(angle_mVSH_orth_t);
angle_mVSH_t_orth_av = mean(angle_mVSH_orth_t);
angle_mVSH_t_svd_min = min(angle_mVSH_svd_t);
angle_mVSH_t_svd_max = max(angle_mVSH_svd_t);
angle_mVSH_t_svd_av = mean(angle_mVSH_svd_t);

%% subpsace angles - SQUID
% angle_sVSH=subspace(phi_0(:,1),[SNin SNout])*180/pi;
% angle_mVSH=subspace(phi_0(:,1),[SNin_tot SNout])*180/pi;
% angle_mVSH_orth = subspace(phi_0(:,1),[SNin_tot_orth SNout])*180/pi;
% angle_mVSH_svd = subspace(phi_0(:,1),[SNin_tot_svd SNout])*180/pi;
% 
% %% subpsace angles - OPM
% angle_sVSH_p=subspace(phi_0p(:,1),[SNin_p SNout_p])*180/pi;
% angle_mVSH_p=subspace(phi_0p(:,1),[SNin_tot_p SNout_p])*180/pi;
% angle_mVSH_orth_p = subspace(phi_0p(:,1),[SNin_tot_orth_p SNout_p])*180/pi;
% angle_mVSH_svd_p = subspace(phi_0p(:,1),[SNin_tot_svd_p SNout_p])*180/pi;
% 
% angle_sVSH_t=subspace(phi_0t(:,1),[SNin_t SNout_t])*180/pi;
% angle_mVSH_t=subspace(phi_0t(:,1),[SNin_tot_t SNout_t])*180/pi;
% angle_mVSH_orth_t = subspace(phi_0t(:,1),[SNin_tot_orth_t SNout_t])*180/pi;
% angle_mVSH_svd_t = subspace(phi_0t(:,1),[SNin_tot_svd_t SNout_t])*180/pi;


%%check mags vs grads
% j=1;
% k=1;
% for i=(1:size(R,2))
%     if mod(i,3)==0 %every third is a magnetometer
%         SNin_mags(j,:)=SNin(i,:);
%         SNout_mags(j,:)=SNout(i,:);
%         phi_mags(j,:)=phi_0(i,:);
%         j=j+1;
%     else
%         SNin_grads(k,:)=SNin(i,:);
%         SNout_grads(k,:)=SNout(i,:);
%         phi_grads(k,:)=phi_0(i,:);
%         k=k+1;
%     end
% end
% %only mags
% pS_mags=pinv([SNin_mags SNout_mags]);
% XN_mags=pS_mags*phi_mags;
% data_rec_vsh_mags=real(SNin_mags*XN_mags(1:size(SNin_mags,2),:));
% 
% %only grads
% pS_grads=pinv([SNin_grads SNout_grads]);
% XN_grads=pS_grads*phi_grads;
% data_rec_vsh_grads=real(SNin_grads*XN_grads(1:size(SNin_grads,2),:));
% 
% 
% %check subspace angles
% angle_single=subspace(phi_0(:,1),SNin)*180/pi;
% angle_single_full=subspace(phi_0(:,1),[SNin SNout])*180/pi;
% angle_single_mags = subspace(phi_mags(:,1),SNin_mags)*180/pi;
% angle_single_grads = subspace(phi_grads(:,1),SNin_grads)*180/pi;
% 
% angle_multi=subspace(phi_0(:,1),[SNin_tot SNout])*180/pi;
% angle_svd=subspace(phi_0(:,1),SNin_tot_svd)*180/pi;
% angle_orth=subspace(phi_0(:,1),SNin_tot_orth)*180/pi;

return
%% plot data to check
%plot data from single channel
chan_num=3; %change to 1 to plot gradiometer
figure(3);
hold on;
plot(data_times, phi_0(chan_num,:))
plot(data_times, data_rec(chan_num,:))
plot(data_times, data_rec_m(chan_num,:))
plot(data_times, data_rec_orth(chan_num,:))
plot(data_times, data_rec_svd(chan_num,:))
title('testing multi-basis combo methods')
xlabel('time')
ylabel('T')
%ylim([-8e-12 8e-12])
legend({'Raw Data','singleVSH','tSSS method', 'tSSS Orth','eSSS method'},'location','northwest')
hold off
