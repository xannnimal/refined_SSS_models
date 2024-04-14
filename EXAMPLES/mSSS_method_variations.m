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
condition_in = cond(SNin);
condition_in_m = cond(SNin_tot);
condition_in_orth = cond(SNin_tot_orth);
condition_in_svd = cond(SNin_tot_svd);
condition_both = cond([SNin SNout]);
condition_both_m = cond([SNin_tot SNout]);
condition_both_orth = cond([SNin_tot_orth SNout]);
condition_both_svd = cond([SNin_tot_svd SNout]);


return
%% simulate current dipole
%current dipole using Samu's implementation of Sarvas
rs=[0,0,0];
q=[0,1,0]; %y direction
r0=[0.05,0,0]; %5cm along x axis
phi_0_squid = dipole_field_sarvas(rs',q',r0',R,EX,EY,EZ,mags)';
phi_0p = dipole_field_sarvas(rs',q',r0',opm_matrix',R_hat',theta_hat',phi_hat',mags_opm)';


%% reconstrct internal data
pS=pinv([SNin SNout]);   
XN=pS*phi_0;
data_rec=real(SNin*XN(1:size(SNin,2),:)); 

pS_m=pinv([SNin_tot SNout]);   
XN_m=pS_m*phi_0;
data_rec_m=real(SNin_tot*XN_m(1:size(SNin_tot,2),:)); 

pS_mo=pinv([SNin_tot_orth SNout]);   
XN_mo=pS_mo*phi_0;
data_rec_orth=real(SNin_tot_orth*XN_mo(1:size(SNin_tot_orth,2),:)); 

pS_svd=pinv([SNin_tot_svd SNout]);   
XN_svd=pS_svd*phi_0;
data_rec_svd=real(SNin_tot_svd*XN_svd(1:size(SNin_tot_svd,2),:)); 

%%check mags vs grads
j=1;
k=1;
for i=(1:size(R,2))
    if mod(i,3)==0 %every third is a magnetometer
        SNin_mags(j,:)=SNin(i,:);
        SNout_mags(j,:)=SNout(i,:);
        phi_mags(j,:)=phi_0(i,:);
        j=j+1;
    else
        SNin_grads(k,:)=SNin(i,:);
        SNout_grads(k,:)=SNout(i,:);
        phi_grads(k,:)=phi_0(i,:);
        k=k+1;
    end
end
%only mags
pS_mags=pinv([SNin_mags SNout_mags]);
XN_mags=pS_mags*phi_mags;
data_rec_vsh_mags=real(SNin_mags*XN_mags(1:size(SNin_mags,2),:));

%only grads
pS_grads=pinv([SNin_grads SNout_grads]);
XN_grads=pS_grads*phi_grads;
data_rec_vsh_grads=real(SNin_grads*XN_grads(1:size(SNin_grads,2),:));


%check subspace angles
angle_single=subspace(phi_0(:,1),SNin)*180/pi;
angle_single_full=subspace(phi_0(:,1),[SNin SNout])*180/pi;
angle_single_mags = subspace(phi_mags(:,1),SNin_mags)*180/pi;
angle_single_grads = subspace(phi_grads(:,1),SNin_grads)*180/pi;

angle_multi=subspace(phi_0(:,1),[SNin_tot SNout])*180/pi;
angle_svd=subspace(phi_0(:,1),SNin_tot_svd)*180/pi;
angle_orth=subspace(phi_0(:,1),SNin_tot_orth)*180/pi;

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
