%% SQUID with all recons
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
%% generate SQUID magnetometers
coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";

%for 102 magnetometers, run this
% [R_mag,EX_mag,EY_mag,EZ_mag,ch_types] = gen_squid_geometry(rawfile, coordsys);
% RT_mag=R_mag';
%ch_types = ones(size(EXT_mag,1),1);

%for 306 channels, run this
[R_mag,EX_mag,EY_mag,EZ_mag] = fiff_getpos(rawfile,coordsys);
RT_mag=transpose(R_mag);
EXT_mag=transpose(EX_mag);
EYT_mag=transpose(EY_mag);
EZT_mag=transpose(EZ_mag);
for i=(1:size(EXT_mag,1))
    if mod(i,3)==0 %every third is a magnetometer
        ch_types(i)=1;
    else
        ch_types(i)=0;
    end
end


%% SSS expansions- multi origin interior
%find semi major and minor
%calculate spheroidal in/out
%find semi major and minor- code assumes mm so change units
[semi_major,semi_minor,origin]=find_ellipse_axis(R_mag');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R_mag',EX_mag',EX_mag',EZ_mag',origin,semi_major,semi_minor,Lin,Lout);

for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end

%calculate multi-vsh in and single-vsh out
[SNin_tot,SNout] = multiVSHin_singleVSHout(center1', center2',R_mag,EX_mag,EY_mag,EZ_mag,ch_types,Lin,Lout);
%calculate single in/out
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R_mag,EX_mag,EY_mag,EZ_mag,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R_mag,EX_mag,EY_mag,EZ_mag,ch_types,Lout);

%% generate single dipole simulated data
% dipole_data = single_dipole_sim(R_mag',EZ_mag',dip_pos,dip_mom);
% %pick a specific channel
% phi_0= dipole_data.trial{1,1}(:,:);

dip_pos = [0.05,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,1,1]; %(size 3xN)
%calculate B field
phi_0 = magneticDipole(RT_mag,EZT_mag,dip_pos,dip_mom);


%% reconstrct internal data
%single in, single out
pS=pinv([SNin SNout]);
XN=pS*phi_0;
data_rec_vsh=real(SNin*XN(1:size(SNin,2),:));
%multi in, vsh out
pS_multi_vsh=pinv([SNin_tot SNout]);   
XN_multi_vsh=pS_multi_vsh*phi_0;
data_rec_multi_vsh=real(SNin_tot*XN_multi_vsh(1:size(SNin_tot,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph=pinv([SNin_spm SNout_spm]);   
XN_sph_sph=pS_sph_sph*phi_0;
data_rec_sph_sph=real(SNin_spm*XN_sph_sph(1:size(SNin_spm,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh=pinv([SNin_spm SNout]);   
XN_sph_vsh=pS_sph_vsh*phi_0;
data_rec_sph_vsh=real(SNin_spm*XN_sph_vsh(1:size(SNin_spm,2),:));

%check condition numbers
cond_vsh_vsh=cond([SNin SNout]);
cond_SNin=cond(SNin);
cond_SNin_tot = cond(SNin_tot);
cond_SNout= cond(SNout);
condition_multi_vsh = cond([SNin_tot SNout]);
cond_SNin_spm=cond(SNin_spm);
cond_SNout_spm= cond(SNout_spm);
condition_sph_sph = cond([SNin_spm SNout_spm]);
condition_sph_vsh = cond([SNin_spm SNout]);

check_data_oidin = subspace(phi_0, SNin_spm)*180/pi;
check_data_multi = subspace(phi_0, SNin_tot)*180/pi;
check_data_single = subspace(phi_0, SNin)*180/pi;

%calculate the subspace angle between the reconstructed and noiseless original data for one time instant
% time=10;
% sub_multi_vsh=subspace(phi_0(:,time),data_rec_multi_vsh(:,time))*180/pi;
% sub_sph_sph=subspace(phi_0(:,time),data_rec_sph_sph(:,time))*180/pi;
% sub_sph_vsh=subspace(phi_0(:,time),data_rec_sph_vsh(:,time))*180/pi;
% sub_vsh_vsh=subspace(phi_0(:,time),data_rec_vsh(:,time))*180/pi;

%% plot data to check
%plot data from single channel
% chan_num=1; %3 corresponds to MEG0111
% data_time=dipole_data.time{1,1};
% data_chan_num=dipole_data.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(phi_0)
plot(data_rec_vsh)
plot(data_rec_multi_vsh)
plot(data_rec_sph_sph)
plot(data_rec_sph_vsh)
title('All SSS Methods, Channel 1, Sandia Helmet phi, dipole 1cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off

%% subsapce angles
sVSH_sVSH_t=[SNin SNout];
mVSH_sVSH_t=[SNin_tot SNout];
oid_oid_t=[SNin_spm,SNout_spm];
oid_sVSH_t=[SNin_spm,SNout];
for i=(1:80)
    %mVSH In and Spheroid In
    angles_mVSH_oid_t(i)=subspace(SNin_tot(:,i),SNin_spm)*180/pi;
%sVSH/sVSH and mVSH/sVSH 
    angles_sVSHsVSH_mVSHsVSH_t(i)=subspace(sVSH_sVSH_t(:,i),mVSH_sVSH_t)*180/pi;
%sVSH/sVSH and Spheroid/Spheroid 
    angles_sVSHsVSH_oidoid_t(i)=subspace(sVSH_sVSH_t(:,i),oid_oid_t)*180/pi;
%sVSH/sVSH and Spheroid/sVSH 
    angles_sVSHsVSH_oidsVSH_t(i)=subspace(sVSH_sVSH_t(:,i),oid_sVSH_t)*180/pi;
%mVSH/sVSH and Spheroid/Spheroid  
    angles_mVSHsVSH_oidoid_t(i)=subspace(mVSH_sVSH_t(:,i),oid_oid_t)*180/pi;
%mVSH/sVSH and Spheroid/sVSH
    angles_mVSHsVSH_oidsVSH_t(i)=subspace(mVSH_sVSH_t(:,i),oid_sVSH_t)*180/pi;
%Spheroid/Spheroid and Spheroid/sVSH
    angles_oidoid_oidsVSH_t(i)=subspace(oid_oid_t(:,i),oid_sVSH_t)*180/pi;
end

%find min/max/average
%mVSH In and Spheroid In
max_mVSH_oid_t=max(angles_mVSH_oid_t);
min_mVSH_oid_t=min(angles_mVSH_oid_t);
av_mVSH_oid_t=mean(angles_mVSH_oid_t);
%sVSH/sVSH and mVSH/sVSH 
max_sVSHsVSH_mVSHsVSH_t = max(angles_sVSHsVSH_mVSHsVSH_t);
min_sVSHsVSH_mVSHsVSH_t = min(angles_sVSHsVSH_mVSHsVSH_t);
av_sVSHsVSH_mVSHsVSH_t = mean(angles_sVSHsVSH_mVSHsVSH_t);
%sVSH/sVSH and Spheroid/Spheroid 
max_sVSHsVSH_oidoid_t = max(angles_sVSHsVSH_oidoid_t);
min_sVSHsVSH_oidoid_t = min(angles_sVSHsVSH_oidoid_t);
av_sVSHsVSH_oidoid_t = mean(angles_sVSHsVSH_oidoid_t);
%sVSH/sVSH and Spheroid/sVSH 
max_sVSHsVSH_oidsVSH_t = max(angles_sVSHsVSH_oidsVSH_t);
min_sVSHsVSH_oidsVSH_t = min(angles_sVSHsVSH_oidsVSH_t);
av_sVSHsVSH_oidsVSH_t = mean(angles_sVSHsVSH_oidsVSH_t);
%mVSH/sVSH and Spheroid/Spheroid  
max_mVSHsVSH_oidoid_t = max(angles_mVSHsVSH_oidoid_t);
min_mVSHsVSH_oidoid_t = min(angles_mVSHsVSH_oidoid_t);
av_mVSHsVSH_oidoid_t = mean(angles_mVSHsVSH_oidoid_t);

%mVSH/sVSH and Spheroid/sVSH
max_mVSHsVSH_oidsVSH_t = max(angles_mVSHsVSH_oidsVSH_t);
min_mVSHsVSH_oidsVSH_t = min(angles_mVSHsVSH_oidsVSH_t);
av_mVSHsVSH_oidsVSH_t = mean(angles_mVSHsVSH_oidsVSH_t);

%Spheroid/Spheroid and Spheroid/sVSH
max_oidoid_oidsVSH_t = max(angles_oidoid_oidsVSH_t);
min_oidoid_oidsVSH_t = min(angles_oidoid_oidsVSH_t);
av_oidoid_oidsVSH_t = mean(angles_oidoid_oidsVSH_t);




