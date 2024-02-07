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


%% SSS expansions- multi origin interior
%find semi major and minor
%calculate spheroidal in/out
%find semi major and minor
[semi_major,semi_minor,origin]=find_ellipse_axis(R');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EX',EX',EZ',origin,semi_major,semi_minor,Lin,Lout);

for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end

%calculate multi-vsh in and single-vsh out
[SNin_tot,~] = multiVSHin_singleVSHout(center1', center2',R,EX,EY,EZ,ch_types,Lin,Lout);
%calculate single in/out
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);


%% generate single dipole simulated data
% for fielf trip generated data
% dipole_data = single_dipole_sim(R_mag',EZ_mag',dip_pos,dip_mom);
% phi_0= dipole_data.trial{1,1}(:,:);

dip_pos = [0.05,0,0]; %[Rx Ry Rz] (size Nx3)
dip_pos_out = [0.15,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,1,1]; %(size 3xN)
%phi_in= magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom',ch_types)';

%add time dependence to dipole moment
f_start = 100; % start frequency
f_end = 50; % end frequency
timestep = 0.0001;
T = 0.05;
rate_of_change = (f_start - f_end)/T;
times = timestep:timestep:T;
for i=(1:3)
    dip_mom_t(i,:) = dip_mom(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
end

%simulate dipoles
for i=(1:size(times,2))
    phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)';
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t(:,i),ch_types)';
end
phi_0=phi_in;

%using Samu's function
% m=dip_mom';
% r0=dip_pos';
% NX=EX;
% NY=EY;
% NZ=EZ;
% magtype='gegen';
% k=1;
% for i=(1:306)
%     if ch_types(i)==1 %every third is a magnetometer
%         mags(k)=i;
%         k=k+1;
%     else
%         k=k;
%     end
% end
% vuovec = mkentta(m,r0,R,NX,NY,NZ,mags,magtype);
% phi_0=vuovec';
%% reconstrct internal data
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

%% check condition numbers
cond_vsh_vsh=cond([SNin SNout]);
cond_SNin=cond(SNin);
cond_SNin_tot = cond(SNin_tot);
cond_SNout= cond(SNout);
condition_multi_vsh = cond([SNin_tot SNout]);
cond_SNin_spm=cond(SNin_spm);
cond_SNout_spm= cond(SNout_spm);
condition_sph_sph = cond([SNin_spm SNout_spm]);
condition_sph_vsh = cond([SNin_spm SNout]);

%% subsapce angles
sVSH_sVSH=[SNin SNout];
mVSH_sVSH=[SNin_tot SNout];
oid_oid=[SNin_spm,SNout_spm];
oid_sVSH=[SNin_spm,SNout];

check_data_vsh_vsh_mags = subspace(phi_mags, [SNin_mags SNout_mags])*180/pi;
check_data_vsh_vsh_grads = subspace(phi_grads, [SNin_grads SNout_grads])*180/pi;
for i=(1:306)
    check_data_vsh_vsh_d(i) = subspace(phi_0(:,i), sVSH_sVSH)*180/pi;
    check_data_mvsh_vsh_d(i) = subspace(phi_0(:,i), mVSH_sVSH)*180/pi;
    check_data_oid_oid_d(i) = subspace(phi_0(:,i), oid_oid)*180/pi;
    check_data_oid_vsh_d(i) = subspace(phi_0(:,i), oid_sVSH)*180/pi;
end
check_data_vsh_vsh_dmin = min(check_data_vsh_vsh_d);
check_data_vsh_vsh_dmax = max(check_data_vsh_vsh_d);
check_data_vsh_vsh_dav = mean(check_data_vsh_vsh_d);

check_data_mvsh_vsh_dmin = min(check_data_mvsh_vsh_d);
check_data_mvsh_vsh_dmax = max(check_data_mvsh_vsh_d);
check_data_mvsh_vsh_dav = mean(check_data_mvsh_vsh_d);

check_data_oid_oid_dmin = min(check_data_oid_oid_d);
check_data_oid_oid_dmax = max(check_data_oid_oid_d);
check_data_oid_oid_dav = mean(check_data_oid_oid_d);

check_data_oid_vsh_dmin = min(check_data_oid_vsh_d);
check_data_oid_vsh_dmax = max(check_data_oid_vsh_d);
check_data_oid_vsh_dav = mean(check_data_oid_vsh_d);

return

%% plot data to check
%plot data from single channel
% chan_num=1; %3 corresponds to MEG0111
% data_time=dipole_data.time{1,1};
% data_chan_num=dipole_data.trial{1,1}(chan_num,:); 

figure(2);
hold on;
plot(times,phi_0(1,:))
plot(times,data_rec_vsh(1,:))
%plot(times,data_rec_multi_vsh(1,:))
%plot(times,data_rec_sph_sph(1,:))
%plot(times,data_rec_sph_vsh(1,:))
title('All SSS Methods, Channel 1, Sandia Helmet phi, dipole 5cm x')
xlabel('time')
ylabel('MEG0121')
%ylim([-8e-12 8e-12])
legend({'Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off


return
%% quantify differences between SSS methods using subspace angles
%calculations independent of the type of raw data

for i=(1:80)
    %mVSH In and Spheroid In
    angles_mVSH_oid(i)=subspace(SNin_tot(:,i),SNin_spm)*180/pi;
%sVSH/sVSH and mVSH/sVSH 
    angles_sVSHsVSH_mVSHsVSH(i)=subspace(sVSH_sVSH(:,i),mVSH_sVSH)*180/pi;
%sVSH/sVSH and Spheroid/Spheroid 
    angles_sVSHsVSH_oidoid(i)=subspace(sVSH_sVSH(:,i),oid_oid)*180/pi;
%sVSH/sVSH and Spheroid/sVSH 
    angles_sVSHsVSH_oidsVSH(i)=subspace(sVSH_sVSH(:,i),oid_sVSH)*180/pi;
%mVSH/sVSH and Spheroid/Spheroid  
    angles_mVSHsVSH_oidoid(i)=subspace(mVSH_sVSH(:,i),oid_oid)*180/pi;
%mVSH/sVSH and Spheroid/sVSH
    angles_mVSHsVSH_oidsVSH(i)=subspace(mVSH_sVSH(:,i),oid_sVSH)*180/pi;
%Spheroid/Spheroid and Spheroid/sVSH
    angles_oidoid_oidsVSH(i)=subspace(oid_oid(:,i),oid_sVSH)*180/pi;
end

%find min/max/average
%mVSH In and Spheroid In
max_mVSH_oid_t=max(angles_mVSH_oid);
min_mVSH_oid_t=min(angles_mVSH_oid);
av_mVSH_oid_t=mean(angles_mVSH_oid);
%sVSH/sVSH and mVSH/sVSH 
max_sVSHsVSH_mVSHsVSH_t = max(angles_sVSHsVSH_mVSHsVSH);
min_sVSHsVSH_mVSHsVSH_t = min(angles_sVSHsVSH_mVSHsVSH);
av_sVSHsVSH_mVSHsVSH_t = mean(angles_sVSHsVSH_mVSHsVSH);
%sVSH/sVSH and Spheroid/Spheroid 
max_sVSHsVSH_oidoid_t = max(angles_sVSHsVSH_oidoid);
min_sVSHsVSH_oidoid_t = min(angles_sVSHsVSH_oidoid);
av_sVSHsVSH_oidoid_t = mean(angles_sVSHsVSH_oidoid);
%sVSH/sVSH and Spheroid/sVSH 
max_sVSHsVSH_oidsVSH_t = max(angles_sVSHsVSH_oidsVSH);
min_sVSHsVSH_oidsVSH_t = min(angles_sVSHsVSH_oidsVSH);
av_sVSHsVSH_oidsVSH_t = mean(angles_sVSHsVSH_oidsVSH);
%mVSH/sVSH and Spheroid/Spheroid  
max_mVSHsVSH_oidoid_t = max(angles_mVSHsVSH_oidoid);
min_mVSHsVSH_oidoid_t = min(angles_mVSHsVSH_oidoid);
av_mVSHsVSH_oidoid_t = mean(angles_mVSHsVSH_oidoid);

%mVSH/sVSH and Spheroid/sVSH
max_mVSHsVSH_oidsVSH_t = max(angles_mVSHsVSH_oidsVSH);
min_mVSHsVSH_oidsVSH_t = min(angles_mVSHsVSH_oidsVSH);
av_mVSHsVSH_oidsVSH_t = mean(angles_mVSHsVSH_oidsVSH);

%Spheroid/Spheroid and Spheroid/sVSH
max_oidoid_oidsVSH_t = max(angles_oidoid_oidsVSH);
min_oidoid_oidsVSH_t = min(angles_oidoid_oidsVSH);
av_oidoid_oidsVSH_t = mean(angles_oidoid_oidsVSH);




