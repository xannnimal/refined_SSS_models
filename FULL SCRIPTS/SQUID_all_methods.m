%% SQUID with all recons
clear
%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% generate SQUID magnetometers
coordsys = 'device'; 
rawfile = 'sample_audvis_raw.fif';
[R,EX,EY,~] = fiff_getpos(rawfile, coordsys);
info = fiff_read_meas_info(rawfile);
% grad = ft_read_sens(rawfile, 'coordsys', 'dewar', 'senstype', 'meg', 'coilaccuracy', 2); % with coilaccuracy being 0, 1 or 2.
% EZ=grad.chanori';
% R=grad.chanpos';

%%% mags only
for i=(1:306)
    ch_types(i)=1;
    mags(i)=i;
end

EZn = load('normal_vectors.mat');
EZ = EZn.EZp;

%%% uncomment for mags/grads
% for i=(1:size(EZ,2))
%     if mod(i,3)==0 %every third is a magnetometer
%         ch_types(i)=1;
%     else
%         ch_types(i)=0;
%     end
% end
% k=1;
% for i=(1:306)
%     if ch_types(i)==1 %every third is a magnetometer
%         mags(k)=i;
%         k=k+1;
%     else
%         k=k;
%     end
% end

%% adjust the height of sensors in top of head
%%% change to move the center of the dipole from 1-7cm away from origin in
%%% xz, pointing in y direction plane
%r0=[0.08,0,0]; %5cm along z axis
% [r0x,r0y,r0z]=sph2cart(0,0.3*pi,0.07);
% r0=[r0x,r0y,r0z];
%%% change to d= 0.01, 0.02. 0.03 to adjust sensors
d=0.00;
%%% change centers of muti-vsh
% center1=[0,0,0];
% [x,y,z]=sph2cart(0,0.3*pi,0.07);
% center2=[x,y,z];

%%% from samu
center1=[0,0,0];
x = 0.07*cos(0.9*pi/2);
z = 0.07*sin(0.9*pi/2);
center2 = [x,0,z];

r2 = 0.04;
x2 = r2*cos(0.9*pi/2);
z2 = r2*sin(0.9*pi/2);
r0 = [x2,0,z2];

for i= (64:84) %(64:84)% (73:84)
    [azimuth,elevation,r] = cart2sph(R(1,i),R(2,i),R(3,i));
    r_new = r-d;
    [vs1,vs2,vs3] = sph2cart(azimuth,elevation,r_new);
    R(1,i)= vs1;
    R(2,i)= vs2;
    R(3,i)= vs3;
end
for i= (112:114)% (73:84)
    [azimuth,elevation,r] = cart2sph(R(1,i),R(2,i),R(3,i));
    r_new = r-d;
    [vs1,vs2,vs3] = sph2cart(azimuth,elevation,r_new);
    R(1,i)= vs1;
    R(2,i)= vs2;
    R(3,i)= vs3;
end
RT=R';
EXT=EX';
EYT=EY';
EZT=EZ';

%% plot sensors
% figure(1)
% grid on;
% hold on;
% rotate3d
% view(135, 20);
% %quiver3(RT(:,1),RT(:,2),RT(:,3),EZT(:,1),EZT(:,2),EZT(:,3))
% %quiver3(RT(:,1),RT(:,2),RT(:,3),EYT(:,1),EYT(:,2),EYT(:,3))
% %quiver3(RT(:,1),RT(:,2),RT(:,3),EXT(:,1),EXT(:,2),EXT(:,3))
% % scatter3(RT(1:63,1),RT(1:63,2),RT(1:63,3),'r')
% % scatter3(RT(68:72,1),RT(68:72,2),RT(68:72,3),'r')
% % scatter3(RT(85:306,1),RT(85:306,2),RT(85:306,3),'r')
% scatter3(RT(1:306,1),RT(1:306,2),RT(1:306,3),'r')
% scatter3(r0(1),r0(2),r0(3),'k*')
% scatter3(RT(64:84,1),RT(64:84,2),RT(64:84,3),'g')
% scatter3(RT(73:84,1),RT(73:84,2),RT(73:84,3),'g')
% scatter3(RT(112:114,1),RT(112:114,2),RT(112:114,3),'g')
% %scatter3(center1(1),center1(2),center1(3),'b*')
% %scatter3(center2(1),center2(2),center2(3),'b*')
% title('Deviation -4cm z, Dipole 7cm from origin')
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% zlabel('z axis (m)')
% hold off;

%% SSS expansions- multi origin interior
%find semi major and minor
%calculate spheroidal in/out
%find semi major and minor
[semi_major,semi_minor,origin]=find_ellipse_axis(R');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EZ',origin,semi_major,semi_minor,Lin,Lout);

for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end


%calculate multi-vsh in and single-vsh out
center1= [-0.00350699, 0.01138051, 0.05947857]; 
center2= [-0.00433911, 0.04081329, 0.05194245]; 
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];
%[SNin_tot,SNout] = multiVSHin_singleVSHout(center1, center2,opm_matrix,R_hat,other_dir,sensing_dir,ch_types,Lin,Lout);
[SNin_tot,~] = multiVSHin_singleVSHout(center1', center2',R,EX,EY,EZ,ch_types,Lin,Lout);

[Sin1,SNin1] = Sin_vsh_vv(center1',R,EX,EY,EZ,ch_types,Lin);
[Sin2,SNin2] = Sin_vsh_vv(center2',R,EX,EY,EZ,ch_types,Lin); %only include l=1
% SNin_tot = orth([SNin1 SNin2]);

%check how close this is to simply performing
% SVD as [S_1,in S_2,in]  = U*Sigma*V^T
%constructing the combined basis from the left singular vectors in U
%as a function of the number of included  u vectors
[U,sigma,Vt] = svd([SNin1 SNin2]);
for i=1:size(SNin2,2)
    check_dif(i) = subspace(U(:,1:120),SNin2(:,i))*180/pi;
    num(i)=i;
end

sig_num = diag(sigma);
ratio = sig_num(120)/sig_num(1); %0.27 percent, could use a threshold of 0.05%

%when we get past 80, U vectos not contained in SNin_tot so subspace jumps
%method for calculating the combined basis is good, matches the SVD
%construct basis from 120 first U vectors - need to motivate some threshold

figure(10)
hold on
plot(num,check_dif)
xlabel('index of SNin2 vector')
ylabel('Subspace angle (degrees)')
title('Between SNin2 vectors and all of U(:,1:120)')
hold off



return


%calculate single in/out
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);



%% generate dependent dipoles
%current dipole using Samu's implementation of Sarvas
rs=[0,0,0];
q=[0,1,0]; %y direction
%r0=[0.05,0,0]; %5cm along x axis
%phi_0 = dipole_field_sarvas(rs',q',r0',R,EX,EY,EZ,mags)';

%add time dependence to dipole moment
dip_mom_out=[1,0,0];
dip_pos_out = [0,0,1.5]; %1.5 meters
f_start = 80; % start frequency
f_end = 30; % end frequency
f_start_out = 50; % start frequency
f_end_out = 30; % end frequency
timestep = 0.0001;
T = 0.1;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;
% 
for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2))*5e8;
end
% 
% %current dipole in, magnetic dipole out
for i=(1:size(times,2))
    %phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
    %phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)'; 
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t_out(:,i),ch_types)'*1e-12;
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)'*1e-12;
end
%phi_0=phi_in+phi_out;
%add gaussian noise at 10 percent of max value of phi_0
noise = randn(size(phi_in,1),size(phi_in,2));
amplitude = 0.15 * phi_in;
%%%% modify this line to do only internal, in+ext, or add noise %%%
phi_0 = phi_in + phi_out; % + amplitude .* noise; %
%%%%
% for i=(1:size(phi_0,1))
%     if mod(i,3)==0 %every third is a magnetometer
%         phi_0(i,:)=phi_0(i,:)*100;
%         phi_in(i,:)=phi_in(i,:)*100;
%         phi_out(i,:)=phi_out(i,:)*100;
%     else
%         phi_0(i,:)=phi_0(i,:);
%         phi_in(i,:)=phi_in(i,:);
%         phi_out(i,:)=phi_out(i,:);
%     end
% end

%% reconstrct internal data
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
% ni=10;
% data_rec_it = xi([SNin_tot,SNout],phi_0,Lin,Lout-1,ni);
% 
% %only mags
% pS_mags=pinv([SNin_mags SNout_mags]);
% XN_mags=pS_mags*phi_mags;
% data_rec_vsh_mags=real(SNin_mags*XN_mags(1:size(SNin_mags,2),:));
% 
% %only grads
% pS_grads=pinv([SNin_grads SNout_grads]);
% XN_grads=pS_grads*phi_grads;
% data_rec_vsh_grads=real(SNin_grads*XN_grads(1:size(SNin_grads,2),:));

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
cond_multi_vsh = cond([SNin_tot SNout]);
cond_SNin_spm=cond(SNin_spm);
cond_SNout_spm= cond(SNout_spm);
cond_sph_sph = cond([SNin_spm SNout_spm]);
cond_sph_vsh = cond([SNin_spm SNout]);

%% subsapce angles
sVSH_sVSH=[SNin SNout];
mVSH_sVSH=[SNin_tot SNout];
oid_oid=[SNin_spm SNout_spm];
oid_sVSH=[SNin_spm SNout];

%check subspace angles- for one time point
% angle_sVSHin=subspace(phi_0(:,1),SNin)*180/pi;
% angle_sVSH_sVSH=subspace(phi_0(:,1),[SNin SNout])*180/pi;
% angle_sVSH_mags = subspace(phi_mags(:,1),SNin_mags)*180/pi;
% angle_sVSH_grads = subspace(phi_grads(:,1),SNin_grads)*180/pi;
% 
% angle_mVSHin = subspace(phi_0(:,1),SNin_tot)*180/pi;
% angle_mVSH_sVSH = subspace(phi_0(:,1),mVSH_sVSH)*180/pi;
% angle_oidin = subspace(phi_0(:,1),SNin_spm)*180/pi;
% angle_oid_oid = subspace(phi_0(:,1),oid_oid)*180/pi;
% angle_oid_sVSH = subspace(phi_0(:,1),oid_sVSH)*180/pi;


%check data for signals with time 
% check_data_vsh_vsh_mags = subspace(phi_mags, [SNin_mags SNout_mags])*180/pi;
% check_data_vsh_vsh_grads = subspace(phi_grads, [SNin_grads SNout_grads])*180/pi;
for i=(1:size(times,2))
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



%% plot data to check
%plot data from single channel
% chan_num=1; %3 corresponds to MEG0111
% data_time=dipole_data.time{1,1};
% data_chan_num=dipole_data.trial{1,1}(chan_num,:); 


chan_num =6;
t_max = 1000;
figure(2);
hold on;
plot(times(:,1:t_max),phi_in(chan_num,1:t_max),'LineWidth',1)
plot(times(:,1:t_max),phi_out(chan_num,1:t_max),'--','LineWidth',1)
%plot(times(:,1:t_max),phi_0(chan_num,1:t_max))
%plot(times,data_rec_vsh_mags(1,:))
%plot(times,data_rec_vsh_grads(1,:))
plot(times(:,1:t_max),data_rec_vsh(chan_num,1:t_max),'LineWidth',1)
plot(times(:,1:t_max),data_rec_multi_vsh(chan_num,1:t_max), 'LineWidth',1)
%plot(times,data_rec_sph_sph(1,:))
%plot(times,data_rec_sph_vsh(1,:))
title('Deviation -4cm, Dipole 5cm')
xlabel('Time (sec)')
ylabel('Dipole Signal (T) chan 6')
%ylim([-8e-12 8e-12])
%legend({'Dip In','VSH Mags','VSH Grads', 'VSH/VSH'},'location','northwest')
legend({'B-in', 'B-ext','sVSH/sVSH','mVSH/sVSH'},'location','northwest')
%legend({'Raw Data','VSH/VSH'},'location','northwest')
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




