%% modified sensor locations to test mSSS
% simulate a source which has difficulties in SSS and show that mSSS helps
% modify sensor locations from 306, simulate dipole somewhere, move some of
% the sensors closer encroaching into the SSS single sphere area 
% "in principle this approach works and is a good solutions to cases where
% we have this problem in on-scalp systems, next step is to use real data"
% dipole on z-axis, move top sensors in

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
filename = "C:/Users/xanmc/mne_data/MNE-sample-data/MEG/sample/sample_audvis_raw.fif";
[R,EX,EY,EZ] = fiff_getpos(filename, coordsys);
info = fiff_read_meas_info(filename);
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
%% adjust the height of sensors in top of head
d=0.03;
for i= (73:84) %(64:84)% (73:84)
    [azimuth,elevation,r] = cart2sph(R(1,i),R(2,i),R(3,i));
    r_new = r-d;
    [vs1,vs2,vs3] = sph2cart(azimuth,elevation,r_new);
    R(1,i)= vs1;
    R(2,i)= vs2;
    R(3,i)= vs3;
end
% for i= (112:114)% (73:84)
%     [azimuth,elevation,r] = cart2sph(R(1,i),R(2,i),R(3,i));
%     r_new = r-d;
%     [vs1,vs2,vs3] = sph2cart(azimuth,elevation,r_new);
%     R(1,i)= vs1;
%     R(2,i)= vs2;
%     R(3,i)= vs3;
% end
RT=R';

r0=[0,0,0.07]; %5cm along z axis

%% plot sensors
% figure(1)
% grid on;
% hold on;
% rotate3d
% view(135, 20);
% % quiver3(RT(:,1),RT(:,2),RT(:,3),EZT(:,1),EZT(:,2),EZT(:,3))
% % quiver3(RT(:,1),RT(:,2),RT(:,3),EYT(:,1),EYT(:,2),EYT(:,3))
% % quiver3(RT(:,1),RT(:,2),RT(:,3),EXT(:,1),EXT(:,2),EXT(:,3))
% % scatter3(RT(1:63,1),RT(1:63,2),RT(1:63,3),'r')
% % scatter3(RT(68:72,1),RT(68:72,2),RT(68:72,3),'r')
% % scatter3(RT(85:306,1),RT(85:306,2),RT(85:306,3),'r')
% scatter3(RT(1:306,1),RT(1:306,2),RT(1:306,3),'r')
% scatter3(r0(1),r0(2),r0(3),'b*')
% %scatter3(RT(64:69,1),RT(64:69,2),RT(64:69,3),'g')
% %scatter3(RT(70:72,1),RT(70:72,2),RT(70:72,3),'g')
% scatter3(RT(73:84,1),RT(73:84,2),RT(73:84,3),'g')
% %scatter3(RT(112:114,1),RT(112:114,2),RT(112:114,3),'g')
% %scatter3(center1(1),center1(2),center1(3),'g*')
% %scatter3(center2(1),center2(2),center2(3),'b*')
% title('306 Neuromag Helmet, Current Dipole=[0,0,0.05m]')
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% zlabel('z axis (m)')
% hold off;

%% SSS expansions- multi origin interior
%spheroidal
[semi_major,semi_minor,origin]=find_ellipse_axis(R');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EZ',origin,semi_major,semi_minor,Lin,Lout);
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


%% generate dependent dipoles
%current dipole using Samu's implementation of Sarvas
k=1;
for i=(1:306)
    if ch_types(i)==1 %every third is a magnetometer
        mags(k)=i;
        k=k+1;
    else
        k=k;
    end
end
rs=[0,0,0];
q=[0,1,0]; %y direction


dip_mom_out=[1,0,0];
dip_pos_out = [0,0,1.5]; %1.5 meters
f_start = 100; % start frequency
f_end = 50; % end frequency
f_start_out = 50; % start frequency
f_end_out = 30; % end frequency
timestep = 0.0001;
T = 0.05;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;

for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2))*1e9;
end

%current dipole in, magnetic dipole out
for i=(1:size(times,2)) 
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t_out(:,i),ch_types)';
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)';
end
%add gaussian noise at 10 percent of max value of phi_0
noise = randn(size(phi_in,1),size(phi_in,2));
% Create an amplitude for that noise that is 10% of the noise-free signal at every element.
amplitude = 0.15 * phi_in;
% Now add the noise-only signal to your original noise-free signal to create a noisy signal.
phi_0 = phi_in + amplitude .* noise + phi_out;

%% reconstruct data
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

figure(1)
hold on
plot(times,phi_0(1,:))
plot(times,data_rec_vsh(1,:))
plot(times,data_rec_multi_vsh(1,:))
title('SQUID, Currrent Dipole [5cm,0,0], Magnetic Dipole [0,0,1.5m]')
xlabel('Time (sec)')
ylabel('Dipole Signal, Chan 1 (T)')
legend({'Raw Data','VSH/VSH','Multi/VSH'},'location','northwest')