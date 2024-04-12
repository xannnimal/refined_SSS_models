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
nchan = size(ch_types,1);

%% SSS expansions- phi
%speficy sensing direction. SQUID=R_hat or EZ, OPM=Theta or phi hat
%single in, single out
[Sin_p,SNin_p] = Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin);
[Sout_p,SNout_p] = Sout_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lout);

%calculate multi-vsh in and single-vsh out
[SNin_tot_p,SNout_p] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin,Lout);
%other methods for mVSH: use svd, first 80 singular values
[~,SNin_1] = Sin_vsh_vv(center1',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin);
[~,SNin_2] = Sin_vsh_vv(center2',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin);
[SNin_tot_svd,sig,~]=svd([SNin_1,SNin_2],'econ');
SNin_tot_svd=SNin_tot_svd(:,1:80); 
%other method: orth of original SNin_tot, keep first 80
SNin_tot_orth=orth(SNin_tot_p);
SNin_tot_orth=SNin_tot_orth(:,1:80);

%calculate spheroidal in/out
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_matrix);
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(opm_matrix,phi_hat,origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm_p,2)
  SNin_spm_p(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm_p(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end

%% simulate dipoles
%current dipole using Samu's implementation of Sarvas
mags = 1:1:nchan;
rs=[0,0,0];
q=[0,1,0]; %y direction
r0=[0.05,0,0]; %5cm along x axis
%phi_0p = dipole_field_sarvas(rs',q',r0',opm_matrix',R_hat',theta_hat',phi_hat',mags)';

dip_mom_out=[1,0,0];
dip_pos_out = [0,5,20];
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
    %phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
    %phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)'; 
    phi_out(:,i) = magneticDipole(opm_matrix',R_hat',theta_hat',phi_hat',dip_pos_out', dip_mom_t_out(:,i),ch_types)'*1e13;
    phi_inp(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',opm_matrix',R_hat',theta_hat',phi_hat',mags)';
end
phi_0p=phi_inp +phi_out;

%field trip dipole sim
%dipole_data_p = single_dipole_sim(opm_matrix,phi_hat,dip_pos,dip_mom);
%pick a specific channel
%phi_0p= dipole_data_p.trial{1,1}(:,:);

%check dipole pos and sensor geometry
% figure(6);
% hold on
% % scatter3(dip_pos(1),dip_pos(2),dip_pos(3), 'r*')
% % scatter3(dip_pos_out(1),dip_pos_out(2),dip_pos_out(3), 'g*')
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3))
% % % scatter3(dip_pos(1),dip_pos(2),dip_pos(3), 'r*')
% % % scatter3(dip_pos_out(1),dip_pos_out(2),dip_pos_out(3), 'g*')
% % scatter3(R(1,:),R(2,:),R(3,:))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), theta_hat(:,1), theta_hat(:,2), theta_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), phi_hat(:,1), phi_hat(:,2), phi_hat(:,3))
% title('SANDIA OPM Geometry')
% %grid on
% rotate3d
% view(135, 20);
% hold off


%% reconstrct internal data
%single in single out
pS_p=pinv([SNin_p, SNout_p]);
XN_p=pS_p*phi_0p;
data_rec_vsh_p=real(SNin_p*XN_p(1:size(SNin_p,2),:));

%multi in, vsh out
pS_multi_vsh_p=pinv([SNin_tot_p SNout_p]);   
XN_multi_vsh_p=pS_multi_vsh_p*phi_0p;
data_rec_multi_vsh_p=real(SNin_tot_p*XN_multi_vsh_p(1:size(SNin_tot_p,2),:)); 

%svd method, truncated at 80
pS_multi_vsh_p_svd=pinv([SNin_tot_svd SNout_p]);   
XN_multi_vsh_p_svd=pS_multi_vsh_p_svd*phi_0p;
data_rec_multi_vsh_svd_p=real(SNin_tot_svd*XN_multi_vsh_p_svd(1:size(SNin_tot_svd,2),:)); 

%orth of SNin_tot, truncated at 80
pS_multi_vsh_p_orth=pinv([SNin_tot_orth SNout_p]);   
XN_multi_vsh_p_orth=pS_multi_vsh_p_orth*phi_0p;
data_rec_multi_vsh_orth_p=real(SNin_tot_orth*XN_multi_vsh_p_orth(1:size(SNin_tot_orth,2),:)); 

%spheroidal in, spheroidal out
pS_sph_sph_p=pinv([SNin_spm_p,SNout_spm_p]);  
XN_sph_sph_p=pS_sph_sph_p*phi_0p;
data_rec_sph_sph_p=real(SNin_spm_p*XN_sph_sph_p(1:size(SNin_spm_p,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_p=pinv([SNin_spm_p SNout_p]);   
XN_sph_vsh_p=pS_sph_vsh_p*phi_0p;
data_rec_sph_vsh_p=real(SNin_spm_p*XN_sph_vsh_p(1:size(SNin_spm_p,2),:));


%check condition numbers
cond_vsh_vsh_p=cond([SNin_p SNout_p]);
cond_SNin_p=cond(SNin_p);
cond_SNin_tot_p = cond(SNin_tot_p);
cond_SNout_p= cond(SNout_p);
condition_multi_vsh_p = cond([SNin_tot_p SNout_p]);
cond_SNin_spm_p=cond(SNin_spm_p);
cond_SNout_spm_p= cond(SNout_spm_p);
condition_sph_sph_p = cond([SNin_spm_p SNout_spm_p]);
condition_sph_vsh_p = cond([SNin_spm_p SNout_p]);

cond_mVSH_svd = cond(SNin_tot_svd);
cond_mVSH_orth = cond(SNin_tot_orth);



%% plot data to check
%plot data from single channel
%chan_num=1; %3 corresponds to MEG0111
%data_time_p=dipole_data_p.time{1,1};
%data_chan_num_p=dipole_data_p.trial{1,1}(chan_num,:); 
% for i=(1:144)
%     chan_num(i)=i;
% end

% figure(2);
hold on;
plot(times, phi_inp(1,:))
plot(times, phi_out(1,:))
plot(times, phi_0p(1,:))
plot(times,data_rec_vsh_p(1,:))
plot(times,data_rec_multi_vsh_p(1,:))
%plot(times,data_rec_sph_sph_p(1,:))
%plot(times,data_rec_sph_vsh_p(1,:))
%title('Raw Data,  Sandia Helmet Phi, dipole 5cm x and 20cm y')
title('Sandia Helmet Phi, dipole 5cm x')
xlabel('time')
ylabel('T')
%ylim([-8e-12 8e-12])
%legend({'B-Dip in','VSH/VSH'},'location','northwest')
legend({'In', 'out','Raw','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off


%%%%%%%%%%%%%%%
%% compare with theta as sensing dir, mVSH
%calculate multi-vsh in and single-vsh out
[SNin_tot_t,SNout_t] = multiVSHin_singleVSHout(center1', center2',opm_matrix',R_hat',phi_hat',theta_hat',ch_types,Lin,Lout);
[Sin_spm_t,Sout_spm_t] = spheroidIN_spheroidOUT(opm_matrix,theta_hat,origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm_t,2)
  SNin_spm_t(:,j) = Sin_spm_t(:,j)/norm(Sin_spm_t(:,j));
end

for j = 1:size(Sout_spm_t,2)
  SNout_spm_t(:,j) = Sout_spm_t(:,j)/norm(Sout_spm_t(:,j));
end


%% generate single dipole simulated data
%phi_0t = dipole_field_sarvas(rs',q',r0',opm_matrix',R_hat',phi_hat',theta_hat',mags)';
for i=(1:size(times,2))
    %phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
    %phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)'; 
    phi_outt(:,i) = magneticDipole(opm_matrix',R_hat',phi_hat',theta_hat',dip_pos_out', dip_mom_t_out(:,i),ch_types)'*1e13;
    phi_int(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',opm_matrix',R_hat',phi_hat',theta_hat',mags)';
end
phi_0t = phi_int+phi_outt;

%old dipole code, broken
%dipole_data_t = single_dipole_sim(opm_matrix,theta_hat,dip_pos,dip_mom);
%pick a specific channel
%phi_0t= dipole_data_t.trial{1,1}(:,:);
%calculate B field
%magneticField = magneticDipole(dip_pos,dip_mom);
%phi_0t = magneticDipole(opm_matrix,theta_hat,phi_hat,dip_pos,dip_mom,ch_types);
%simulate dipoles
% for i=(1:size(times,2))
%     phi_in_t(:,i) = magneticDipole_pointMags(opm_matrix',theta_hat',dip_pos', dip_mom_t(:,i))';
%     phi_out_t(:,i) = magneticDipole_pointMags(opm_matrix',theta_hat',dip_pos_out', dip_mom_t_out(:,i))';
% end
% phi_0t=phi_in_t; %+ phi_out_t;

%% reconstrct internal data
%single in, single out
[Sin_t,SNin_t] = Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',phi_hat',theta_hat',ch_types,Lin);
[Sout_t,SNout_t] = Sout_vsh_vv([0,0,0]',opm_matrix',R_hat',phi_hat',theta_hat',ch_types,Lout);
pS_t=pinv([SNin_t, SNout_t]);
XN_t=pS_t*phi_0t;
data_rec_vsh_t=real(SNin_t*XN_t(1:size(SNin_t,2),:));
%multi in, vsh out
pS_multi_vsh_t=pinv([SNin_tot_t SNout_t]);   
XN_multi_vsh_t=pS_multi_vsh_t*phi_0t;
data_rec_multi_vsh_t=real(SNin_tot_t*XN_multi_vsh_t(1:size(SNin_tot_t,2),:)); 
%spheroidal in, spheroidal out
%pS_sph_sph_t=pinv([Sin_spm_t Sout_spm_t]); 
pS_sph_sph_t=pinv([SNin_spm_t,SNout_spm_t]);  
XN_sph_sph_t=pS_sph_sph_t*phi_0t;
data_rec_sph_sph_t=real(SNin_spm_t*XN_sph_sph_t(1:size(SNin_spm_t,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_t=pinv([SNin_spm_t SNout_t]);   
XN_sph_vsh_t=pS_sph_vsh_t*phi_0t;
data_rec_sph_vsh_t=real(SNin_spm_t*XN_sph_vsh_t(1:size(SNin_spm_t,2),:));

%check condition numbers
cond_vsh_vsh_t=cond([SNin_t SNout_t]);
cond_SNin_t=cond(SNin_t);
cond_SNin_tot_t = cond(SNin_tot_t,2);
cond_SNout_t= cond(SNout_t);
condition_multi_vsh_t = cond([SNin_tot_t SNout_t]);
cond_SNin_spm_t=cond(SNin_spm_t);
cond_SNout_spm_t= cond(Sout_spm_t);
condition_sph_sph_t = cond([SNin_spm_t SNout_spm_t]);
condition_sph_vsh_t = cond([SNin_spm_t SNout_t]);



%% plot data to check
%plot data from single channel
%chan_num=1; %3 corresponds to MEG0111
%data_time_t=dipole_data_t.time{1,1};
%data_chan_num_t=dipole_data_t.trial{1,1}(chan_num,:); 

figure(3);
hold on;
plot(times, phi_int(1,:))
plot(times, phi_outt(1,:))
plot(times, phi_0t(1,:))
plot(times,data_rec_vsh_t(1,:))
plot(times,data_rec_multi_vsh_t(1,:))
%plot(times,data_rec_sph_sph_t(1,:))
%plot(times,data_rec_sph_vsh_t(1,:))
% title('Raw Data, Sandia Helmet Theta, dipole 5cm x and 20cm')
title('Sandia Helmet Theta, dipole 5cm x')
xlabel('time')
ylabel('T')
%legend({'B-Dip In','VSH/VSH'},'location','northwest')
legend({'In','Out','Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off


%% compare sig values
% [U_t,sig_t,Vt_t]=svd(SNin_tot_t,"econ");
% [U_p,sig_p,Vt_p]=svd(SNin_tot_p,"econ");
% sig_t=diag(sig_t);
% sig_p=diag(sig_p);
% figure(4);
% hold on;
% plot(sig_t)
% plot(sig_p)
% title('singlular values multi-VSH in/out')
% legend({'theta','phi'},'location','northwest')
% hold off

% cond_SNin_tot_p_check = max(sig_p)/sig_p(81,1);
% cond_SNin_tot_t_check = max(sig_t)/min(sig_t);



%% calculate subspace angles between bases
% do the collumns of one basis compared to the collumns of the other 
%theta
sVSH_sVSH_t=[SNin_t SNout_t];
mVSH_sVSH_t=[SNin_tot_t SNout_t];
oid_oid_t=[Sin_spm_t,Sout_spm_t];
oid_sVSH_t=[Sin_spm_t,SNout_t];
%phi 
sVSH_sVSH_p=[SNin_p SNout_p];
mVSH_sVSH_p=[SNin_tot_p SNout_p];
oid_oid_p=[Sin_spm_p,Sout_spm_p];
oid_sVSH_p=[Sin_spm_p,SNout_p];

%% compare data reconstructions
%for one time point
% angles_mVSH_svd = subspace(phi_0p, SNin_tot_svd)*180/pi;
% angles_mVSH_orth = subspace(phi_0p, SNin_tot_orth)*180/pi;
% 
% angles_SNin_t = subspace(phi_0t, SNin_t)*180/pi;
% angles_SNin_tot_t= subspace(phi_0t, SNin_tot_t)*180/pi;
% angles_vsh_vsh_t = subspace(phi_0t, sVSH_sVSH_t)*180/pi;
% angles_mvsh_vsh_t = subspace(phi_0t, mVSH_sVSH_t)*180/pi;
% angles_oid_oid_t = subspace(phi_0t, oid_oid_t)*180/pi;
% angles_oid_vsh_t = subspace(phi_0t, oid_sVSH_t)*180/pi;
% 
% angles_SNin_p = subspace(phi_0p, SNin_p)*180/pi;
% angles_SNin_tot_p= subspace(phi_0p, SNin_tot_p)*180/pi;
% angles_vsh_vsh_p = subspace(phi_0p, sVSH_sVSH_p)*180/pi;
% angles_mvsh_vsh_p = subspace(phi_0p, mVSH_sVSH_p)*180/pi;
% angles_oid_oid_p = subspace(phi_0p, oid_oid_p)*180/pi;
% angles_oid_vsh_p = subspace(phi_0p, oid_sVSH_p)*180/pi;


%for data with more than 1 time point
%phi direction
for i=(1:size(times,2))
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
check_data_oid_vsh_pav = mean(check_data_oid_vsh_p);

%theta direction
for i=(1:size(times,2))
    check_data_vsh_vsh_t(i) = subspace(phi_0t(:,i), sVSH_sVSH_t)*180/pi;
    check_data_mvsh_vsh_t(i) = subspace(phi_0t(:,i), mVSH_sVSH_t)*180/pi;
    check_data_oid_oid_t(i) = subspace(phi_0t(:,i), oid_oid_t)*180/pi;
    check_data_oid_vsh_t(i) = subspace(phi_0t(:,i), oid_sVSH_t)*180/pi;
end
check_data_vsh_vsh_tmin = min(check_data_vsh_vsh_t);
check_data_vsh_vsh_tmax = max(check_data_vsh_vsh_t);
check_data_vsh_vsh_tav = mean(check_data_vsh_vsh_t);

check_data_mvsh_vsh_tmin = min(check_data_mvsh_vsh_t);
check_data_mvsh_vsh_tmax = max(check_data_mvsh_vsh_t);
check_data_mvsh_vsh_tav = mean(check_data_mvsh_vsh_t);

check_data_oid_oid_tmin = min(check_data_oid_oid_t);
check_data_oid_oid_tmax = max(check_data_oid_oid_t);
check_data_oid_oid_tav = mean(check_data_oid_oid_t);

check_data_oid_vsh_tmin = min(check_data_oid_vsh_t);
check_data_oid_vsh_tmax = max(check_data_oid_vsh_t);
check_data_oid_vsh_tav = mean(check_data_oid_vsh_t);




return
for i=(1:80)
    %mVSH In and Spheroid In
    angles_mVSH_oid_t(i)=subspace(SNin_tot_t(:,i),Sin_spm_t)*180/pi;
    angles_mVSH_oid_p(i)=subspace(SNin_tot_p(:,i),Sin_spm_p)*180/pi;
%sVSH/sVSH and mVSH/sVSH 
    angles_sVSHsVSH_mVSHsVSH_t(i)=subspace(sVSH_sVSH_t(:,i),mVSH_sVSH_t)*180/pi;
    angles_sVSHsVSH_mVSHsVSH_p(i)=subspace(sVSH_sVSH_p(:,i),mVSH_sVSH_p)*180/pi;
%sVSH/sVSH and Spheroid/Spheroid 
    angles_sVSHsVSH_oidoid_t(i)=subspace(sVSH_sVSH_t(:,i),oid_oid_t)*180/pi;
    angles_sVSHsVSH_oidoid_p(i)=subspace(sVSH_sVSH_p(:,i),oid_oid_p)*180/pi;
%sVSH/sVSH and Spheroid/sVSH 
    angles_sVSHsVSH_oidsVSH_t(i)=subspace(sVSH_sVSH_t(:,i),oid_sVSH_t)*180/pi;
    angles_sVSHsVSH_oidsVSH_p(i)=subspace(sVSH_sVSH_p(:,i),oid_sVSH_p)*180/pi;
%mVSH/sVSH and Spheroid/Spheroid  
    angles_mVSHsVSH_oidoid_t(i)=subspace(mVSH_sVSH_t(:,i),oid_oid_t)*180/pi;
    angles_mVSHsVSH_oidoid_p(i)=subspace(mVSH_sVSH_p(:,i),oid_oid_p)*180/pi;
%mVSH/sVSH and Spheroid/sVSH
    angles_mVSHsVSH_oidsVSH_t(i)=subspace(mVSH_sVSH_t(:,i),oid_sVSH_t)*180/pi;
    angles_mVSHsVSH_oidsVSH_p(i)=subspace(mVSH_sVSH_p(:,i),oid_sVSH_p)*180/pi;
%Spheroid/Spheroid and Spheroid/sVSH
    angles_oidoid_oidsVSH_t(i)=subspace(oid_oid_t(:,i),oid_sVSH_t)*180/pi;
    angles_oidoid_oidsVSH_p(i)=subspace(oid_oid_p(:,i),oid_sVSH_p)*180/pi;
end

%find min/max/average
%mVSH In and Spheroid In
max_mVSH_oid_t=max(angles_mVSH_oid_t);
min_mVSH_oid_t=min(angles_mVSH_oid_t);
av_mVSH_oid_t=mean(angles_mVSH_oid_t);
max_mVSH_oid_p=max(angles_mVSH_oid_p);
min_mVSH_oid_p=min(angles_mVSH_oid_p);
av_mVSH_oid_p=mean(angles_mVSH_oid_p);
%sVSH/sVSH and mVSH/sVSH 
max_sVSHsVSH_mVSHsVSH_t = max(angles_sVSHsVSH_mVSHsVSH_t);
min_sVSHsVSH_mVSHsVSH_t = min(angles_sVSHsVSH_mVSHsVSH_t);
av_sVSHsVSH_mVSHsVSH_t = mean(angles_sVSHsVSH_mVSHsVSH_t);
max_sVSHsVSH_mVSHsVSH_p = max(angles_sVSHsVSH_mVSHsVSH_p);
min_sVSHsVSH_mVSHsVSH_p = min(angles_sVSHsVSH_mVSHsVSH_p);
av_sVSHsVSH_mVSHsVSH_p = mean(angles_sVSHsVSH_mVSHsVSH_p);
%sVSH/sVSH and Spheroid/Spheroid 
max_sVSHsVSH_oidoid_t = max(angles_sVSHsVSH_oidoid_t);
min_sVSHsVSH_oidoid_t = min(angles_sVSHsVSH_oidoid_t);
av_sVSHsVSH_oidoid_t = mean(angles_sVSHsVSH_oidoid_t);
max_sVSHsVSH_oidoid_p = max(angles_sVSHsVSH_oidoid_p);
min_sVSHsVSH_oidoid_p = min(angles_sVSHsVSH_oidoid_p);
av_sVSHsVSH_oidoid_p = mean(angles_sVSHsVSH_oidoid_p);
%sVSH/sVSH and Spheroid/sVSH 
max_sVSHsVSH_oidsVSH_t = max(angles_sVSHsVSH_oidsVSH_t);
min_sVSHsVSH_oidsVSH_t = min(angles_sVSHsVSH_oidsVSH_t);
av_sVSHsVSH_oidsVSH_t = mean(angles_sVSHsVSH_oidsVSH_t);
max_sVSHsVSH_oidsVSH_p = max(angles_sVSHsVSH_oidsVSH_p);
min_sVSHsVSH_oidsVSH_p = min(angles_sVSHsVSH_oidsVSH_p);
av_sVSHsVSH_oidsVSH_p = mean(angles_sVSHsVSH_oidsVSH_p);
%mVSH/sVSH and Spheroid/Spheroid  
max_mVSHsVSH_oidoid_t = max(angles_mVSHsVSH_oidoid_t);
min_mVSHsVSH_oidoid_t = min(angles_mVSHsVSH_oidoid_t);
av_mVSHsVSH_oidoid_t = mean(angles_mVSHsVSH_oidoid_t);
max_mVSHsVSH_oidoid_p = max(angles_mVSHsVSH_oidoid_p);
min_mVSHsVSH_oidoid_p = min(angles_mVSHsVSH_oidoid_p);
av_mVSHsVSH_oidoid_p = mean(angles_mVSHsVSH_oidoid_p);
%mVSH/sVSH and Spheroid/sVSH
max_mVSHsVSH_oidsVSH_t = max(angles_mVSHsVSH_oidsVSH_t);
min_mVSHsVSH_oidsVSH_t = min(angles_mVSHsVSH_oidsVSH_t);
av_mVSHsVSH_oidsVSH_t = mean(angles_mVSHsVSH_oidsVSH_t);
max_mVSHsVSH_oidsVSH_p = max(angles_mVSHsVSH_oidsVSH_p);
min_mVSHsVSH_oidsVSH_p = min(angles_mVSHsVSH_oidsVSH_p);
av_mVSHsVSH_oidsVSH_p = mean(angles_mVSHsVSH_oidsVSH_p);
%Spheroid/Spheroid and Spheroid/sVSH
max_oidoid_oidsVSH_t = max(angles_oidoid_oidsVSH_t);
min_oidoid_oidsVSH_t = min(angles_oidoid_oidsVSH_t);
av_oidoid_oidsVSH_t = mean(angles_oidoid_oidsVSH_t);
max_oidoid_oidsVSH_p = max(angles_oidoid_oidsVSH_p);
min_oidoid_oidsVSH_p = min(angles_oidoid_oidsVSH_p);
av_oidoid_oidsVSH_p = mean(angles_oidoid_oidsVSH_p);


