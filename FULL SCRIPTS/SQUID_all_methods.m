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
rawfile = 'sample_audvis_raw.fif';
[R,EX,EY,EZ] = fiff_getpos(rawfile, coordsys);
info = fiff_read_meas_info(rawfile);
grad = ft_read_sens(rawfile, 'coordsys', 'dewar', 'senstype', 'meg', 'coilaccuracy', 2); % with coilaccuracy being 0, 1 or 2.
EZ=grad.chanori';
R=grad.chanpos';

for i=(1:size(EZ,2))
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
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EZ',origin,semi_major,semi_minor,Lin,Lout);

for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end

%calculate multi-vsh in and single-vsh out
%[SNin_tot,SNout] = multiVSHin_singleVSHout(center1, center2,opm_matrix,R_hat,other_dir,sensing_dir,ch_types,Lin,Lout)
[SNin_tot,~] = multiVSHin_singleVSHout(center1', center2',R,EX,EY,EZ,ch_types,Lin,Lout);
%calculate single in/out
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
%[Sin,SNin] = Sin_basic(rawfile,[0,0,0]',100,coordsys,Lin);
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
r0=[0.05,0,0]; %5cm along x axis
%phi_0 = dipole_field_sarvas(rs',q',r0',R,EX,EY,EZ,mags)';

%add time dependence to dipole moment
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
% 
for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2))*1e9;
end
% 
% %current dipole in, magnetic dipole out
for i=(1:size(times,2))
    %phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
    %phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)'; 
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t_out(:,i),ch_types)';
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)';
end
%phi_0=phi_in+phi_out;
%add gaussian noise at 10 percent of max value of phi_0
noise = randn(size(phi_in,1),size(phi_in,2));
amplitude = 0.15 * phi_in;
%%%% modify this line to do only internal, in+ext, or add noise %%%
phi_0 = phi_in + phi_out + amplitude .* noise; %
%%%%
for i=(1:size(phi_0,1))
    if mod(i,3)==0 %every third is a magnetometer
        phi_0(i,:)=phi_0(i,:)*100;
    else
        phi_0(i,:)=phi_0(i,:);
    end
end

%% use FieldTrip leadfield
% mri = ft_read_mri('Subject01.mri');
% mri.coordsys = 'ctf'; % just to be sure, it could be that this has been already added by the reading function
% mri = ft_convert_coordsys(mri, 'neuromag');
% cfg           = [];
% cfg.output    = 'brain';
% segmentedmri  = ft_volumesegment(cfg, mri);
% save segmentedmri segmentedmri
% seg = load("segmentedmri.mat");
% segmentedmri = seg.segmentedmri;
% cfg = [];
% cfg.method='singleshell';
% headmodel = ft_prepare_headmodel(cfg, segmentedmri);
% headmodel = ft_convert_units(headmodel, 'cm');
% dip_pos = [5,0,7]; %[Rx Ry Rz] (size Nx3)
% dip_mom = [1,1,0]; %(size 3xN
% 
% grad = ft_read_sens(rawfile, 'senstype', 'meg', 'coilaccuracy', 2); % with coilaccuracy being 0, 1 or 2.
% 
% cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz] (size Nx3)
% cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] (size 3xN)
% cfg.sourcemodel.unit       = 'cm'; %string, can be 'mm', 'cm', 'm' (default is automatic)
% cfg.sourcemodel.inside = true(size(cfg.sourcemodel.pos,1),1);
% cfg.unit='cm';
% cfg.reducerank      = 2;
% cfg.headmodel     = headmodel; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
% cfg.grad          = grad; %structzure with gradiometer definition or filename, see FT_READ_SENS
% sim_data = ft_prepare_leadfield(cfg);
% lead_field = sim_data.leadfield{1, 1};
% %multiply the ‘leadfield matrix’, consisting of N-source-component columns 
% % with a matrix of N-source-component rows time courses of activation.
% phi_0 = lead_field*dip_mom';
%check geometry and dip pos
% figure(7);
% hold on
% % scatter3(dip_pos(1),dip_pos(2),dip_pos(3), 'r*')
% % scatter3(dip_pos_out(1),dip_pos_out(2),dip_pos_out(3), 'g*')
% scatter3(R(1,:),R(2,:),R(3,:))
% quiver3(R(1,:),R(2,:),R(3,:), EX(1,:), EX(2,:), EX(3,:))
% %quiver3(R(1,:),R(2,:),R(3,:), EY(1,:), EY(2,:), EY(3,:))
% title('SQUID- EX')
% grid on
% rotate3d
% view(135, 20);
% hold off

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
ni=10;
data_rec_it = xi([SNin_tot,SNout],phi_0,Lin,Lout-1,ni);

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

%check subspace angles- for one time point
angle_sVSHin=subspace(phi_0(:,1),SNin)*180/pi;
angle_sVSH_sVSH=subspace(phi_0(:,1),[SNin SNout])*180/pi;
angle_sVSH_mags = subspace(phi_mags(:,1),SNin_mags)*180/pi;
angle_sVSH_grads = subspace(phi_grads(:,1),SNin_grads)*180/pi;

angle_mVSHin = subspace(phi_0(:,1),SNin_tot)*180/pi;
angle_mVSH_sVSH = subspace(phi_0(:,1),mVSH_sVSH)*180/pi;
angle_oidin = subspace(phi_0(:,1),SNin_spm)*180/pi;
angle_oid_oid = subspace(phi_0(:,1),oid_oid)*180/pi;
angle_oid_sVSH = subspace(phi_0(:,1),oid_sVSH)*180/pi;


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

% figure(2);
chan_num =3;
hold on;
%plot(times(:,1:250),phi_in(chan_num,1:250))
%plot(times(:,1:250),phi_out(chan_num,1:250))
plot(times(:,1:250),phi_0(chan_num,1:250))
%plot(times,data_rec_vsh_mags(1,:))
%plot(times,data_rec_vsh_grads(1,:))
plot(times(:,1:250),data_rec_vsh(chan_num,1:250))
plot(times(:,1:250),data_rec_multi_vsh(chan_num,1:250))
%plot(times,data_rec_sph_sph(1,:))
%plot(times,data_rec_sph_vsh(1,:))
title('SQUID, Currrent Dipole [5cm,0,0], Magnetic Dipole [0,0,1.5m]')
xlabel('Time (sec)')
ylabel('Dipole Signal, Chan 1 (T)')
%ylim([-8e-12 8e-12])
%legend({'Dip In','VSH Mags','VSH Grads', 'VSH/VSH'},'location','northwest')
legend({'Raw Data','VSH/VSH','Multi/VSH'},'location','northwest')
%legend({'Raw Data','VSH/VSH'},'location','northwest')
hold off



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




