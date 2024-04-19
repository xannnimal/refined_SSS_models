%% EXAMPLE: SQUID sensors, multi-origin VSH in, single VSH out
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


%%%%%%%%%%%% example from fieldtrip
%--------------------------------------------------------------------------------------
% making a leadfield using the single-sphere headmodel that is
% produced with CTF software
%--------------------------------------------------------------------------------------

% read header, which contains the gradiometer description
hdr  = ft_read_header('Subject01.ds');
grad = hdr.grad;

% read headshape
shape = ft_read_headshape('Subject01.shape');
shape = rmfield(shape, 'fid'); %remove the fiducials->these are stored in MRI-voxel

% read in the single sphere models produced with CTF software
ctf_ss = ft_read_headmodel('Subject01.hdm');

% plotting the head model together with the head shape
ft_plot_sens(grad);
ft_plot_headmodel(ctf_ss, 'facecolor', 'cortex');
ft_plot_headshape(shape);

% prepare the leadfield for the single sphere model
cfg                  = [];
cfg.grad             = grad;
cfg.headmodel        = ctf_ss;
cfg.resolution       = 1;
cfg.unit             = 'cm';
sourcemodel_ctf_ss   = ft_prepare_leadfield(cfg);

% use the same geometry for the grid in what is to follow
sourcemodel = removefields(sourcemodel_ctf_ss, {'leadfield', 'leadfielddimord', 'label'});


return






%% generate SQUID magnetometers
coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";
filename = "C:/Users/xanmc/mne_data/MNE-sample-data/MEG/sample/sample_audvis_raw.fif";
info = fiff_read_meas_info(filename);
[raw] = fiff_setup_read_raw(filename);
[data,times] = fiff_read_raw_segment(raw);
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


%% SSS expansions
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

%% current dipole using Samu's implementation of Sarvas
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
phi_0 = dipole_field_sarvas(rs',q',r0',R,EX,EY,EZ,mags)';

%with fieldtrip
dip_pos = [5,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,1,0]; %(size 3xN
% dip_pos_out = [0,0.25,0]; %[Rx Ry Rz] (size Nx3)
% dip_mom_out = [1,1,1];%(size 3xN)

% for field trip generated data
freq=2;
%dipole_data = single_dipole_sim(R',EZ',dip_pos,dip_mom,freq);
% specify grad
grad = [];
grad.chanpos=R'*100; %convert to cm
grad.coilpos = R'*100;
grad.coilori= EZ'*100; 
grad.senstype = 'meg';
grad.tra= eye(size(R',1));
grad.label = info.ch_names(1,1:306);
% for i=1:size(R',1)
%   grad.label{i} = sprintf('OPM%03d', i);
% end

% create a concentric 3-sphere volume conductor, the radius is the same as for the electrodes
% vol   = [];
% vol.r = 12 * [0.88 0.92 1.00]; % radii of spheres, the head radius is 12 cm
% vol.c = [1 1/80 1];            % conductivity
% vol.o = [0 0 0];               % center of sphere
% vol.r = 10;
% vol.o = [0 0 0];
%The dipoles position and orientation have to be specified with
cfg=[];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, grad);
cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz] (size Nx3)
cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] (size 3xN)
cfg.sourcemodel.unit       = 'cm'; %string, can be 'mm', 'cm', 'm' (default is automatic)
cfg.sourcemodel.frequency = freq;
cfg.unit='cm';
cfg.headmodel     = headmodel; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
cfg.grad          = grad; %structzure with gradiometer definition or filename, see FT_READ_SENS
%dipole_data = ft_dipolesimulation(cfg);
dipole_data = ft_prepare_leadfield(cfg);
% phi_0= dipole_data.trial{1,1}(:,:);
% times = dipole_data.time{1, 1};

return

%old dipole code
% dip_pos_out = [0,0.25,0]; %[Rx Ry Rz] (size Nx3)
% dip_mom = [0,1,1]; %(size 3xN
% dip_mom_out = [1,1,1];%(size 3xN)
% %normalize moments so they have magnitude 1
% dip_mom = dip_mom/norm(dip_mom);
% dip_mom_out = dip_mom_out/norm(dip_mom_out);
% 
% %add time dependence to dipole moment
% f_start = 100; % start frequency
% f_end = 50; % end frequency
% f_start_out = 50; % start frequency
% f_end_out = 30; % end frequency
% timestep = 0.0001;
% T = 0.05;
% rate_of_change = (f_start - f_end)/T;
% rate_of_change_out=(f_start_out-f_end_out)/T;
% times = timestep:timestep:T;
% for i=(1:3)
%     dip_mom_t(i,:) = dip_mom(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
%     dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2));
% end
% 
% %current dipole in, magnetic dipole out
% for i=(1:size(times,2))
%     phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
%     phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)';
%     %phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t(:,i),ch_types)';
% end
% %phi_0=phi_in;

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

%check condition numbers
condition_in = cond(SNin);
condition_in_m = cond(SNin_tot);
condition_in_orth = cond(SNin_tot_orth);
condition_in_svd = cond(SNin_tot_svd);
condition_both = cond([SNin SNout]);
condition_both_m = cond([SNin_tot SNout]);
condition_both_orth = cond([SNin_tot_orth SNout]);
condition_both_svd = cond([SNin_tot_svd SNout]);
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
