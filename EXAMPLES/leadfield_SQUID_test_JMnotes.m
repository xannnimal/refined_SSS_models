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


% mri = ft_read_mri('Subject01.mri');
% 
% cfg           = [];
% cfg.output    = 'brain';
% segmentedmri  = ft_volumesegment(cfg, mri);
% 
% cfg = [];
% cfg.method='singleshell';
% headmodel = ft_prepare_headmodel(cfg, segmentedmri);
% headmodel = ft_convert_units(headmodel, 'cm');

% +++++++++++++++++++++
% Note that this MRI volume is defined in CTF-convention coordinate space (ALS), rather than the fif-file’s RAS coordinate system (for the head coordinates).
% I would do the following:

mri = ft_read_mri(Subject01.mri);
mri.coordsys = 'ctf'; % just to be sure, it could be that this has been already added by the reading function
mri = ft_convert_coordsys(mri, 'neuromag');
%++++
cfg           = [];
cfg.output    = 'brain';
segmentedmri  = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method='singleshell';
headmodel = ft_prepare_headmodel(cfg, segmentedmri);
headmodel = ft_convert_units(headmodel, 'cm');



dip_pos = [5,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,1,0]; %(size 3xN
freq=2;
grad = [];
% +++++++++
% w.r.t. the above, I would recommend to do this:
% 
coilaccuracy =0;
grad = ft_read_sens(filename, 'coilaccuracy', coilaccuracy); % with coilaccuracy being 0, 1 or 2.
% 
% This will return a fieldtrip-style definition of the sensor array, and avoid a lot of low-level hassle downstream.
% 
% ================
% grad.chanpos=R'*100; %convert to cm
% grad.coilpos = R'*100;
% grad.coilori= EZ';
% grad.senstype = 'meg';
% grad.tra= eye(size(R',1));
% grad.label = info.ch_names(1,1:306);


cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz] (size Nx3)
cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] (size 3xN)
cfg.sourcemodel.unit       = 'cm'; %string, can be 'mm', 'cm', 'm' (default is automatic)
cfg.sourcemodel.frequency = freq; 
cfg.unit='cm';
cfg.headmodel     = headmodel; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
cfg.grad          = grad; %structzure with gradiometer definition or filename, see FT_READ_SENS
dipole_data = ft_prepare_leadfield(cfg);


%+++++++++++++++++++++
%I would stick with using ft_prepare_leadfield for now. If you would use ft_compute_leadfield, you would need to ensure that all objects are consistently defined w.r.t. each other (mainly with respect to the geometrical units, and some preparatory work with respect to the headmodel handling (e.g. for the ’singleshell method’ the forwpar are computed, as you have already noticed this is not the case when calling ft_compute leadfield directly). As a side note: cfg.sourcemodel.freq is not needed here. Also, I’d recommend to use a dip_pos which is a bit ‘higher’ in the helmet, e.g. [5 0 7]. Fieldtrip’s automatic behavior is to check whether the dipole positions are within the compartment of the headmodel, I am pretty sure that this will be the case, even for the [5 0 7] dipole. Otherwise you may want to explicitly state cfg.sourcemodel.inside = true(size(cfg.sourcemodel.pos,1),1); (% make the line future proof, in case you define more than 1 position….)
