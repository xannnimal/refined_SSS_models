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


mri = ft_read_mri('Subject01.mri');
cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas = [ 87   60  116];
cfg.fiducial.lpa = [ 29  145  155];
cfg.fiducial.rpa = [144  142  158];
cfg.coordsys = 'ctf'; % the desired coordinate system
mri_realigned = ft_volumerealign(cfg, mri);

cfg = [];
cfg.method = 'flip';
mri_resliced = ft_volumereslice(cfg, mri_realigned);

cfg           = [];
cfg.output    = {'gray', 'white', 'csf', 'skull', 'scalp'};
segmentedmri  = ft_volumesegment(cfg, mri_resliced);

cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg, segmentedmri);

cfg        = [];
cfg.method = 'simbio';
cfg.conductivity = [1.79 0.33 0.43 0.01 0.14];   % the order follows mesh.tissuelabel, which is 'csf', 'gray', 'scalp', 'skull', 'white'
headmodel  = ft_prepare_headmodel(cfg, mesh);
return


dip_pos = [5,0,0]; %[Rx Ry Rz] (size Nx3)
dip_mom = [0,1,0]; %(size 3xN
freq=2;
grad = [];
grad.chanpos=R'*100; %convert to cm
grad.coilpos = R'*100;
grad.coilori= EZ'*100; 
grad.senstype = 'meg';
grad.tra= eye(size(R',1));
grad.label = info.ch_names(1,1:306);

cfg=[];
cfg.method = 'localspheres';
cfg.grad=grad;
headmodel = ft_prepare_headmodel(cfg);
cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz] (size Nx3)
cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] (size 3xN)
cfg.sourcemodel.unit       = 'cm'; %string, can be 'mm', 'cm', 'm' (default is automatic)
cfg.sourcemodel.frequency = freq;
cfg.unit='cm';
cfg.headmodel     = headmodel; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
cfg.grad          = grad; %structzure with gradiometer definition or filename, see FT_READ_SENS
dipole_data = ft_prepare_leadfield(cfg);