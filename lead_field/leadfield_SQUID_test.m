%% use FieldTrip leadfield
rawfile = 'sample_audvis_raw.fif'; % from MNE-Python
% need some mri anatomy, can load from an example, but I have the egmented
% mri saved as a file so you don't have to uncomment this part
% mri = ft_read_mri('Subject01.mri');
% mri.coordsys = 'ctf'; % just to be sure, it could be that this has been already added by the reading function
% mri = ft_convert_coordsys(mri, 'neuromag');
% cfg           = [];
% cfg.output    = 'brain';
% segmentedmri  = ft_volumesegment(cfg, mri);

% load the preprocessed mri
seg = load("segmentedmri.mat");
segmentedmri = seg.segmentedmri;
cfg = [];
cfg.method='singleshell';
headmodel = ft_prepare_headmodel(cfg, segmentedmri);
headmodel = ft_convert_units(headmodel, 'cm');

% for multiple dipoles, make these a num_dipole x 3 matrix
dip_pos = [5,0,7]; %in centimeters
dip_mom = [1,1,0]; 

% load the sensor poitions from 'sample_audvis_raw.fif'
grad = ft_read_sens(rawfile, 'coordsys', 'dewar', 'senstype', 'meg', 'coilaccuracy', 0); % with coilaccuracy being 0, 1 or 2.
%EZ=grad.chanori'; %for reference, these are the data structures used in the SSS calcs
%R=grad.chanpos';

% create FieldTrip configuration
cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz]
cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] 
cfg.sourcemodel.unit       = 'cm'; %string, can be 'mm', 'cm', 'm'
cfg.sourcemodel.inside = true(size(cfg.sourcemodel.pos,1),1); %assert the dipoles are inside the head
cfg.unit='cm';
cfg.reducerank      = 2;
cfg.headmodel     = headmodel; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
cfg.grad          = grad; %structzure with gradiometer definition or filename, see FT_READ_SENS

%generate simulated leadfield
sim_data = ft_prepare_leadfield(cfg);
dipole_data = sim_data.leadfield{1, 1};

% I'm waiting to hear back from Jan to confirm this, but I think the output
% "dipole_data" needs to be dotted into the sensing direction of each
% sensor to get the actual magnetic flux, so like: dipole_data(i,:)*EZ(i,:)