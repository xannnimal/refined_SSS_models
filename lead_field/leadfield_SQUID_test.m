%% use FieldTrip leadfield
clear
rawfile = 'sample_audvis_raw.fif'; % from MNE-Python
% need some mri anatomy, can load from an example, but I have the segmented
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
dip_mom = [0,1,0]; 

vol = [];
vol.r = 1.00; %[0.88 0.92 1.00]; % radii of spheres
vol.cond = 1; %[1 1/80 1];       % conductivity
vol.o = [0 0 0];          % center of sphere

% load the sensor poitions from 'sample_audvis_raw.fif'
grad = ft_read_sens(rawfile, 'coordsys', 'dewar', 'senstype', 'meg', 'coilaccuracy', 0); % with coilaccuracy being 0, 1 or 2.
EZ=grad.chanori'; %for reference, these are the data structures used in the SSS calcs
%R=grad.chanpos';

% The positions of the sources can be specified as a regular 3-D
% sourcemodel that is aligned with the axes of the head coordinate system
cfg.xgrid      = 'auto';   %vector (e.g. -20:1:20) or 'auto' (default = 'auto')
cfg.ygrid      = 'auto';   %vector (e.g. -20:1:20) or 'auto' (default = 'auto')
cfg.zgrid      = 'auto';     %vector (e.g.   0:1:20) or 'auto' (default = 'auto')
cfg.singleshell.batchsize = 'all';
%cfg.resolution = 1; %number (e.g. 1 cm) for automatic sourcemodel generation
% or specify a few specific sources of interest
% cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz]
% cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] 
% cfg.sourcemodel.unit       = 'cm'; %string, can be 'mm', 'cm', 'm'
% cfg.sourcemodel.inside = true(size(cfg.sourcemodel.pos,1),1); %assert the dipoles are inside the head
cfg.unit='cm';
cfg.reducerank      = 'no'; %2;
cfg.headmodel     = vol; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
cfg.grad          = grad; %structure with gradiometer definition or filename, see FT_READ_SENS

%generate simulated leadfield
%ft_prepare_leadfield, this calls ft_prepare_vol_sens and ft_compute_leadfield
% if no dipole position/orientation is given, this returns
% 642 cells of nchanx3 lead fields
[leadfield] = ft_prepare_leadfield(cfg);

% need to dot each (nchan x 3) matrix cell into the dipole moment
for i=1:size(leadfield.inside,1)
    for j=1:size(EZ,2)
        sourcemodel(j,i)= dot(leadfield.leadfield{1,i}(j,:), dip_mom');
    end
end

%reduce dimensions using SVD
[Usource,S,V] = svd(sourcemodel,'econ');


return
%[lf] = ft_compute_leadfield(dippos, sens, headmodel, varargin)
sim_data_check = ft_compute_leadfield(dip_pos,grad,vol);

[data] = ft_dipolesimulation(cfg);


% I'm waiting to hear back from Jan to confirm this, but I think the output
% "dipole_data" needs to be dotted into the sensing direction of each
% sensor to get the actual magnetic flux, so like: dipole_data(i,:)*EZ(i,:)