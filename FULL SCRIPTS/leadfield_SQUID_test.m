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
lf = ft_prepare_leadfield(cfg);
lf_data = lf.leadfield{1, 1};

%% code from ft_dipolesimulation that uses output of ft_prepare_leadfield
% for i=1:Ndipoles
%       dipsignal{iTr}(i,:) = cos(cfg.sourcemodel.frequency(i)*diptime{iTr}*2*pi + cfg.sourcemodel.phase(i)) * cfg.sourcemodel.amplitude(i);
% end
% for i = 1:3
%     data.trial{trial} = data.trial{trial} + ...
%     lf(:,i:3:end) * (repmat(dipmom{trial}(i:3:end),1,nsamples) .* dipsignal{trial});
% end
dip_mom = dip_mom';
for i=1:3
    sim_data = lf_data(:,i:3:end)*(repmat(dip_mom(i:3:end),1,1));
end

% I'm waiting to hear back from Jan to confirm this, but I think the output
% "dipole_data" needs to be dotted into the sensing direction of each
% sensor to get the actual magnetic flux, so like: dipole_data(i,:)*EZ(i,:)