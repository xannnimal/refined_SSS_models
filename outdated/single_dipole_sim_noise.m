function dipole_data = single_dipole_sim_noise(opm_matrix,sensing_dir,dip_pos,dip_mom, noise)
%% simulate data from a single dipole using cfg from field trip
% Xan McPherson, 2023
%   opm_matrix: nchanx3 matrix of (x,y,z) channel corrdinates
%   sensing_dir: nchanx3 matrix of coilpos for sensing, R theta or phi hat
%   dip_pos: [Rx Ry Rz] size Nx3 position of dipole
%   dip_mom: [Qx Qy Qz] size 3xN moment of dipole
%   cfg.relnoise (1e-4)    = add noise= relnoise#*rms of data
%   cfg.absnoise    = add noise with absolute level = absnoise#
%   cfg.randomseed  = 'yes' or a number or vector with the seed value (default = 'yes')

% specify grad
grad = [];
grad.coilpos = opm_matrix;
grad.coilori= sensing_dir; 
grad.senstype = 'meg';
grad.tra= eye(size(opm_matrix,1));
for i=1:size(opm_matrix,1)
  grad.label{i} = sprintf('OPM%03d', i);
end

% specify cfg using "sourcemodel"
vol.r = 10;
vol.o = [0 0 0];
%The dipoles position and orientation have to be specified with
cfg.sourcemodel.pos        = dip_pos; %[Rx Ry Rz] (size Nx3)
cfg.sourcemodel.mom        = dip_mom; %[Qx Qy Qz] (size 3xN)
cfg.sourcemodel.unit       = 'm'; %string, can be 'mm', 'cm', 'm' (default is automatic)
cfg.relnoise = noise;
cfg.randomseed             = 1; %set for reprodicibility 
cfg.headmodel     = vol; %structure with volume conduction model, see FT_PREPARE_HEADMODEL
cfg.grad          = grad; %structure with gradiometer definition or filename, see FT_READ_SENS

dipole_data = ft_dipolesimulation(cfg);

end