function [SNin_tot,SNout] = multiVSHin_singleVSHout(center1, center2,opm_matrix,R_hat,other_dir,sensing_dir,ch_types,Lin,Lout)
%% Calculate multi-origin vsh interior and single-origin vsh out
% Xan McPherson, 2023
% give two centers, code calculates two SSS expansions, then combines them
% into one basis using "combine_basis.m"
% INPUT
%   center1, center2: 3x1 (x,y,z) locations of expansion centers
%   opm_matrix: 3xnchan matrix of sensor locations
%   R_hat,other_dir,sensing_dir: three 1xnchan normal vectors for coil
%       orientations
%   ch_types: 1xnchan vector of 1's for magnetometers
%   Lin,Lout: vsh truncation order, typically (8,3)
% OUTPUT 
%   SNin_tot: nchan x order matrix of norm combined interior expansion
%   SNout: nchan x 15 matrix norm exterior expansion

% third input is the sensing direction. should be theta or phi
[Sin_1,SNin_1] = Sin_vsh_vv(center1,opm_matrix,R_hat,other_dir,sensing_dir,ch_types,Lin);
[Sin_2,SNin_2] = Sin_vsh_vv(center2,opm_matrix,R_hat,other_dir,sensing_dir,ch_types,Lin);
%combine using "combine_basis"
SNin_tot = combine_basis(SNin_1,SNin_2);
%external using single origin expansion
[Sout,SNout] = Sout_vsh_vv([0,0,0]',opm_matrix,R_hat,other_dir,sensing_dir,ch_types,Lout);

end