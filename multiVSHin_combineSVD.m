function [SNin_tot] = multiVSHin_combineSVD(center1, center2,chanpos,ori1,ori2,sensing_dir,ch_types,Lin)
%% Calculate multi-origin vsh interior and single-origin vsh out
% Xan McPherson, 2023
% give two centers, code calculates two SSS expansions, then combines them
% into one basis using "combine_basis.m"
% INPUT
%   center1, center2: 3x1 (x,y,z) locations of expansion centers
%   chanpos: 3xnchan matrix of sensor locations
%   ori1,ori2,sensing_dir: three 3xnchan normal vectors for coil
%       orientations. For example, input (ori1=EX,ori2=EY,sensing_dir=EZ) for SQUID
%   ch_types: 1xnchan vector of 1's for magnetometers
%   Lin,Lout: vsh truncation order, typically (8,3)
% OUTPUT 
%   SNin_tot: nchan x order matrix of norm combined interior expansion

% third input is the sensing direction.
[~,SNin_1] = Sin_vsh_vv(center1,chanpos,ori1,ori2,sensing_dir,ch_types,Lin);
[~,SNin_2] = Sin_vsh_vv(center2,chanpos,ori1,ori2,sensing_dir,ch_types,Lin);
%combine using SVD technique from eSSS paper 
%if matrix is square, the orthonormal basis calculated by orth(A) matches the matrix U calculated in the singular value decomposition [U,S] = svd(A,"econ")
[U,sig,Vt]=svd([SNin_1,SNin_2]);


SNin_tot=U;
end