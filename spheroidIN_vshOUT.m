function [SNin_spm,SNout] = spheroidIN_vshOUT(matrix,R_hat,other_dir,sensing_dir,semi_major,semi_minor,Lin,Lout,ch_types)
%% use Tierney et al functions for in spheroidal harmonic expansions
% exterior is single origin vsh
%INPUT
%   matrix: nchanx3 channel position matrix
%   R_hat,phi_hat,sensing_dir: three nchanx3 matricies coil oirientations
%   Semi major/minor: establish spheroidal axis, need to double check good
%   values
%   Lin/Lout: interior and exterior expansion order, AMM paper used (9,2)
%   but we will use (8,3)
%OUT: normalized in spheroidal harmonic, out vsh expansions

%interior spheroidal
chan_ori= [R_hat,other_dir,sensing_dir]; %nchanx9
interior_spm = spm_ipharm(matrix,chan_ori,semi_major,semi_minor,Lin);
Sin_spm=interior_spm;
for j = 1:size(interior_spm,2)
  SNin_spm(:,j) = interior_spm(:,j)/norm(interior_spm(:,j));
end

%exterior single origin vsh
[Sout,SNout] = Sout_vsh_vv([0,0,0]',matrix',R_hat',other_dir',sensing_dir',ch_types,Lout);
