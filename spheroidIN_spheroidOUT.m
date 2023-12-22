function [Sin_spm,Sout_spm] = spheroidIN_spheroidOUT(matrix,R_hat,other_dir,sensing_dir,semi_major,semi_minor,Lin,Lout)
%% use Tierney et al functions for spheroidal harmonic expansions
%INPUT
%   matrix: nchanx3 channel position matrix
%   Chan_ori: nchanx9 matrix of three nchanx3 coil position matricies
%   Semi major/minor: establish spheroidal axis, need to double check good
%   values
%   Lin/Lout: interior and exterior expansion order, AMM paper used (9,2)
%   but we will use (8,3)
%OUT: normalized in/out spheroidal harmonic expansions
chan_ori=[R_hat,other_dir,sensing_dir];
interior_spm = spm_ipharm(matrix,chan_ori,semi_major,semi_minor,Lin);
Sin_spm=interior_spm;
% for j = 1:size(interior_spm,2)
%   SNin_spm(:,j) = interior_spm(:,j)/norm(interior_spm(:,j));
% end
exterior_spm = spm_epharm(matrix,chan_ori,semi_major,semi_minor,Lout);
Sout_spm=exterior_spm;
% for j = 1:size(exterior_spm,2)
%   SNout_spm(:,j) = exterior_spm(:,j)/norm(exterior_spm(:,j));
% end

end