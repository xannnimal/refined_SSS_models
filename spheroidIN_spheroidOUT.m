function [Sin_spm,Sout_spm] = spheroidIN_spheroidOUT(matrix,R_hat,other_dir,sensing_dir,o,semi_major,semi_minor,Lin,Lout)
%% use Tierney et al functions for spheroidal harmonic expansions
%INPUT
%   matrix: nchanx3 channel position matrix
%   Chan_ori: nchanx9 matrix of three nchanx3 coil position matricies
%   Semi major/minor: establish spheroidal axis, need to double check good
%   values
%   o: origin from find_ellipse
%   Lin/Lout: interior and exterior expansion order, AMM paper used (9,2)
%   but we will use (8,3)
%OUT: normalized in/out spheroidal harmonic expansions

chan_ori=[R_hat,other_dir,sensing_dir];
vtest = double(bsxfun(@minus,matrix,o'));
Sin_spm = spm_ipharm(vtest,chan_ori,semi_major,semi_minor,Lin);
Sout_spm = spm_epharm(vtest,chan_ori,semi_major,semi_minor,Lout);


end