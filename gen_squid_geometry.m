function [R_mag,EX_mag,EY_mag,EZ_mag, ch_types] = gen_squid_geometry(rawfile, coordsys)
%% from .fiff file generate SQUID magnetometer positions and orientations
%   rawfile: any meg SQUID .fif from MNE, string
%   coordsys: device (or head), string
%   OUTPUT: matrix of R positions (x,y,z) and three sets of nx3 normal coil
%   orientations

% Specify sensor positions/orientations
[R,EX,EY,EZ] = fiff_getpos(rawfile,coordsys);
RT=transpose(R);
EXT=transpose(EX);
EYT=transpose(EY);
EZT=transpose(EZ);

% select every third channel to avoid redundancy
nchan=size(R,2);
mags = 3:3:nchan; % Magnetometer channels
RT_mag = RT(mags,:);
EXT_mag = EXT(mags,:);
EYT_mag = EYT(mags,:);
EZT_mag = EZT(mags,:);
R_mag=RT_mag';
EX_mag=EXT_mag';
EY_mag=EYT_mag';
EZ_mag=EZT_mag';

ch_types = ones(size(EXT_mag,1),1); % 1 = magnetometers,

end