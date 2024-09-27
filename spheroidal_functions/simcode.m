
%- matrix of positions provided 
%--------------------------------------------------------------------------

posmat = randn(100,6);
S=[];
S.positions = posmat;
S.wholehead=0;
S.offset=6.5; 
S.nSamples = 1000;
[D] = spm_opm_sim(S);


%- simulate radial point mag, 6.5mm offset and ~30mm spacing
%--------------------------------------------------------------------------

S=[];
S.space=30;
S.wholehead=0;
S.offset=6.5; 
S.nSamples = 1000;
S.axis=1;
[D] = spm_opm_sim(S);


%- simulate dual point mag, 6.5mm offset and ~30mm spacing
%--------------------------------------------------------------------------
S=[];
S.space=30;
S.wholehead=0;
S.offset=6.5; 
S.nSamples = 1000;
S.axis=2;
[D] = spm_opm_sim(S);


%- simulate dual point mag, 6.5mm offset and ~30mm spacing
%--------------------------------------------------------------------------
S=[];
S.space=30;
S.wholehead=0;
S.offset=6.5; 
S.nSamples = 1000;
S.axis=3;
[D] = spm_opm_sim(S);

