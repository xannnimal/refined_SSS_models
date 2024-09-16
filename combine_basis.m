function SNin_tot = combine_basis(SNin_1,SNin_2)
%% combine two matricies into one basis using spatial overlap
% Xan McPherson, 2023
% overlap method based on tSSS temporal intersection
%INPUT
%   SNin_1,SNin_2: nchanx80 single origin interior normalized VSH
%OUTPUT
%   SNin_tot: (nchan x order) combined interior basis

%find intersection of two bases
Ein_1 = orth(SNin_1); %orth returns an orthonormal basis for range of Bin
Ein_2 = orth(SNin_2);
[Q1,ignore] = qr(Ein_1,0); % upper triangular matrix 'ignore' of the same dimension as Ein and a unitary matrix QA so that X = Q*R.
[Q2,ignore] = qr(Ein_2,0);
[U,S,V] = svd(Q1'*Q2);
Up = Q1*U'; %spatial overlap, already orth
Vp = Q2*V;
%create orthogonal projection operator
P_orth= eye(length(Up))-Up*transpose(Up);
%project out the overlap in SNin_1 and SNin_2
SNin_2_proj = P_orth*SNin_2;
[U,sig,Vt]=svd(SNin_2_proj,"econ"); %econ gives square diag, lowers the amplitude of the reconstruction

% Find the U vectors that significantly differ from SNin_1
for j = 1:size(U,2) 
    kulma_u(j) = subspace(U(:,j),SNin_1);
end
inds = find(kulma_u*180/pi>45); %was 0.05 in radians, change to 40deg
%combined matrix
SNin_tot = [SNin_1 U(:,inds)];
end