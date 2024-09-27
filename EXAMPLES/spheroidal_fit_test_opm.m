%% spheroid fit test with OPM helmet
%after meeting with Tim, we need to transpose the Sandia matrix so the y
% axis is the longest, so this code does that. Second, we need to substract
% the origin of the optimized spheroid from the matrix geometry to center
% it and pass the centered, optimized coordinated into the spheroidal basis
% calculation functions (My function "SpheroidIN_spheroidOUT" does this
% automatically)
clear
%% opm geometry 
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,EZ,EX,EY,ch_types] = gen_opm_geometry(filename);
nchan = size(ch_types,1);

%% SSS expansions- multi origin interior
%speficy sensing direction
sensing_dir=EZ;
other_dir2 = EY;
other_dir1 = EX;
%find major and minor axis of spheroidal ellipse
%Y must be the longest axis of opm_matrix, coords given have
%helmet flipped upside down
opm_mod(:,1)= opm_matrix(:,1);
opm_mod(:,2)= opm_matrix(:,3);
opm_mod(:,3)= opm_matrix(:,2);
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_mod);

%%%%%%%%%%%%%%%%%%%%
function [a,b,o] = find_ellipse_axis(chanpos)
%code from Tierney et al "spm_opm_amm.m" https://github.com/spm/spm/blob/main/spm_opm_amm.m
% modified by Xan McPherson 2023 to calculate the semi major and minor axis of an ellipse given
% channel positions for point magnetometers
%INPUT
%   chanpos= channel positions (nc x 3 matrix)   
%OUTPUT
%   a= semi-major axis
%   b= semi-minor axis
v = chanpos;
vrange = abs((max(v)-min(v)));
[~,ind]=max(vrange);
if ind==1
    [ o, r]=spheroid_fit(v,1);
end

if ind==2
    [ o, r ]=spheroid_fit(v,2);
end

if ind==3
    [ o, r]=spheroid_fit(v,3);
end

% if (ind~=2)
%     error('Y is not longest axis.... fix please')
% end


inside = v(:,1).^2/r(1)^2+v(:,2).^2/r(2)^2+v(:,3).^2/r(3)^2;
c = sum(inside<1);
while c>0
  rt = r-0.001; %-0.001
  inside = v(:,1).^2/rt(1)^2+v(:,2).^2/rt(2)^2+v(:,3).^2/rt(3)^2;
  cc = sum(inside<1);
  if(cc<=c)
    r = r-0.001;  %changed from -1 for meter input, not mm -0.001
    c = cc;
  end 
end

[X,Y,Z]=ellipsoid(o(1),o(2),o(3),r(1),r(2),r(3),10);
figure(5)
hold on 
plot3(v(:,1),v(:,2),v(:,3),'.k')
plot3(X(:),Y(:),Z(:),'.g')
daspect([1,1,1])
hold off

%-construct the projectors
%--------------------------------------------------------------------------
a = max(r); %major
b = min(r); %minor
%use function output to run spm_ipharm(vtest,n,a,b,Lin);
end

%%%%%%%%
function [ o, r] = spheroid_fit( X, ax )
%Tierney et al https://github.com/spm/spm/blob/main/spm_opm_amm.m
%  [o, r] = ellipsoid_fit( X, ax);
%
% Parameters:
%  X   - Coordinates  n x 3 matrix
%  ax  - numeric indicating longer axis
%
% Output:
% o   - origin
% r   - radii


x =X(:,1);
y =X(:,2);
z =X(:,3);
on = ones(size(x,1),1);
b = x.^2 + y.^2 + z.^2;
if ax==1
    A = [ y.^2 + z.^2 - 2*x.^2, 2*x,2*y,2*z,on];
    beta = pinv(A) *b;
    v(1) = -2 * beta(1) - 1;
    v(2) = beta(1) - 1;
    v(3) = beta(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

if ax==2
    A = [ x.^2 + z.^2 - 2*y.^2, 2*x,2*y,2*z,on];
    beta = pinv(A)*b;
     v(1) = beta(1)-1;
    v(2) = -2*beta(1)-1;
    v(3) = beta(1)-1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

if ax==3
    A = [ x.^2 + y.^2 - 2*z.^2, 2*x,2*y,2*z,on];
    beta = pinv(A) *b; 
    v = beta;
    v(1) = beta(1) - 1;
    v(2) = beta(1) - 1;
    v(3) = -2*beta(1) - 1;
    v = [ v(1) v(2) v(3) 0 0 0 beta( 2 : 5 )' ];
end

A = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
    

o = -A( 1:3, 1:3 ) \ v( 7:9 )';
T = eye( 4 );
T( 4, 1:3 ) = o';
R = T * A * T';
[ vec, s ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
r = sqrt( 1 ./ diag( abs( s ) ) );
sgns = sign( diag( s ) );
r = r .* sgns;
r =vec*r;

end
