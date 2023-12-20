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