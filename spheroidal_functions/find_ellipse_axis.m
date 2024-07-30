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

if (ind~=2)
    error('Y is not longest axis.... fix please')
end


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

% [X,Y,Z]=ellipsoid(o(1),o(2),o(3),r(1),r(2),r(3),10);
% figure(5)
% hold on 
% plot3(v(:,1),v(:,2),v(:,3),'.k')
% plot3(X(:),Y(:),Z(:),'.g')
% daspect([1,1,1])
% hold off

%-construct the projectors
%--------------------------------------------------------------------------
a = max(r); %major
b = min(r); %minor
%use function output to run spm_ipharm(vtest,n,a,b,Lin);
end
