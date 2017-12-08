function H=vgg_projective_basis_2d(p1,p2,p3,p4);
%
% function H=vgg_projective_basis_2d(p1,p2,p3,p4);
%
% Computes 3x3 homography matrix H such that 
%
% H * [p1 p2 p3 p4] ~ [ 1 0 0 1
%                       0 1 0 1
%                       0 0 1 1 ]
%
% fsm@robots
%

if ( (~ exist('p1')) | (~ exist('p2')) | (~ exist('p3')) | (~ exist('p4')) )
  error('missing input argument(s)');
end

if ( prod(size(p1)) == 3 )
  p1=reshape(p1,[3 1]);
else
  error('p1 must have 3 entries');
end

if ( prod(size(p2)) == 3 )
  p2=reshape(p2,[3 1]);
else
  error('p2 must have 3 entries');
end

if ( prod(size(p3)) == 3 )
  p3=reshape(p3,[3 1]);
else
  error('p3 must have 3 entries');
end

if ( prod(size(p4)) == 3 )
  p4=reshape(p4,[3 1]);
else
  error('p4 must have 3 entries');
end

%
tmp=[p1 p2 p3] \ p4;
p1 = p1*tmp(1);
p2 = p2*tmp(2);
p3 = p3*tmp(3);

% now p4=p1+p2+p3

H=inv([p1 p2 p3]);

return;
