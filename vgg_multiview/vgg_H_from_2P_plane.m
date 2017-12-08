function [H] = vgg_H_from_2P_plane(P1,P2, L)
%
%	[H] = vgg_H_from_2P_plane(P, L)
%	[H] = vgg_H_from_2P_plane(P1, P2, L)
%
%USAGE 1:
% Returns the homography matrix from the plane whos normal is L to the
% image whos projection matrix is P. This is useful for mapping two images
% of a plane.
%USAGE 2:
% Given two camera matrices and a plane, returns the homography matrix that
% maps points from the first camera onto the second.
%
%IN:
%	P,P1,P2 - 3x4 Camera matrix
%	L - 1x4 Plane normal
%
%OUT:
%	H - 3x3 projective homography matrix mapping points from the image
%	to the world-points on the plane. Note that 
%
%EXAMPLE:
%	Let P1, P2 be two projection matrices and let L be the normal to
%	the plane then the homography H=H1*inv(H2) maps points in P2 to
%	points in P1 where:
%		H1=vgg_H_from_2P_plane(P1, L)
%		H2=vgg_H_from_2P_plane(P2, L)
%		% And this will be the resulting image:
%		[u,v]=homflow(H, size(i1,1), size(i1,2));
%		t=imgwarp(i1, u,v);
%
% An alternative is to use H = P1*vgg_H_from_P_plane(P2,L).

% $Id: vgg_H_from_2P_plane.m,v 1.2 2002/02/22 22:30:55 werner Exp $
% Yoni, Fri Apr  6 12:44:54 2001

if nargin==2
   P=P1; L=P2;
   
   H = P(:,1:3)*L(4) - P(:,4)*L(1:3);

elseif nargin==3
   nL=null(L);
   H=P2*nL*inv(P1*nL);

else
   error('Wrong number of arguments');
end
