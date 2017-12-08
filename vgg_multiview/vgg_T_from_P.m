% T = vgg_T_from_P(P)  Trifocal tensor from 3 camera matrices.
%
% P ... cell (3), camera matrices
% T ... double (3,3,3) 
%
% For 3 corresponding lines l1..l3 (each of size (1,3)) in cameras P1..P3 it is
%         for i=1:3, l1(1,i) = l2*T(:,:,i)*l3'; end
% up to scale.
%
% T is obtained by with unique absolute scale. In detail, for any cameras P{1:3} and 
% scene line L (4-by-4 Pluecker matrix, L+L'=0, vgg_contreps(L)*L=0), if
%
%   for i=1:3, l{i} = vgg_contreps(P{i}*L*P{i}'); end
%   T = vgg_T_from_P(P);
%   for i=1:3, l1(1,i) = l{2}*T(:,:,i)*l{3}'; end
%
% then
%
%   l1 ==  l{1}*(l{3}*e32) == -l{1}*(l{2}*e23)
%
% where e32 = P{3}*vgg_C_from_P(P{2}), e23 = P{2}*vgg_C_from_P(P{3}) are epipoles.

function T = vgg_T_from_P(P)

T(:,:,1) = P{2}*vgg_contreps(P{1}(2,:)'*P{1}(3,:)-P{1}(3,:)'*P{1}(2,:))*P{3}';
T(:,:,2) = P{2}*vgg_contreps(P{1}(3,:)'*P{1}(1,:)-P{1}(1,:)'*P{1}(3,:))*P{3}';
T(:,:,3) = P{2}*vgg_contreps(P{1}(1,:)'*P{1}(2,:)-P{1}(2,:)'*P{1}(1,:))*P{3}';

return