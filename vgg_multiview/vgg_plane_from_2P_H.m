% A = vgg_plane_from_2P_H(P1,P2,H)  From two cameras and inter-image homography, it gets the scene plane.
%
% The method is linear and there is no preconditioning - works only for consistent triplet (P1,P2,H).
%
% P1, P2 ... double (3,4)
% H ... double (3,3)
% A ... double (4,1)
%
% It is H = P2*vgg_H_from_2P_plane(A,P1). Try the example at the end of this file.

function A = vgg_plane_from_2P_H(P,P2,H)

Q = [P2*vgg_contreps(P(2,:)'*P(3,:) - P(3,:)'*P(2,:))
     P2*vgg_contreps(P(3,:)'*P(1,:) - P(1,:)'*P(3,:))
     P2*vgg_contreps(P(1,:)'*P(2,:) - P(2,:)'*P(1,:))];
A = (Q\H(:))';

return


P = randn(3,4);
P2 = randn(3,4);
A = randn(1,4);
H = P2*vgg_H_from_P_plane(A,P);
vgg_plane_from_2P_H(P,P2,H) ./ A

