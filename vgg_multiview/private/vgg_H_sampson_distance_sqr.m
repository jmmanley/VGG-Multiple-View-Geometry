function d = vgg_H_sampson_distance_sqr(H,X1,X2)

% d = vgg_H_sampson_distance_sqr(H,X1,X2)
%
% Returns an approximation of the squared geometric distance from the
% 4d joint-space point [X1;X2] to the H manifold. See Torr CVIU 97.

if (size(X1) ~= size(X2))
  error('Point sets not same size!');
end

p1 = X1 ./ repmat(X1(3,:),3,1);
p2 = X2 ./ repmat(X2(3,:),3,1);

alg = vgg_H_algebraic_distance(H,p1,p2);

N = size(X1,2);

h = reshape(H',9,1);

G1 = [ H(1,1) - p2(1,:) * H(3,1) ; ...
  H(1,2) - p2(1,:) * H(3,2) ; ...      
  -p1(1,:) * H(3,1) - p1(2,:) * H(3,2) - H(3,3) ; ...
  zeros(1,N) ];

G2 = [ H(2,1) - p2(2,:) * H(3,1) ; ...
  H(2,2) - p2(2,:) * H(3,2) ; ...
  zeros(1,N) ; ...
  -p1(1,:) * H(3,1) - p1(2,:) * H(3,2) - H(3,3) ];

magG1 = sqrt(sum(G1 .* G1));
magG2 = sqrt(sum(G2 .* G2));
magG1G2 = sum(G1 .*  G2);

alpha = acos( magG1G2 ./ (magG1 .* magG2) );

D1 = alg(1,:) ./ magG1;
D2 = alg(2,:) ./ magG2;

d = (D1.*D1 + D2.*D2 - 2 * D1 .* D2 .* cos(alpha)) ./ sin(alpha);


