function X = vgg_intersect_quadrics(A, B, C);
%
% function X = vgg_intersect_quadrics(A, B, C);
%
% Purpose:
%   Compute intersections of three quadrics in projective 3-space.
%
% Input:
%   A, B, C are symmetric 4x4 matrices.
%
% Output:
%   X is 4-by-8 complex. Each column is a root.
%
% The function works by converting the problem into an 8x8 generalized
% eigen-problem.
%
% fsm@robots.ox.ac.uk
%

if any(size(A) ~= [4 4])
  error('bad A')
end
if any(size(B) ~= [4 4])
  error('bad B')
end
if any(size(C) ~= [4 4])
  error('bad C')
end

% symmetrize and normalize.
A = A + A.'; A = A / norm(A, 'fro'); a = [A(1, 1) 2*A(1, 2) 2*A(1, 3) 2*A(1, 4) A(2, 2) 2*A(2, 3) 2*A(2, 4) A(3, 3) 2*A(3, 4) A(4, 4)].';
B = B + B.'; B = B / norm(B, 'fro'); b = [B(1, 1) 2*B(1, 2) 2*B(1, 3) 2*B(1, 4) B(2, 2) 2*B(2, 3) 2*B(2, 4) B(3, 3) 2*B(3, 4) B(4, 4)].';
C = C + C.'; C = C / norm(C, 'fro'); c = [C(1, 1) 2*C(1, 2) 2*C(1, 3) 2*C(1, 4) C(2, 2) 2*C(2, 3) 2*C(2, 4) C(3, 3) 2*C(3, 4) C(4, 4)].';

% make multiplication tables (for degrees 1+2=3, 2+2=4 and 3+1=4)
if 1
  ijk = [
    1
    22
    43
    64
    82
    105
    126
    147
    163
    186
    208
    229
    244
    267
    289
    310
    325
    351
    372
    393
    406
    432
    454
    475
    487
    513
    535
    556
    568
    594
    617
    638
    649
    675
    698
    719
    730
    756
    779
    800
    ];
  X12 = zeros(20, 4, 10); X12(ijk) = 1;
  
  ijk = [
    1
    37
    73
    109
    145
    181
    217
    253
    289
    325
    352
    390
    426
    462
    501
    537
    573
    609
    645
    681
    703
    741
    778
    814
    852
    889
    925
    962
    998
    1034
    1054
    1092
    1129
    1165
    1203
    1240
    1276
    1313
    1349
    1385
    1405
    1446
    1482
    1518
    1561
    1597
    1633
    1669
    1705
    1741
    1756
    1797
    1834
    1870
    1912
    1949
    1985
    2022
    2058
    2094
    2107
    2148
    2185
    2221
    2263
    2300
    2336
    2373
    2409
    2445
    2458
    2499
    2537
    2573
    2614
    2652
    2688
    2726
    2762
    2798
    2809
    2850
    2888
    2924
    2965
    3003
    3039
    3077
    3113
    3149
    3160
    3201
    3239
    3275
    3316
    3354
    3390
    3428
    3464
    3500
    ];
  X22 = zeros(35, 10, 10); X22(ijk) = 1;
  
  ijk = [
    1
    37
    73
    109
    145
    181
    217
    253
    289
    325
    361
    397
    433
    469
    505
    541
    577
    613
    649
    685
    702
    740
    776
    812
    851
    887
    923
    959
    995
    1031
    1071
    1107
    1143
    1179
    1215
    1251
    1287
    1323
    1359
    1395
    1403
    1441
    1478
    1514
    1552
    1589
    1625
    1662
    1698
    1734
    1772
    1809
    1845
    1882
    1918
    1954
    1991
    2027
    2063
    2099
    2104
    2142
    2179
    2215
    2253
    2290
    2326
    2363
    2399
    2435
    2473
    2510
    2546
    2583
    2619
    2655
    2692
    2728
    2764
    2800
    ];
  X31 = zeros(35, 20, 4); X31(ijk) = 1;
  
  clear ijk
end

% compose ideal in degree 2.
I2 = [a b c];

% compute ideal in degrees 3 and 4.
I3 = reshape(reshape(X12, [20* 4 10]) * I2, [20 3* 4]);
I4 = reshape(reshape(X22, [35*10 10]) * I2, [35 3*10]);

% the columns of C3 and C4 span complements of I3 and I4, respectively.
[U,S,V] = svd(I3.'); C3 = V(:, 13:20); %W3 = diag(S).'
[U,S,V] = svd(I4.'); C4 = V(:, 28:35); %W4 = diag(S).'

% make the four multiplication operators.
M = cell(1, 4);
for k=1:4
  M{k} = C4.' * X31(:, :, k) * C3;
end

% FIXME: here we only use M{1} and M{2} to compute the Schur
% decomposition although we use all four multiplication operators to
% extract the root coordinates.
[AA, BB, Q, Z] = qz(M{1}, M{2});
for k=1:4
  M{k} = Q*M{k}*Z;
end
%M{:}

X = zeros(4, 8);
X(1, :) = diag(M{1}).';
X(2, :) = diag(M{2}).';
X(3, :) = diag(M{3}).';
X(4, :) = diag(M{4}).';

% normalize
X = X ./ repmat(max(abs(X), [], 1), [4 1]);

return;

