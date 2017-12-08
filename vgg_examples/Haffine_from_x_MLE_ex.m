% Example of using vgg_Haffine_from_x_MLE

% data is square, rotated square

exem1 = [ 0 1 1 0 ; 0 0 1 1 ];
exem2 = [ 0.5 1 0.5 0; 0 0.5 1 0.5 ];

H = vgg_Haffine_from_x_MLE(vgg_get_homg(exem1) ,vgg_get_homg(exem2));
disp('Affine H relating two noise free point sets:');
disp(H);

% add noise

dat = 0.05 * randn(2,4);
exem1 = exem1 + dat;

dat = 0.05 * randn(2,4);
exem2 = exem2 + dat;

H = vgg_Haffine_from_x_MLE(vgg_get_homg(exem1) ,vgg_get_homg(exem2));
disp(' ')
disp('MLE estimate of affine H relating the point sets with added noise:');
disp(H);
