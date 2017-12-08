% VGG MultiView Compute Library
%
% Conversions 
%   vgg_KR_from_P         - extract K, R from P such that P = K*R*[eye(3) -t]
%   vgg_F_from_P          - fundamental matrix from 2 cameras
%   vgg_P_from_F          - 2 camera matrices from fundamental matrix
%   vgg_T_from_P          - trifocal tensor from 3 cameras
%   vgg_H_from_2P_plane   - inter-image homography from 2 cameras and 3D plane
%   vgg_H_from_P_plane    - projection matrix from image onto 3D plane
%   vgg_plane_from_2P_H   - 3D plane from 2 cameras and inter-image homography
%
% Multiview tensors from image correspondences
%   vgg_H_from_x_lin             - homography from points in 2 images, linear method
%   vgg_H_from_x_nonlin          - MLE of the above, by nonlinear method
%   vgg_Haffine_from_x_MLE       - MLE of affine transformation from points in 2 images, linear
%   vgg_F_from_7pts_2img         - fundamental matrix from 7 points in 2 images
%   vgg_PX_from_6pts_3img        - cameras and world points from 6 points in 3 images
%
% Preconditioning for estimation
%   vgg_conditioner_from_image - conditioning shift+scaling from image dimensions
%   vgg_conditioner_from_pts   - conditioning shift+scaling from image points
%
% Self-calibration and similar
%   vgg_signsPX_from_x         - swaps signs of P and X so that projection scales are positive
%   vgg_selfcalib_qaffine      - quasi-affine from projective reconstruction
%   vgg_selfcalib_metric_vansq - metric from projective and 3 orthogonal principal directions and square pixels
%
% Estimation
%   vgg_X_from_xP_lin          - 3D point from image projections and cameras, linear
%   vgg_X_from_xP_nonlin       - MLE of that, non-linear method
%   vgg_line3d_from_lP_lin     - 3D line segment from image line segments and cameras, linear
%   vgg_line3d_from_lP_nonlin  - MLE of that, non-linear method
%
% 3D lines representations
%   vgg_line3d_pv_from_XY      - Pluecker vector from 2 points on the line
%   vgg_line3d_pv_from_pm      - Pluecker matrix from Pluecker vector
%   vgg_line3d_pm_from_pv      - Pluecker vector from Pluecker matrix
%   vgg_line3d_Ppv             - rearrange camera matrix to project Pluecker vector to image line
%   vgg_line3d_pv_from_2planes - Pluecker vector from 2 planes meeting in the line
%   vgg_line3d_XY_from_pm      - 2 points on 3D line from Pluecker matrix
%   vgg_line3d_XY_from_pv      - 2 points on 3D line from Pluecker vector
%  (vgg_contreps               - dual of Pluecker matrix of 3D line)
%
% Auxiliary & miscellaneous
%   vgg_get_homg    - adding row of ones
%   vgg_get_nonhomg - dividing by the final coordinates
%   vgg_projective_basis_2d
%   vgg_rms_error
%   vgg_scatter_plot_homg
%   vgg_scatter_plot
