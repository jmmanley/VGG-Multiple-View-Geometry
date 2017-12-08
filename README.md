# MATLAB Functions for Multiple View Geometry

Obtained from [http://www.robots.ox.ac.uk/~vgg/hzbook/code/](http://www.robots.ox.ac.uk/~vgg/hzbook/code/).

Please report any bugs to [Andrew Zisserman](removethisifyouarehuman-az@robots.ox.ac.uk).


### Acknowledgements:
These functions are written by: David Capel, Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler, and Andrew Zisserman



## VGG MultiView Compute Library

### Conversions:
* `vgg_KR_from_P.m` extract K, R from P such that P = K*R*[eye(3) -t]
* `vgg_F_from_P.m` fundamental matrix from 2 cameras
* `vgg_P_from_F.m` 2 camera matrices from fundamental matrix
* `vgg_T_from_P.m` trifocal tensor from 3 cameras
* `vgg_H_from_2P_plane.m` inter-image homography from 2 cameras and 3D plane
* `vgg_H_from_P_plane.m` projection matrix from image onto 3D plane
* `vgg_plane_from_2P_H.m` 3D plane from 2 cameras and inter-image homography

### Multiview tensors from image correspondences:
* `vgg_H_from_x_lin.m` homography from points in 2 images, linear method
* `vgg_H_from_x_nonlin.m` MLE of the above, by nonlinear method
* `vgg_Haffine_from_x_MLE.m` MLE of affine transformation from points in 2 images, linear
* `vgg_F_from_7pts_2img.m` fundamental matrix from 7 points in 2 images
* `vgg_PX_from_6pts_3img.m` cameras and world points from 6 points in 3 images

### Preconditioning for estimation:
* `vgg_conditioner_from_image.m` conditioning shift+scaling from image dimensions
* `vgg_conditioner_from_pts.m` conditioning shift+scaling from image points

### Self-calibration and similar:
* `vgg_signsPX_from_x.m` swaps signs of P and X so that projection scales are positive
* `vgg_selfcalib_qaffine.m` quasi-affine from projective reconstruction
* `vgg_selfcalib_metric_vansq.m` metric from projective and 3 orthogonal principal directions and square pixels

### Estimation:
* `vgg_X_from_xP_lin.m` 3D point from image projections and cameras, linear
* `vgg_X_from_xP_nonlin.m` MLE of that, non-linear method
* `vgg_line3d_from_lP_lin.m` 3D line segment from image line segments and cameras, linear
* `vgg_line3d_from_lP_nonlin.m` MLE of that, non-linear method



## VGG User Interface Library

### GUI’s:
* `vgg_gui_F.m` Visualizes epipolar geometry between two views
* `vgg_gui_H.m` Visualizes a homography between two views



## Examples

These examples use images and matrices included in the directory `vgg_examples`. Change to that directory before running the example functions.

* `view_homog_ex.m` Example of using `vgg_gui_H`
* `view_fund_ex.m` Example of using `vgg_gui_F`
* `Haffine_from_x_MLE_ex.m` Example of using `vgg_Haffine_from_x_MLE`
* `F_from_Ps_ex.m` Example on computing F from two camera matrices using `vgg_F_from_P`
* `H_from_image_corr_ex.m` Example on computing H from points using `vgg_H_from_x_lin`
* `testhomog_vgg.m` Example of computing H from two images from a rotating camera. This example also requires `ransacfithomography_vgg.m` and Peter Kovesi's functions (such as `matchbycorrelation.m` and `ransac.m`). See link below.



## Links to other highly recommended Computer Vision software

* [Peter Kovesi's Matlab Functions for Computer Vision and Image Analysis](http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/index.html)
* [Jean-Yves Bouguet's Matlab Calibration Software](http://www.vision.caltech.edu/bouguetj/calib_doc/)



## Note on release versions

November 2012 updates to

* `ransacfithomography_vgg.m`
* `testhomog_vgg.m`
* `vgg_H_from_x_nonlin.m`

To maintain compatibility with Peter Kovesi's functions and for Matlab R2012a compatibility.

Thanks to: Relja Arandjelovic, Peter Corke and Alexander Khanin.