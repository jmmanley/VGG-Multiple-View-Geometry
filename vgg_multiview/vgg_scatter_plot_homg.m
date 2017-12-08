function vgg_scatter_plot_homg(homg_points, style)
% vgg_scatter_plot_homg(homg_points, style)
%
% Scatter plot homogeneous points

if nargin == 1, style = '+'; end

vgg_scatter_plot(vgg_get_nonhomg(homg_points), style);
