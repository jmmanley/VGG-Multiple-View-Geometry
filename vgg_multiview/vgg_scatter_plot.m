function vgg_scatter_plot(points,style)
%vgg_scatter_plot(points,style)
%
%Plots nonhomg points as a scatter plot

if nargin == 1, style = '+'; end
[r,c]=size(points);
if r == 2
  h = plot(points(1,:),points(2,:),style);
elseif r == 3
  h = plot3(points(1,:), points(2,:), points(3,:), style);
else
  error('rows != 2 or 3');
end
