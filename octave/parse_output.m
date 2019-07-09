clear;
[var, _, val] = textread('input.txt', '%s %s %f', 'headerlines', 1);
grid_size = val(1);
dx = val(3);
num_grains = val(5);
num_steps = val(6);
n_step = val(7);
dt = val(8);
x_max = grid_size * dx;
y_max = grid_size * dx;
##num_steps = num_steps(end);
####n_step = num_steps;
[x,y] = meshgrid(linspace(0,grid_size - 1,grid_size),linspace(0,grid_size - 1,grid_size));
z=load('initial_structure.txt') - 1;
% Color matrix containing indices into the colormap. The values in C map colors
% in the colormap array to the vertices surrounding each face. The color of a face
% depends on the color at one of its four vertices. Of the four vertices, the one
% that come first in X and Y determines the color of the face. If you do not specify
% X and Y, MATLAB uses X=1:n and Y=1:m, where [m,n] = size(C). Because of this
% relationship between the vertex colors and face colors, none of the values
% i the last row and column of C are represented in the plot.
h = figure(1, 'position', [1200,500,350,350]);
colormap jet;
h1 = pcolor(x.*dx,y.*dx,z);
colorbar();
set(h1, 'EdgeColor', 'none');
title('Initial Condition');
axis 'square';
##saveas(h, 'initial_condition.jpg')
for i = n_step:n_step:num_steps
  z=load(['out_' int2str(i) '.txt']);
  h = figure();
  colormap gray;
  h1 = pcolor(x.*dx,y.*dx,z);
  caxis([0,1]);
  colorbar();
  set(h1, 'EdgeColor', 'none');
  title(sprintf('Timestep %d (%4.2fs)',i,dt*i));
  axis 'square';
%  saveas(h, sprintf('out_%d_aniso.jpg',i))
end
