
%%This function shows how to use MagTense for a grid of PM cubes where some
%%cubes have a large permeability (modeling Fe with a constant relative
%%permeability)
function [] = MagTense_Example002_PM_Fe_const_mur()

%make sure to source the right path for the generic Matlab routines
addpath(genpath('../util/'));
addpath('../../MEX_files/');
%define the vacuum permeability
mu0 = 4*pi*1e-7;

mur = 3;
Ms = 2.1;

%%Get a default tile from MagTense
tile = getDefaultMagTile();
    
%ensure the tile is a permanent magnet
tile.magnetType = getMagnetType('hard');

%set the geometry to be a rectangular prism
tile.tileType = getMagTileType('prism');

%set the dimensions of the prism (5 mm on either side)
d = 0.1;
tile.abc = [d,d,d];

%set the center position of the prism (centered at Origin)
tile.offset = [0.45,0.55,d/2];

%set the easy axis of the cube. This is expected to be with respect to the global
%coordinate system basis
tile.u_ea = [0,1,0];
%ensure the two hard axes are perpendicular and normalized
tile.u_oa1 = [1,0,0];
tile.u_oa2 = [0,0,1];

%set the relative permeability for the easy axis (NdFeB magnet)
tile.mu_r_ea = 1; %1.06;
%and for the two hard axes
tile.mu_r_oa = 1; %1.17;

%set the remanence of the magnet (1.2 T converted to A/m)
tile.Mrem = 1.2 / mu0;

%iron
iron = tile;
iron.color = [0,0,1];
iron.magnetType = getMagnetType('soft');
iron.offset = [0.65,0.55,d/2];
%set the relative permeability for the easy axis (iron)
iron.mu_r_ea = mur;
%and for the two hard axes
iron.mu_r_oa = mur;
iron.Mrem = 0; 
tiles = [tile, refineTiles( iron, res )];

% %build a 2D grid of 5x5 identical cubes
% nx = 10;
% ny = 10;
% 
% %make space for the struct array
% tiles(1) = tile;
% tiles(nx*ny) = tile;

% for i=1:nx
%     for j=1:ny
%         ind = (i-1)*ny+j;
%         tiles(ind) = tile;
%         %set the center of the magnet at i,j
%         tiles(ind).offset(1) = -floor(nx/2) * d + (i-1) * d;
%         tiles(ind).offset(2) = -floor(ny/2) * d + (j-1) * d;
%         
%         [ph,th,r] = cart2sph(tiles(ind).offset(1),tiles(ind).offset(2),tiles(ind).offset(3));
%         tiles(ind).u_ea(1) = cos(ph) * cos(ph) - sin(ph) * sin(ph);
%         tiles(ind).u_ea(2) = cos(ph) * sin(ph) + sin(ph) * cos(ph);
%         
%         tiles(ind).u_oa1(1) = -tiles(ind).u_ea(2);
%         tiles(ind).u_oa1(2) = tiles(ind).u_ea(1);
%     end
% end
% n_r = 5;
% ind_lst = zeros(n_r,1);
% tiles_tmp = tiles;
% res = struct('nx',2,'ny',2,'nz',1);
% for i=1:n_r
%     %insert Fe in n (random) tiles
%     fnd = false;
%     while fnd == false
%         ind = randi(nx*ny);
%         if ~any( ind-ind_lst == 0 )
%             fnd = true;
%             ind_lst(i) = ind;
%         end
%     end  
%     tiles_tmp(ind).color = [0,0,1];
%     
%     tiles_tmp(ind).magnetType = getMagnetType('soft');
%     
%     
%     tiles = [tiles, refineTiles( tiles_tmp(ind), res )];
% end
% 
% tiles(ind_lst) = [];

%%Let MagTense find the magnetization vector of the cube by iterating to a
%%self-consistent solution. The two empty arguments ( [] ) ensure that
%%default values are used and they are in any case not relevant for this
%%example. The 1e-6 argument is the maximum relative error (change in
%%magnetization from one iteration to the next) allowed before convergence
%%is reached. Finally, 100 reflects the max. no. of allowed iterations

%load the const mur state function
stFcn=MakeMH_Fe_const_mur( mur, Ms );

tiles = IterateMagnetization( tiles, stFcn, 300, 1e-3, 100 );
%tiles = IterateMagnetization( tiles, [], [], 1e-3, 100 );

close all;%plotTiles(tiles,true);axis equal;alpha 0.3
%%Now find the field in a set of points
%define a range of points spanning the xy plane at z=0
x = linspace( 0,1, 101);
y = linspace( 0,1, 101);
z = d/2;

%use meshgrid to fill out the span
[X,Y,Z] = meshgrid(x,y,z);


%get the field
[H,Hnorm] = getHMagTense( tiles, X, Y, Z );
H = H .* mu0;
Hnorm = Hnorm .* mu0;
%Plot the solution in a contour map
surf_and_con( X,Y,Hnorm);
xlabel('x');
ylabel('y');

%%plot the magnetic field (H). Notice how the field opposes the
%%magnetization inside the cube.
getFigure();
quiver(X,Y,H(:,:,1),H(:,:,2),2,'linewidth',2);
contour(X,Y,Hnorm,'linewidth',1);
plotTiles(tiles,true);
alpha 0.3;
axis equal;
xlabel('x');
ylabel('y');

%%Now find the field in a set of points
x = 0:0.001:1;
y = 0:0.001:1;
z = -0.5:0.001:0.5;

offset = [0.45,0.55,d/2];

pts = zeros( numel(x)+numel(y)+numel(z), 3 );
pts(1:numel(x),:) = [x; zeros(1,numel(x))+offset(2); zeros(1,numel(x))+offset(3)]';
pts((numel(x)+1):(numel(x)+numel(y)),:) = [zeros(1,numel(y))+offset(1); y; zeros(1,numel(y))+offset(3)]';
pts((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z)),:) = [zeros(1,numel(z))+offset(1); zeros(1,numel(z))+offset(2); z]';

%get the field from the Fortran implementation
tic
H = getHFromTiles_mex( tiles, pts, int32(length(tiles)), int32(length(pts(:,1))) );
toc

%Find the norm of the field
Hnorm = squeeze( sqrt( sum(H.^2,2) ) );

%Plot the solution
getFigure();
plot(x,Hnorm(1:numel(x)),'r');
getFigure();
plot(y,Hnorm((numel(x)+1):(numel(x)+numel(y))),'g');
getFigure();
plot(z,Hnorm((numel(x)+1+numel(y)):(numel(x)+numel(y)+numel(z))),'b');

% %define a range of points spanning the xy plane at z=0
% x = linspace( -0.05,0.05, 20);
% y = 0.025 + 5e-3; %height above the cube
% z = linspace( -0.06,0.06, 21);
% 
% %use meshgrid to fill out the span
% [X,Y,Z] = meshgrid(x,y,z);
% 
% 
% %get the field
% [H,Hnorm] = getHMagTense( tiles, X, Y, Z );
% H = H .* mu0;
% Hnorm = Hnorm .* mu0;
% 
% %Plot the solution in a contour map
% surf_and_con( squeeze(X),squeeze(Z),Hnorm);
% xlabel('x');
% ylabel('z');
% 
% 
% getFigure();
% quiver(squeeze(X),squeeze(Z),H(:,:,1),H(:,:,3),2,'linewidth',2);
% contour(squeeze(X),squeeze(Z),Hnorm,'linewidth',2);
% 
% alpha 0.3;
% axis equal;
% xlabel('x');
% ylabel('z');
% 
% getFigure();contour(squeeze(X),squeeze(Z),H(:,:,1),'linewidth',2,'showtext','on');
% xlabel('x');
% ylabel('z');
end