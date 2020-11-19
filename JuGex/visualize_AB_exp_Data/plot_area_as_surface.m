function patch_hdl=plot_area_as_surface(volume,color,alpha)
data.volume=volume;
hold on
data=struct;
data.volume=volume;
% Calculate min and max value to scale the intensities
data.vmax=double(max(data.volume(:)));
data.vmin=double(min(data.volume(:)));
% Get or Calculate the staritng iso value
data.iso=double(0.25);  % set iso to 0.2
% Default input values
data.scaling=1;
data.posx=0.5;
data.posy=0.5;
data.posz=0.5;
data.scales=[1 1 1];
data.color=[.95 .87 .73];
data.sizes=size(data.volume);
data.isorender=false;
data.slicerender=false;
data.h1=[];
data.h2=[];
data.h3=[];
% make downsampled volumes, for fast iso surface rendering.
data.small1=imresize3d(data.volume,0.25,[],'linear');
data.small2=imresize3d(data.volume,0.5,[],'linear');
% create isosurface
data=make_isosurface(data);
% plot isosurface to figure
data=show_isosurface(data);
set(data.handle_patch,'Facecolor',color);
set(data.handle_patch,'Facealpha',alpha);
patch_hdl=data.handle_patch;