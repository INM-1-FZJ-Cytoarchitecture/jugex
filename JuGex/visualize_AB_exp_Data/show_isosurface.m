function data=show_isosurface(data)
%figure
%cla(data.handles.axes3d);
data.isorender=true;
data.handle_patch= patch(data.FV,'FaceLighting', 'gouraud', 'Facecolor', data.color,'Tag','brain_surface','EdgeColor', 'none','VertexNormals',data.N);
material dull
grid on
%data.handel_light=camlight('headlight');
axis equal; 
daspect(data.scales); 
%axis off; 
