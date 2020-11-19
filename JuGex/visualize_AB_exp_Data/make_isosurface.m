% Make the vertex, face and normals list of the iso-surface
function data=make_isosurface(data)

if(data.scaling==0.25)
    data.FV = isosurface(data.small1, data.iso);
    data.N  = isonormals(smooth3(data.small1),data.FV.vertices);
elseif(data.scaling==0.5)
    data.FV = isosurface(data.small2, data.iso);
    data.N  = isonormals(smooth3(data.small2),data.FV.vertices);
else
    data.FV = isosurface(data.volume, data.iso);
    data.N  = isonormals(smooth3(data.volume),data.FV.vertices);
end
data.FV.vertices=data.FV.vertices*(1/data.scaling);
data.FVmean=mean(data.FV.vertices,1);
