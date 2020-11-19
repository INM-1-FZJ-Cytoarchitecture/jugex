function spm_merge_maps(input_maps_cell,save_str)
% spm_merge_maps  Merge two or more maps into one file.
% !!!!!!!!!! only merge maps if they are in the same space !!!!!!
% spm_vol and spm_read_vol is ussed to read map files
%
%   spm_merge_maps(input_maps_cell,save_str)
%
%   input_maps_cell:    cell array including file names of the maps. 
%   save_str:           filename of merged file. Will be saved in pwd.
%                       Saves as .nii or img/hdr file depending on extension 
%
%   Example:
%   maps_2_combine={'map1.nii','map2.nii','map3.nii'};
%   spm_merge_maps(maps_2_combine,'merged_maps.nii');
% 
%   Will save the maps 'map1.nii','map2.nii','map3.nii' merged in pwd as merged_maps.nii nifti file

    

V1=spm_vol(input_maps_cell{1});
V2=spm_vol(input_maps_cell{2});
Y1=spm_read_vols(V1);
Y2=spm_read_vols(V2);
Y3=Y1+Y2; % only ever do this if they're in the same space
if size(input_maps_cell,2)>2  % merge more than 2 maps
for i=3:size(input_maps_cell,2) % start with index 3 since first 2 are allready merged in Y3
    V_add=spm_vol(input_maps_cell{i});
    Y_add=spm_read_vols(V_add);
    Y3=Y3+Y_add; % only ever do this if they're in the same space
end
end
% copy first header to combined header. Possible since both are in the same
% space
V3=V1;
V3.fname=save_str;
spm_write_vol(V3,Y3);


end