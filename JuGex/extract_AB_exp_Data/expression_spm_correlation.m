% Get a list of probes that are differentially 
% expressed between contrast vs. target structures 
function [samples,probes,zscores_of_validated_corrs,data_2_plot] = expression_spm_correlation(activation_file, ~, ~, specimen, gene_list,map_threshold,search_mode,struct_id,disp_url,disp_detail,parameter)

display_TB_details=disp_detail;


% Copyright 2013 Allen Institute for Brain Science
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

[ids,gene_symbol,entrez_id] = read_gen_file(gene_list);
%%%

[samples, explevels,zscores, probes] = download_expression(ids, [specimen.donor_id],search_mode,struct_id,disp_url);  %% An der Stelle die GenIDs aus Schizo und Depr. nehmen
%%% Aufbau samples  
%%% ----------------------------------------------------------
%(Wichtig: Enthält Ortsinfos der Samples) 
%( size(samples,2) = Anzahl der Samples an der die gesuchte ID untersucht
% wurde )
%
% samples{1,1}
%                   donor: [1x1 struct]
%                  sample: [1x1 struct]
%               structure: [1x1 struct]
%     top_level_structure: [1x1 struct]
% samples{1,1}.sample
%                            well: 4668
%                         polygon: 993106
%                             mri: {[126]  [79]  [57]}


%%% Aufbau probes  
%%% ----------------------------------------------------------
%(Wichtig: Enthält Exp Level der Probes) 
%( size(explevels) = Anzahl der Samples X Anzahl der Probe_ids


%%% Aufbau explevels  
%%% ----------------------------------------------------------
%   /Redundante Rückgabe aus Funktion. Zur Überprüfung der Zuordnung von
%   Exp Level zur Probe_Id

if search_mode==2
    corrs='';
    zscores_of_validated_corrs=zscores;
    spmhd = spm_vol(activation_file);
    spmMNI = spmhd.mat;
    % Concatenate the transform from MNI to SPM image onto original MNI
    % transform. Now we can transform a T1 coordinate directly to an SPM voxel.
    aibsToSPM = inv(spmMNI)*specimen.alignment3d;

    coords = transform_samples(samples,aibsToSPM);
    coords = int32(round(coords));  % Koordinaten triplet für alle samples
    
    data_2_plot={coords,zscores,explevels};
    
    if display_TB_details==1
    for i=1:size(samples,2)
    fprintf('Donor: %s\tMap: ---\tStructure: %s\tWell_ID: %d\tMRI: %d %d %d\n',samples{1,i}.donor.name,samples{1,i}.structure.name,samples{1,i}.sample.well,samples{1,i}.sample.mri{1,1},samples{1,i}.sample.mri{1,2},samples{1,i}.sample.mri{1,3});
    end
    end
else

[zscores_of_validated_corrs,data_2_plot] = filter_sample_positions(activation_file, samples, explevels,zscores, specimen.alignment3d, map_threshold,display_TB_details,parameter);
end

% Compute the correlation between expression levels and spm voxels.
function [zscores,data_2_plot] = filter_sample_positions(activation_file,samples, explevels, zscores, MNI,map_threshold,display_TB_details,parameter)

% Load the activation and mask images.   
spmhd = spm_vol(activation_file);  %%% evtl. Header nicht kompatibel????
spmimg = spm_read_vols(spmhd);
spmMNI = spmhd.mat;

% Concatenate the transform from MNI to SPM image onto original MNI
% transform. Now we can transform a T1 coordinate directly to an SPM voxel.
aibsToSPM = inv(spmMNI)*MNI;

coords = transform_samples(samples,aibsToSPM);
coords = int32(round(coords));  % Koordinaten triplet für alle samples

% Find samples inside the image.
for i=1:size(coords,2)
    coord = coords(:,i);   
    if ((sum(coord>0) ~= 3) || spmimg(coord(1),coord(2),coord(3)) <= str2num(map_threshold)/10 || spmimg(coord(1),coord(2),coord(3)) ==0)    
        coords(:,i) = -1;
    end
end

% Filter the exp levels and coordinates to only include the valid ones.
% extract well ids
try
    sf_config = jsondecode(fileread('sf.json'));
    display_TB_details=str2num(sf_config.display_detail);
    
    disp('sf.json file in pwd detected');
    disp('continue with following parameter');
    disp(' ');
    disp(['Display mode: ' sf_config.display_detail]);
    disp(' ');
    disp(['SF active for VOI1: ' sf_config.activate_voi_1_sf]);
    disp(['VOI1 SF: ' sf_config.voi_1_sf]);
    disp(['SF active for VOI2: ' sf_config.activate_voi_2_sf]);
    disp(['VOI2 SF: ' sf_config.voi_2_sf]);  
    disp('press any key to continue');
    pause()
    % prüfen ob gerade VOI1 oder VOI2 ausgewertet wird
    if strcmp(activation_file,parameter.pmap{1,1}) %VOI1 wird gerade bearbeitet
        if strcmp(sf_config.activate_voi_1_sf,'1')%prüfe ob SF für erste VOI aktiviert werden soll
            filter_strings=sf_config.voi_1_sf;% setzte SF auf Wert aus JSON
        else %filter soll für VOI1 nicht aktiviert werden
            filter_strings={''};
        end
    end
    % prüfen ob gerade VOI1 oder VOI2 ausgewertet wird
    
    if strcmp(activation_file,parameter.pmap{1,2})%VOI2 wird gerade bearbeitet
        if strcmp(sf_config.activate_voi_2_sf,'1')%prüfe ob SF für zweites VOI aktiviert werden soll
            filter_strings=sf_config.voi_2_sf;% setzte SF auf Wert aus JSON
        else %filter soll für VOI2 nicht aktiviert werden
            filter_strings={''};
        end
    end
    
    
catch
    disp('No sf.json file in pwd detected. Continue without sematic filtering');
    pause();
    filter_strings={''};
end




for i=1:size(samples,2)
    well_id(i)=samples{1,i}.sample.well;
    structure_name{i}=samples{1,i}.structure.name;
    donor_name{i}=samples{1,i}.donor.name;
    mri_raw(i,1:3)=[samples{1,i}.sample.mri{1,1},samples{1,i}.sample.mri{1,2},samples{1,i}.sample.mri{1,3}];
    
    % include semantic Filter => here hardcoded 'insu'  should be cell
    % array in order to allow for more than user chosen structure
    %filter_strings={'anterior orbital gyrus','lateral orbital gyrus','inferior frontal gyrus, orbital part'};
    %filter_strings={''};
    
    %%%%%%%%%%%filter string muss ein cell sein!!!!
    
    sf_true(i)=0;
    for check_filter_strings=1:size(filter_strings,2)
        if isempty(filter_strings{1})
        sf_true(i)=1;
        else
            if ~isempty(strfind(samples{1,i}.structure.name,filter_strings{check_filter_strings})) && sf_true(i)==0
                sf_true(i)=1;
                %continue;
            end
        end
    end
    %sprintf('Donor: %s   Structure: %s    Well_ID: %d      MRI: %d %d %d',samples{1,i}.donor.name,samples{1,i}.structure.name,samples{1,i}.sample.well,samples{1,i}.sample.mri{1,1},samples{1,i}.sample.mri{1,2},samples{1,i}.sample.mri{1,3})
end


valid = sum(coords>0)==3;
% include semantic filter => pairwise multiply valid with sf_true to eliminate samples which are inside the map but not on semantic filter structure 
sprintf('%03d orig filtered TBs; %03d sf filtered TBs   => %03d remaining TBs',sum(valid),sum(sf_true),sum(valid.*sf_true))
valid=logical(valid.*sf_true);
well_id=well_id';
well_id=well_id(valid,:);
structure_name=structure_name';
structure_name=structure_name(valid,:);
donor_name=donor_name';
donor_name=donor_name(valid,:);

mri_raw=mri_raw(valid,:);


explevels = explevels(valid,:);
zscores = zscores(valid,:);
coords = coords(:,valid);
nsamples = size(coords,2);
data_2_plot={coords,zscores,explevels,well_id};
[~,name,~] = fileparts(activation_file);  % nicht geprüft ob es geht, wenn zwei maps genommen werden. ausgabe geprüft mit map gegen aba macro label
if display_TB_details==1 
for t=1:size(explevels,1)
fprintf('Donor: %s\tMap: %s\tStructure: %s\tWell_ID: %d\tMRI: %d %d %d\n',donor_name{t},name,structure_name{t},well_id(t),mri_raw(t,:))
end
end

