function [main_r, TS]=pval_filter(main_r)
counter=1;

variable_names_types = [["area", "categorical"]; ...
            ["x","double"];...
            ["y","double"];...
            ["z","double"];...
            ["i","double"];...
            ["j","double"];...
            ["well_id","double"];...
            ["pval","double"]];
% Make table using fieldnames & value types from above
TS = table('Size',[0,size(variable_names_types,1)],... 
	'VariableNames', variable_names_types(:,1),...
	'VariableTypes', variable_names_types(:,2));




for i=1:size(main_r,2)
    %%%%hier muss der Fall abgefangen werden das gegen ABA Label gerechnet
    %%%%wird, dann ist main_r(i).pmap leer!!!
    if isempty(main_r(i).pmap)
        continue;
    else
        map=spm_vol(main_r(i).pmap);
        map_vol=spm_read_vols(map);
        for j=1:size(main_r(i).validated_zscores,1)
            
            x=main_r(i).data2plot{1, 1}(1,j);
            y=main_r(i).data2plot{1, 1}(2,j);
            z=main_r(i).data2plot{1, 1}(3,j);
            well_id=main_r(i).data2plot{1, 4}(j);     % ABA Label Extraktion hat keine well_id  muss eingebaut werden
            %sprintf('Area: %s    Koordinate: [%d %d %d]    well_id:%d',main_r(i).name,x,y,z,well_id)
            warning off
            TS.area(counter)=main_r(i).name;
            TS.x(counter)=x;
            TS.y(counter)=y;
            TS.z(counter)=z;
            TS.i(counter)=i;
            TS.j(counter)=j;
            TS.well_id(counter)=well_id;
            TS.pval(counter)=map_vol(x,y,z);
            warning on;
            counter=counter+1;
        end
    end
    
end

% indices to unique values in column 3
[~, ind] = unique(TS.well_id);
% duplicate indices
duplicate_ind = setdiff(1:size(TS.well_id, 1), ind);
% duplicate values
duplicate_value = TS.well_id(duplicate_ind);

to_del=[];

for ind_non_unique=1:size(duplicate_ind,2)
    
      [~,pos]=ismember(duplicate_value(ind_non_unique),TS.well_id);
      
      echo_debug=1;

      
    if TS(pos,:).pval    >    TS(duplicate_ind(ind_non_unique),:).pval
        if echo_debug==1
            TS(pos,:)
            disp ('is larger than')
            TS(duplicate_ind(ind_non_unique),:)
            disp(['so delete i=' num2str(TS(duplicate_ind(ind_non_unique),:).i)]);
            disp(['so delete j=' num2str(TS(duplicate_ind(ind_non_unique),:).j)]);
        end
        % get i and j index of main_r
        to_del=[to_del;TS(duplicate_ind(ind_non_unique),:).i,TS(duplicate_ind(ind_non_unique),:).j];
    else
         if echo_debug==1
            TS(pos,:)
            disp ('is smaller than')
            TS(duplicate_ind(ind_non_unique),:)
            disp(['so delete i=' num2str(TS(pos,:).i)]);
            disp(['so delete j=' num2str(TS(pos,:).j)]);
        end
        % get i and j index of main_r
        to_del=[to_del;TS(pos,:).i,TS(pos,:).j];
    end
end


%all rows of main_r which should be deleted are in to_del, this matrix have
%to be flipped to start deleting from the end of the matrix so that ind. do
%not shift
to_del = sortrows(to_del,'descend');

% delete corresponding lines in main_r
for n_2_del=1:size(to_del,1)
        main_r(to_del(n_2_del,1)).validated_zscores(to_del(n_2_del,2),:)=[];  % delete vaidated zscores
        main_r(to_del(n_2_del,1)).data2plot{1, 1}(:,to_del(n_2_del,2))=[];% delete coordinates
        main_r(to_del(n_2_del,1)).data2plot{1, 2}(to_del(n_2_del,2),:)=[];% delete zscores
        main_r(to_del(n_2_del,1)).data2plot{1, 3}(to_del(n_2_del,2),:)=[];% delete raw data
        main_r(to_del(n_2_del,1)).data2plot{1, 4}(to_del(n_2_del,2),:)=[];% delete well_id
end

clear TS

variable_names_types = [["area", "categorical"]; ...
    ["x","double"];...
    ["y","double"];...
    ["z","double"];...
    ["well_id","double"];...
    ["pval","double"]];
% Make table using fieldnames & value types from above
TS = table('Size',[0,size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));

counter=1;

ontology = readtable('tmp_ontology.csv');

for i=1:size(main_r,2)
    %%%%hier muss der Fall abgefangen werden das gegen ABA Label gerechnet
    %%%%wird, dann ist main_r(i).pmap leer!!!
    if ~isempty(main_r(i).pmap)
        map=spm_vol(main_r(i).pmap);
        map_vol=spm_read_vols(map);
    end
    for j=1:size(main_r(i).validated_zscores,1)
        
        x=main_r(i).data2plot{1, 1}(1,j);
        y=main_r(i).data2plot{1, 1}(2,j);
        z=main_r(i).data2plot{1, 1}(3,j);
        if ~isempty(main_r(i).pmap)
        well_id=main_r(i).data2plot{1, 4}(j);
        end
        %sprintf('Area: %s    Koordinate: [%d %d %d]    well_id:%d',main_r(i).name,x,y,z,well_id)
        warning off
        if isempty(main_r(i).pmap)
            %TS.area(counter)='ABA Label';
            TS.area(counter)=ontology.name{ontology.id==str2num(main_r(i).name)};
        else
            TS.area(counter)=main_r(i).name;
        end
        TS.x(counter)=x;
        TS.y(counter)=y;
        TS.z(counter)=z;
        if isempty(main_r(i).pmap)
            TS.well_id(counter)=0;
            TS.pval(counter)=0;
        else
            TS.well_id(counter)=well_id;
            TS.pval(counter)=map_vol(x,y,z);
        end
        TS.pval(counter)=map_vol(x,y,z);
        warning on;
        counter=counter+1;
    end
end