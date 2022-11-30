function main_r = apply_post_sf(main_r,delete_status,checked_unique_ABA_label)
%apply_post_sf filter structure names from main_r var and select only
%needed ones
%   Detailed explanation goes here

if delete_status==1 && istable(checked_unique_ABA_label)
    disp('do remove in if block');
    fig_sf=findall(groot,'Tag','fig_sf');
    disp('close and continue');
    
    label_2del=checked_unique_ABA_label.structure_name(checked_unique_ABA_label.use==false,:);
    
    for i=1:size(main_r,2)%loop over all rows of main_r struct (eq. maps)
        if size(main_r(i).validated_zscores,1)>0 %only do something if there were TBs detected
            del_vector=zeros(1,size(main_r(i).validated_zscores,1));
            for j=1:size(main_r(i).validated_zscores,1)
                if any(strcmp(label_2del,main_r(i).data2plot{1,5}{j,1}))
                % Do Something
                    disp([main_r(i).data2plot{1,5}{j,1} ' should be deleted  (row: ' num2str(i) ' ; 2nd row: ' num2str(j)]);
                    del_vector(j)=1;
                else
                % Do Something else                    
                    disp([main_r(i).data2plot{1,5}{j,1} ' should NOT be deleted  (row: ' num2str(i) ' ; 2nd row: ' num2str(j)]);
                end
                
            end
            main_r(i).validated_zscores(find(del_vector),:)=[];
        end    
        
    end   
    %close(fig_sf);
    disp('deleted');
else
    
    % extract all structure names from main_r struct
    extracted_struct_names=table();
    for i=1:size(main_r,2)
        extracted_struct_names=[extracted_struct_names;cell2table(main_r(i).data2plot{1, 5})];
    end
    % create unique table
    unique_extracted_struct_names=unique(extracted_struct_names);
    n_unique_structure_names=size(unique_extracted_struct_names,1);
    
    unique_extracted_struct_names.Properties.VariableNames{1} = 'structure_name';
    unique_extracted_struct_names(1:n_unique_structure_names,2)=array2table(logical(ones(n_unique_structure_names,1)));
    unique_extracted_struct_names.Properties.VariableNames{2} = 'use';
    
    fig_sf = uifigure('Position',[100 100 600 500]);
    fig_sf.Tag='fig_sf';
    uit = uitable(fig_sf,'Data',unique_extracted_struct_names,'Position',[20 40 560 440]);
    uit.ColumnEditable=[false,true];
    
    
    % Create a push button
    btn = uibutton(fig_sf,'push',...
        'Position',[420, 218, 100, 22],...
        'ButtonPushedFcn', @(btn,event) plotButtonPushed(btn,uit.Data,main_r,fig_sf));
    uiwait(fig_sf);
end
end

% Create the function for the ButtonPushedFcn callback
function plotButtonPushed(btn,data,main_r,fig_sf)
%uiresume(fig_sf);
disp('#################################');
disp('Should tissue blocks with the following ABA Labels be ignored?'); 
disp('   (Tisse Blocks with those label will be removed and the remainig tissue blocks');
disp('    will be saved and used for subsequent analysis)');
disp('#################################');
data(data.use==false,:)
prompt = "Remove tissue blocks with those label ? Y/N [Y]: ";
txt = input(prompt,"s");
if isempty(txt)
    txt = 'Y';
end
if txt=='Y'
    disp('do remove');
    %btn.UserData=data;
    close(fig_sf)
    main_r=apply_post_sf(main_r,1,data);
    %main_r='neu';
    %return
elseif txt=='N'
    disp('No AllenBrain Label were removed!');
    close(fig_sf)
    %return
else
    disp('Unknown selection! Please try Y/N [Y]:');
    close(fig_sf)
    main_r=apply_post_sf(main_r,0,data);
    %return
end
end