clc;
clear all;
%%% check if SPM available
try
    [v,r] = spm('Ver','',1); r = str2double(r);
    str =sprintf('==> %-40s % s',v,'installed');
    disp(str);
    disp('-----------------------------------------');
catch
    error('SPM cannot be found in MATLAB path. Please install SPM12!');
end


% further dependencies removed

% fex_functions2check={'export_fig'};
% for i=1:size(fex_functions2check,2)
%     try
%         str = which(fex_functions2check{i});
%         if isempty(str)
%             error(['FEX function "' fex_functions2check{i} '" found in MATLAB path.']);
%         else
%             str =sprintf('==> %-40s % s',fex_functions2check{i},'installed');
%             disp(str);
%         end
%     catch
%         error(['FEX function "' fex_functions2check{i} '" found in MATLAB path.']);
%     end
% end
%disp('-----------------------------------------');


%%% Select root folder
        answer = inputdlg('Name of JuGex folder to store all project relevant data (without path)','Name of new working folder?');
        folder_name = uigetdir(pwd,'Select parent folder for the choosen project folder');
        try
            [s,mess,messid] = mkdir(folder_name,answer{1});
            folder_name = fullfile(folder_name,answer{1});
            if strcmp('Directory already exists.',mess)
                str =sprintf('==> %-40s % s',answer{1},'already exist! No new folder was created');
                disp(str);
            else
                str =sprintf('==> %-40s % s',answer{1},'created');
                disp(str);
            end
            disp('-----------------------------------------');
        catch
            error('Could not create folder! Check write permisson.');
        end

%%% change PWD to root folder

cd(folder_name);

disp(['Switching to new created project folder.']);
disp(['Current working directory:' pwd]);
disp(' ');
disp('-----------------------------------------');
disp(['Creating needed child folder in project folder.']);
disp(' ');
%%% create/check needed Folder
folder2create={'gene_list','maps','output',['output' filesep 'img'],['output' filesep 'extracted_data'],['output' filesep 'analyzed_data']};

    for i=1:size(folder2create,2)
        try
            [s,mess,messid] = mkdir(folder2create{i});
            %             if strcmp('Directory already exists.',mess)
            %                 disp([folder2create{i} ' already exist! No new folder was created']);
            %             else
            %                 disp([folder2create{i} ' created']);
            %             end
        catch
            err_str=['Could not create folder' folder2create{i} ' Check write permisson'];
            error(err_str);
        end
    end
    if exist('err_str','var')==0
        %disp(['==> ' v ' installed']);
        for i=1:size(folder2create,2)
            %disp(['==> ' folder2create{i} ' created']);
            str =sprintf('==> %-40s % s',folder2create{i},'created');
            disp(str);
        end
        disp('-----------------------------------------');
        %disp(['==> ' pwd ' is current working directory']);
    end
    
% elseif strcmp('Select',choice)
%     % check if correct folder exist
%     for i=1:size(folder2create,2)
%         if isdir(folder2create{i}) ~= 1
%             err_str='Folder missing! Restart programm and select "Create new folder"';
%             error('ReGen:FolderChk', err_str)
%         end
%     end
%     if exist('err_str','var')==0
%         %disp(['==> ' v ' installed']);
%         for i=1:size(folder2create,2)
%             %disp(['==> ' folder2create{i} ' exist']);
%             str =sprintf('==> %-40s % s',folder2create{i},'exist');
%             disp(str);
%         end
%         disp('-----------------------------------------');
%         disp(['==> ' pwd ' is current working directory']);
%     end
% end
disp('Calling Configuration script to start Jugex tools.');
disp('-----------------------------------------');
fprintf('If you have allready called the JuGex script and would\nlike to continiue working at a project please navigate\nto the specific project folder and type\n\n');
disp('Configuration		to configure new analysis and extract ABA gene expression data');
disp('Analysis		    to analyze the extracted data');
disp('Visualization		to visualize the results'); 	
disp('-----------------------------------------');
Configuration;