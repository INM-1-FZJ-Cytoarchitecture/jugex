function extract_exp_lvl(mode,parameter,output_mode)

if output_mode==0
    disp_url=0;
    disp_details=0;
elseif output_mode==1
    disp_url=1;
    disp_details=0;   
elseif output_mode==2
    disp_url=0;
    disp_details=1;
elseif output_mode==3
    disp_url=1;
    disp_details=1;
end

%close all;
%clear;
%clc;
tic;
if strcmp(mode,'man')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EDIT VALUES
    % MAPS 
    pmap={'../maps/example1','../maps/example2'};
    %mask={''}; obsolete
    % AREA NAME
    name={'Fp2 R','Fp1 R'};
    % ALLEN BRAIN SPECIMENS TO ANALYZE
    % All 6 AllenBrains (timepoint 12.2014) specimens={ 'H0351.2001', 'H0351.2002', 'H0351.1009', 'H0351.1012', 'H0351.1015', 'H0351.1016'};
    specimens={'H0351.2001', 'H0351.2002', 'H0351.1009', 'H0351.1012', 'H0351.1015', 'H0351.1016'};
    % CANDIDATE GENES
    gene_list='../gene_list/MDD_Gene_List.csv';
    % OUTPUT
    project_name='pro_1';
    output_folder='./output/';
    map_threshold=2;
    search_mode=2; % [1: VOI vs VOI; 2: VOI vs AHBA Label
    struct_id=4009; %ABI Structure Ids
    %%% DO NOT EDIT VALUES BELOW THAT LINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'gui')
    [pmap,name,specimens,gene_list,project_name,output_folder,map_threshold,search_mode,struct_id]=convert_parameter(parameter);
elseif strcmp(mode,'latest')
    load('latest_settings.mat');
    parameter=output_parameter;
    clear output_paramter;
    [pmap,name,specimens,gene_list,project_name,output_folder,map_threshold,search_mode,struct_id]=convert_parameter(parameter);
else
    error('No running mode selected. [man:manual; gui:call via config_gui; latest: using latest settings]');
end

% CORRELATION PLOT (OBSOLET)
plot_correlation=0;
bar_plot_color={'blue','red'};





% maps = struct('pmap',pmap,...
%     'mask',mask,...
%     'name',name,...
%     'bar_plot_color',bar_plot_color,...
%     'validated_zscores','',...
%     'correlation_coeff','',...
%     'specimen','');
% 
% main_r= struct('pmap','',...
%     'mask','',...
%     'name','',...
%     'bar_plot_color','',...
%     'validated_zscores','',...
%     'correlation_coeff','',...
%     'specimen','',...
%     'data2plot','');


maps = struct('pmap',pmap,...
    'name',name,...
    'bar_plot_color',bar_plot_color,...
    'validated_zscores','',...
    'correlation_coeff','',...
    'specimen','');

main_r= struct('pmap','',...
    'name','',...
    'bar_plot_color','',...
    'validated_zscores','',...
    'correlation_coeff','',...
    'specimen','',...
    'data2plot','');



n=25; %egal....

nturns=size(specimens,2)*size(maps,2);
hwait = waitbar(0,'Getting microarray data from ABA Server...');
step=1;

for i=1:size(specimens,2) % loop over all specimens
    waitbar(step / nturns)
    % figure;
    hold on;
    if search_mode==1  % es sollen exp lvl von area gegen area analysiert werden (search mode==1)
        for j=1:size(maps,2) % loop over all maps per specimen
            
            activation_file =  maps(j).pmap;
            %mask_file = maps(j).mask;
            name=maps(j).name;
            c=maps(j).bar_plot_color;
            %n = 25;
            specimen_name = specimens{i};
            
            % clear vars to be sure there is no conflict
            varlist = {'top_corrs','samples','probes','x_label_cell','x_label_str'};
            clear(varlist{:})
            
            
            %% Do the computation.
                %disp_url=0;
                %disp_details=0;
            specimen = download_specimen(specimen_name,disp_url); %Offline Zeile
            [samples,probes,zscores_of_validated_corrs,data_2_plot] = expression_spm_correlation(activation_file,[],n,specimen,gene_list,map_threshold,search_mode,struct_id,disp_url,disp_details,parameter);
            %%% Add zscores to main struct
            maps(j).validated_zscores=zscores_of_validated_corrs;
            maps(j).correlation_coeff=[];
            maps(j).specimen=specimen_name;
            maps(j).data2plot=data_2_plot;
            step=step+1;
        end
    else  % es sollen exp lvl von area gegen ref. struktur analysiert werden (search mode==2)
        % Structure Ids zB in ontology.csv von Allen Datasets
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%  Gleich search mode 1, aber nur für ein Areal
        % check das richtiges areal genommen (anderes feld kann leer
        % sein (sollte leer sein :-) ) GUI sollte bei anwählen einer
        % Radiobox anderes Maps und names Feld löschen und deaktivieren
        activation_file =  maps(1).pmap;    % j durch 1 ersetzt, da nur eine Schleife durchlaufen wird
        %mask_file = maps(1).mask;
        name=maps(1).name;
        c=maps(1).bar_plot_color;
        %n = 25;
        specimen_name = specimens{i};
        
        tmp_search_mode=1;
        
        % clear vars to be sure there is no conflict
        varlist = {'top_corrs','samples','probes','x_label_cell','x_label_str'};
        clear(varlist{:})
        
        
        %% Do the computation.
        specimen = download_specimen(specimen_name,disp_url); 
        [samples,probes,zscores_of_validated_corrs,data_2_plot] = expression_spm_correlation(activation_file,[],n,specimen,gene_list,map_threshold,tmp_search_mode,struct_id,disp_url,disp_details,parameter);
        maps(1).validated_zscores=zscores_of_validated_corrs;
        maps(1).correlation_coeff=[];
        maps(1).specimen=specimen_name;
        maps(1).data2plot=data_2_plot;
        step=step+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %main_r=[main_r,maps];   %%% zwischen result schon mal an main_r übergeben, da beide schritte in else zweig sind
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% exp... noch mal aufrufen aber anderen search query übergeben
        activation_file =  maps(1).pmap;    % j durch 1 ersetzt, da nur eine Schleife durchlaufen wird
       % mask_file = maps(1).mask;
        maps(2).name=num2str(struct_id);
        c=maps(2).bar_plot_color;
        %n = 25;
        specimen_name = specimens{i};
        % clear vars to be sure there is no conflict
        varlist = {'top_corrs','samples','probes','x_label_cell','x_label_str'};
        clear(varlist{:})
        specimen = download_specimen(specimen_name,disp_url); 
        [samples,probes,zscores_of_validated_corrs,data_2_plot] = expression_spm_correlation(activation_file,[],n,specimen,gene_list,map_threshold,search_mode,struct_id,disp_url,disp_details,parameter);
        maps(2).validated_zscores=zscores_of_validated_corrs;
        maps(2).correlation_coeff=[];
        maps(2).specimen=specimen_name;
        maps(2).data2plot=data_2_plot;
        step=step+1;
    end
    main_r=[main_r,maps];
    
    %[h,p_mean_zscores]=ttest2(mean(maps(1).validated_zscores),mean(maps(2).validated_zscores));
    %disp('ttest der correlation coeff (Ortsinfo enthalten, da pmaps Values mit Genexp. korreliert wurde');
end
try
    close(hwait) ;
    %close all;
end
% aufbau main_r

%data_2_plot={coords,zscores,explevels,well_id};

main_r(1) = [];
str=[output_folder filesep 'extracted_data' filesep project_name '_th_' map_threshold '.mat'];
save(str,'main_r')
str=[output_folder filesep 'extracted_data' filesep project_name '_th_' map_threshold '_parameter.mat'];
save(str,'parameter')
if exist(str,'file')==2
    % Construct a questdlg with two options
    choice = questdlg('Extracted data was saved!Would you like to proceed?', ...
        'Extraction of microarray data succesfull', ...
        'Yes','','Later','Yes');
    % Handle response
    switch choice
        case 'Yes'
            disp('Extracted data was saved');
            disp(['path: ' str ]);
            Analysis('gui');
            
        case 'Later'
            disp('Extracted data was saved');
            disp(['path: ' str ]);
            disp('bye');
    end
else
    error(['Could not save downloaded microarray data to .' filesep output_folder filesep 'extracted_data' filesep '. Check write permissions! Please try again']);
end

overview=create_sample_overview(main_r);
str=[output_folder filesep 'extracted_data' filesep 'overview_' project_name '_th_' map_threshold '.mat'];
save(str,'overview');
disp('Overview extracted sample positions');
T = struct2table(overview)
toc;
end

function [pmap,name,specimens,gene_list,project_name,output_folder,map_threshold,search_mode,struct_id]=convert_parameter(parameter)
pmap={parameter.pmap{1},parameter.pmap{2}};
%mask={'colin27T1_forebrain_bothMask.hdr','colin27T1_forebrain_bothMask.hdr'};
name={parameter.name{1},parameter.name{2}};
specimens=parameter.donors;
gene_list=parameter.gene_list{1};
project_name=parameter.project_name{1};
output_folder=parameter.output_folder{1};
map_threshold=parameter.map_threshold{1};
search_mode=parameter.search_mode{1};
struct_id=parameter.struct_id{1};
end