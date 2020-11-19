function do_boxplot(gen_list,mode,parameter,p_r,varargin)
% Some function that requires 2 inputs and has some optional inputs.

% only want 3 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:sdo_boxplot:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

%n_probes_or_genes=285; %für probe_ids
n_probes_or_genes=25; %für gene_level

% set defaults for optional inputs
optargs = {1 n_probes_or_genes};

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[start_g_p, end_g_p] = optargs{:};







if strcmp(mode,'gui')
    % read gene list
    [probe_id,gene_symbol,entrez_id] = read_gen_file(gen_list);
    % read main_r structure [output of spm correlation demo]
    [path_result_mat,gen_list,n_rep,area1_name,area2_name,search_mode]=convert_parameter(parameter);
    load(path_result_mat);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EXAMPLE VALUES
    % path_result_mat='../Output_extract_AB_exp_Data/main_r_Fp1_Fp2_L_plus_sero_02th_5donors_nogauss.mat';
    % gen_list='../gene_list/MDD_Gene_List.csv';
    % n_rep=100;
    % area1_name='Fp1 L'; % identic+case sensitive to main_r.name
    % area2_name='Fp2 L'; % identic+case sensitive to main_r.name
    % search_mode=1;  %[1= gensymbol mode; 2= probe_id mode]
    %% DO NOT EDIT CODE BELOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    area1_zscores=[];
    area2_zscores=[];
    n_probes=size(area1_zscores,2);
    
    for i=1:size(main_r,2)
        if strcmp(main_r(i).name,main_r(1).name)
            area1_zscores=[area1_zscores;main_r(i).validated_zscores];
        else
            area2_zscores=[area2_zscores;main_r(i).validated_zscores];
        end
    end
    
    combined_zscores=[area1_zscores;area2_zscores];
    labels=probe_id(start_g_p:end_g_p);
    
    if parameter.analysis_mode==1
        [ combined_zscores,area1_zscores,area2_zscores,unique_entrez_id  ] = switch2gensymbol_lvl( entrez_id,combined_zscores,size(area1_zscores,1),size(area2_zscores,1) );
        
        labels=unique_entrez_id(start_g_p:end_g_p);
    end
    
    
    
    
    
    
    boxplot_fig_h=figure;
    
    fig = gcf;
    hold on
    
    position_fp1 = 1:1:n_probes;
    h_sub3=subplot (1,1,1);
    
    box_O = boxplot([area1_zscores(:,start_g_p:end_g_p)],'colors','r','positions',position_fp1,'width',0.18,'orientation','horizontal'); %,'orientation','horizontal'
    set(gca,'XTickLabel',{' '})  % Erase xlabels
    hold on  % Keep
    position_S = 1.3:1:n_probes+.3;  % Define position for 12 Month_S boxplots
    box_S = boxplot([area2_zscores(:,start_g_p:end_g_p)],'colors','b','positions',position_S,'width',0.18,'orientation','horizontal');
    xlim([-2 2])
    
    title_str=['zscores explvl  ' area1_name '(red) ' area2_name ' (blue)'];
    title_str= strrep(title_str, '_', '\_');
    title(title_str);
    set(h_sub3,'ytick',1:n_probes_or_genes, 'yticklabel',num2str(labels))
    lims=ylim;
    for ttt=1:size(labels,1)
        labels_gene_symbols(ttt)=gene_symbol(min(find(entrez_id(:,1)==unique_entrez_id(ttt))))
    end
    labels_gene_symbols=labels_gene_symbols';
    %     for zzz=1:size(labels,1)
    %     status=0;
    % while status<1
    %     try
    % url=['http://api.brain-map.org/api/v2/data/query.json?criteria=model::Probe,rma::criteria,[id$eq''' num2str(labels(zzz)) ''']'];
    % [str,status] = urlread(url);
    %     catch
    %         disp('lost connection. I start a new attempt!');
    %     end
    % end
    %
    % json = parse_json(str);
    % cust_label{zzz}=json.msg{1,1}.label
    % json.msg{1,1}.id;
    % json.msg{1,1}.label;
    
    yyaxis right
    set(h_sub3,'ytick',1:31, 'yticklabel',cust_label)
    
%end

fig = gcf;
set(boxplot_fig_h, 'PaperType', 'A4');
set(boxplot_fig_h, 'PaperOrientation', 'portrait');
set(boxplot_fig_h, 'PaperUnits', 'centimeters');
set(boxplot_fig_h, 'PaperPositionMode', 'manual');
set(boxplot_fig_h,'PaperPosition',[0.63 0.63 19.72 28.41]);
set(boxplot_fig_h,   'PaperSize',[20.98 29.68]);
%print -dtiffnocompression -r300 boxplot_p
print(boxplot_fig_h,'boxplot_rrr','-dtiffnocompression','-r300')
%export_fig 'boxplot.tif' -r300

else
    disp('boxplot can not be plotted in manual mode. please call the function via the GUI.')
end
end % ende do_boxplot function

function [path_result_mat,gen_list,n_rep,area1_name,area2_name,search_mode]=convert_parameter(parameter)

path_result_mat=parameter.extracted_data_file{1};
gen_list=parameter.gene_list{1};
n_rep=parameter.nrepetitions;
area1_name=parameter.name{1}; % identic+case sensitive to main_r.name
area2_name=parameter.name{2}; % identic+case sensitive to main_r.name
search_mode=parameter.analysis_mode;  %[1= gensymbol mode; 2= probe_id mode]

end