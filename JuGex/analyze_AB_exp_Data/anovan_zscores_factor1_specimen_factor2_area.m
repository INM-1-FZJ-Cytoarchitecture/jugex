function anovan_zscores_factor1_specimen_factor2_area(mode,parameter)

tic;
p2use=0.05;
%clearvars;
if strcmp(mode,'man')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EDIT VALUES
    path_result_mat='../Output_extract_AB_exp_Data/main_r_Fp1_Fp2_L_plus_sero_02th_5donors_nogauss.mat';
    gen_list='../gene_list/MDD_Gene_List.csv';
    n_rep=100;
    area1_name='Fp1 L'; % identic+case sensitive to main_r.name
    area2_name='Fp2 L'; % identic+case sensitive to main_r.name
    search_mode=1;  %[1= gensymbol mode; 2= probe_id mode]
    %%% DO NOT EDIT CODE BELOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'gui')
    [path_result_mat,gen_list,n_rep,area1_name,area2_name,search_mode]=convert_parameter(parameter);
else
    error('No running mode selected. [man:manual; gui:call via config_gui]');
end
% read gene list
[probe_id,gene_symbol,entrez_id] = read_gen_file(gen_list);

load(path_result_mat);

p_value_based_sample_allocation=1;

if ~isempty(area1_name) && ~isempty(area2_name)  % If one map name is empty we compare map vs ABA label, so we can not assign sample by pval since there is no pval for ABA label ;-)
    
    if p_value_based_sample_allocation==1
        % function pval_filter is saved in utilities and removes duplicate
        % entries of main_r only the TS with the higher pval will be kept
        [main_r,TS]=pval_filter(main_r);
    end
end

[~,name_res_mat,~]=fileparts(path_result_mat);
th = name_res_mat(end-3:end);
[~,gene_liste_name,~]=fileparts(gen_list);
if search_mode==1
    sm_str='all-probes_mode';
elseif search_mode==2
    sm_str='single-probe_mode';
end
diary_name=sprintf('%s_vs_%s_genelist_%s_nrep_%s_%s_%s.txt',area1_name,area2_name,gene_liste_name,num2str(n_rep),th,sm_str);
diary(diary_name);


area1_zscores=[];
area2_zscores=[];
area1_specimen=[];
area2_specimen=[];
area1_area=[];
area2_area=[];
area1_koords=[];
area2_koords=[];

% create struct & groups for anovan
for i=1:size(main_r,2)
    if strcmp(main_r(i).name,area1_name)
        area1_zscores=[area1_zscores;main_r(i).validated_zscores];
        for j=1:size(main_r(i).validated_zscores,1)
            area1_specimen{end+1}=main_r(i).specimen;
            area1_area{end+1}=main_r(i).name;
            area1_koords(end+1,:)=[main_r(i).data2plot{1, 1}(1,j),main_r(i).data2plot{1, 1}(2,j),main_r(i).data2plot{1, 1}(3,j)];
        end
    else
        area2_zscores=[area2_zscores;main_r(i).validated_zscores];
        for j=1:size(main_r(i).validated_zscores,1)
            area2_specimen{end+1}=main_r(i).specimen;
            area2_area{end+1}=main_r(i).name;
            area2_koords(end+1,:)=[main_r(i).data2plot{1, 1}(1,j),main_r(i).data2plot{1, 1}(2,j),main_r(i).data2plot{1, 1}(3,j)];
        end
    end
end
%clc;

factor_specimen=[area1_specimen area2_specimen];
factor_area=[area1_area area2_area];
combined_zscores=[area1_zscores;area2_zscores];

n_samples=size(combined_zscores,1);
n_samples_area1=size(area1_area,2);
n_samples_area2=size(area2_area,2);
n_probeids=size(combined_zscores,2);
%%%% create additional factors (age, gender, race) based on specimen_info
specimen_info=build_specimen_info();
for counter=1:n_samples
    info_index=find(strcmp({specimen_info.name},factor_specimen(counter)));
    diff_specimen=unique(factor_specimen);
    for r=1:size(diff_specimen,2)
        if strcmp(factor_specimen{counter},diff_specimen{r})
            specimen_n=r;
        end
    end
    factor_specimen_numeric(counter)=specimen_n;
    
    diff_area=unique(factor_area);
    for r=1:size(diff_area,2)
        if strcmp(factor_area{counter},diff_area{r})
            area_n=r;
        end
    end
    factor_area_numeric(counter)=area_n;
    
    factor_age{counter}=num2str(specimen_info(info_index).age);
    factor_age_numeric(counter)=specimen_info(info_index).age;
    switch specimen_info(info_index).age  % continous für age ging nicht???
        case num2cell(1:10)
            fg=1;
        case num2cell(11:20)
            fg=2;
        case num2cell(21:30)
            fg=3;
        case num2cell(31:40)
            fg=4;
        case num2cell(41:50)
            fg=5;
        case num2cell(51:60)
            fg=6;
        case num2cell(61:70)
            fg=7;
        case num2cell(71:80)
            fg=8;
        case num2cell(81:90)
            fg=9;
    end
    factor_age_grouped{counter}=num2str(fg);
    factor_age_grouped_numeric(counter)=fg;
    factor_gender{counter}=specimen_info(info_index).gender;
    if strcmp(specimen_info(info_index).gender,'M')
        gb=0;
    elseif strcmp(specimen_info(info_index).gender,'F')
        gb=1;
    end
    factor_gender_binary(counter)=gb;
    factor_race{counter}=specimen_info(info_index).race;
    diff_race=unique({specimen_info.race});
    for r=1:size(diff_race,2)
        if strcmp(specimen_info(info_index).race,diff_race{r})
            race_n=r;
        end
    end
    factor_race_numeric(counter)=race_n;
end
%%%% switch Search mode; 1== combine genes using winzorized_mean; 2== using
%%%% individual probe_id values
if search_mode==1
    if p_value_based_sample_allocation~=1
    % function pval_filter is saved in utilities and removes duplicate
    % entries of main_r only the TS with the higher pval will be kept
    % without it the debug info will not be visible
    clc;
    end
    
    disp(' ');
    disp('########################################');
    disp('####      Gene-level Analysis      #####');
    disp('####                               #####');
    disp('####        All-probes mode        #####');
    disp('########################################');
    [ combined_zscores,area1_zscores,area2_zscores,unique_entrez_id  ] = switch2gensymbol_lvl( entrez_id,combined_zscores,n_samples_area1,n_samples_area2 );
    n_genes=size(combined_zscores,2);
    %%%% Referenz Anovan
    Reference_Anovan_p=zeros(1,n_genes);
    Reference_Anovan_eta2=zeros(1,n_genes);
    Reference_Anovan_CI_l=zeros(1,n_genes);
    Reference_Anovan_CI_h=zeros(1,n_genes);
    Reference_Anovan_diff_mean=zeros(1,n_genes);
    F_vec_ref_anovan=zeros(1,n_genes);
elseif search_mode==2
    if p_value_based_sample_allocation~=1
    % function pval_filter is saved in utilities and removes duplicate
    % entries of main_r only the TS with the higher pval will be kept
    % without it the debug info will not be visible
    clc;
    end
    disp(' ');
    disp('########################################');
    disp('####      Gene-level Analysis      #####');
    disp('####                               #####');
    disp('####       Single-probe mode       #####');
    disp('########################################');
    %%%% Referenz Anovan
    Reference_Anovan_p=zeros(1,n_probeids);
    Reference_Anovan_eta2=zeros(1,n_probeids);
    Reference_Anovan_CI_l=zeros(1,n_probeids);
    Reference_Anovan_CI_h=zeros(1,n_probeids);
    Reference_Anovan_diff_mean=zeros(1,n_probeids);
    F_vec_ref_anovan=zeros(1,n_probeids);
end

if search_mode==2
    n_genes=size(combined_zscores,2);
end

% filter TS which are included in both maps, the TS will be used only for
% the area with the higher p val at the TS position



%%%% define factors to analyze
% factors={factor_area_numeric factor_specimen_numeric}; varnames = {'Area';'Specimen'};
factors={factor_area factor_specimen factor_age_numeric factor_race}; varnames = {'Area';'Specimen';'Age';'Race'};
%factors={factor_area_numeric factor_specimen_numeric factor_gender_binary factor_age_numeric}; varnames = {'Area';'Specimen';'Gender';'Age'};
% factors={factor_area_numeric factor_specimen_numeric factor_gender_binary}; varnames = {'Area';'Specimen';'Gender'};
% factors={factor_area factor_specimen factor_race}; varnames = {'Area';'Specimen';'Race'};
% factors={factor_area factor_specimen factor_age factor_gender}; varnames = {'Area';'Specimen';'Age';'Gender'};
% factors={factor_area factor_specimen factor_age factor_race}; varnames = {'Area';'Specimen','Age','Race'};
% factors={factor_area factor_specimen factor_age factor_race factor_gender}; varnames = {'Area';'Specimen';'Age';'Race';'Gender'};

for i=1:n_genes
    [p,tab,stats] = anovan(combined_zscores(:,i),factors,'varnames',varnames,'continuous',3,'display','off');
    %[p,tab,stats] = anovan(combined_zscores(:,i),factors,'varnames',varnames,'display','off','continuous',4);
    %[p,tab,stats] = anovan(combined_zscores(:,i),factors,'varnames',varnames,'display','off');
    F_vec_ref_anovan(i) = tab{2,6};
    
    ss_total=tab{7,2};
    ss_between_group_area=tab{2,2};
    Reference_Anovan_eta2(i)=ss_between_group_area/ss_total;
    
    % A confidence interval on the difference between means is computed using the following formula:
    % Lower Limit = M1 - M2 -(tCL)(S m1-m2)
    % Upper Limit = M1 - M2 +(tCL)(S m1-m2)
    %
    % where M1 - M2 is the difference between sample means,
    % tCL is the t for the desired level of confidence, and S m1-m2
    % is the estimated standard error of the difference between sample means.
    
    mse=(var(combined_zscores(factor_area_numeric==1,i))+ var(combined_zscores(factor_area_numeric==2,i)))/2;
    sm1m2=sqrt((2*mse)/n_genes);
    % t for 95% with 48 degree of freedom confidence interval	2.011
    % (http://onlinestatbook.com/2/calculators/inverse_t_dist.html)
    Reference_Anovan_CI_l(i)=(mean(combined_zscores(factor_area_numeric==1,i))- mean(combined_zscores(factor_area_numeric==2,i)))-2.011*sm1m2;
    Reference_Anovan_CI_h(i)=(mean(combined_zscores(factor_area_numeric==1,i))- mean(combined_zscores(factor_area_numeric==2,i)))+2.011*sm1m2;
    
    
    
    Reference_Anovan_diff_mean(i)=mean(combined_zscores(factor_area_numeric==1,i))-mean(combined_zscores(factor_area_numeric==2,i));
    Reference_Anovan_p(i)=p(1);
end
%%%% FDR correction using function based on orig Benjamini paper (fex)
% [h_anovan, crit_p_anovan, adj_p_anovan]=fdr_bh(Reference_Anovan_p,.05,'pdep','yes');

specimens=unique(factor_specimen);
F_mat_perm_anovan=zeros(n_rep,n_genes);
p_mat_perm_anovan=zeros(n_rep,n_genes);
F_mat_perm_anovan(1,:)=F_vec_ref_anovan;
h_waitbar = waitbar(0,'busy...');
%%%% n_rep Anovan repetitions with pairwised permutated data
for rep=2:n_rep
    % Anovan mit perm area1/area2
    shuffle = randperm(n_samples); % es werden einmal(pro schleife) die label (area name) geshuffelt und dann für die Anovan benutzt
    F_vec_perm_anovan=zeros(1,n_genes);
    p_vec_perm_anovan=zeros(1,n_genes);
    
    % Check that user has the distrib_toolbox Toolbox installed.
    hasDCT = license('test', 'distrib_computing_toolbox');
    if ~hasDCT
        % User does not have the toolbox installed.
        %disp('No Parallel Computing Toolbox installed. Will perform in normal mode!')
        for i=1:n_genes
            % an der Stelle kann age, gender, ethiks als factor eingebracht
            % werden. Nur area wird geschuffelt!
            % Es können auch mehr als 3 areas eingelesen werden, da nur die lables geschuffelt werden
            
            %[p,tab,stats] = anovan(combined_zscores(:,i),{factor_area_numeric(shuffle) factor_specimen_numeric factor_gender_binary factor_age_numeric},'varnames',varnames,'display','off','continuous',4);
            [p,tab,stats] = anovan(combined_zscores(:,i),{factor_area(shuffle) factor_specimen factor_age_numeric factor_race},'varnames',varnames,'continuous',3,'display','off');
            F_vec_perm_anovan(i) = tab{2,6}; % n(Genes or probe_ids) F Werte als Vector
            % test nur fürs Histogramm
            p_vec_perm_anovan(i) = p(1);
            % test ende
            % n(Genes or probe_ids) F Werte mal n rep in matrix schreiben
        end
    else
        %disp('Parallel Computing Toolbox installed. Starting parallel pool!')
        parfor i=1:n_genes
            % an der Stelle kann age, gender, ethiks als factor eingebracht
            % werden. Nur area wird geschuffelt!
            % Es können auch mehr als 3 areas eingelesen werden, da nur die lables geschuffelt werden
            
            %[p,tab,stats] = anovan(combined_zscores(:,i),{factor_area_numeric(shuffle) factor_specimen_numeric factor_gender_binary factor_age_numeric},'varnames',varnames,'display','off','continuous',4);
            [p,tab,stats] = anovan(combined_zscores(:,i),{factor_area(shuffle) factor_specimen factor_age_numeric factor_race},'varnames',varnames,'continuous',3,'display','off');
            F_vec_perm_anovan(i) = tab{2,6}; % n(Genes or probe_ids) F Werte als Vector
            % test nur fürs Histogramm
            p_vec_perm_anovan(i) = p(1);
            % test ende
            % n(Genes or probe_ids) F Werte mal n rep in matrix schreiben
        end
    end
    %     parfor i=1:n_genes
    %         % an der Stelle kann age, gender, ethiks als factor eingebracht
    %         % werden. Nur area wird geschuffelt!
    %         % Es können auch mehr als 3 areas eingelesen werden, da nur die lables geschuffelt werden
    %
    %         %[p,tab,stats] = anovan(combined_zscores(:,i),{factor_area_numeric(shuffle) factor_specimen_numeric factor_gender_binary factor_age_numeric},'varnames',varnames,'display','off','continuous',4);
    %         [p,tab,stats] = anovan(combined_zscores(:,i),{factor_area(shuffle) factor_specimen factor_age_numeric factor_race},'varnames',varnames,'continuous',3,'display','off');
    %         F_vec_perm_anovan(i) = tab{2,6}; % n(Genes or probe_ids) F Werte als Vector
    %         % test nur fürs Histogramm
    %         p_vec_perm_anovan(i) = p(1);
    %         % test ende
    %         % n(Genes or probe_ids) F Werte mal n rep in matrix schreiben
    %     end
    F_mat_perm_anovan(rep,:)=F_vec_perm_anovan;
    p_mat_perm_anovan(rep,:)=p_vec_perm_anovan;
    waitbar(rep / n_rep)
end
close(h_waitbar);
disp(' ');

%%%% F_vec_ref_anovan   referenz F Werte der Anova
%%%% F_mat_perm_anovan  Matrix der permutierten F Werte dim(rep x n probes)

ref = max(F_mat_perm_anovan,[],2);
Uncorrected_permuted_p=zeros(1,n_genes);

for i=1:n_genes
    Uncorrected_permuted_p(i)=sum(F_mat_perm_anovan(:,i)>=F_vec_ref_anovan(:,i))/n_rep;
    FWE_corrected_p(i)=mean(ref>=F_vec_ref_anovan(i));  % mean als kurzform für summe(werte)/anzahl (siehe zeile davor)
end

% FDR auf Uncorrected_permuted_p (besser als obere FDR, da wir oben p Werte aus nicht
% permutierten Anovans korrigiert haben. Wegen der geringen sample size
% sollten die Werte aber auf der annahme der lable austauschbarkeit
% permutiert werden.
%  [h_anovan, crit_p_anovan, FDR_corrected_p]=fdr_bh(Uncorrected_permuted_p,p2use,'pdep','yes');
%
% figure(99), subplot(2,2,1); scatter(Reference_Anovan_p,Uncorrected_permuted_p), title('Classical vs. Permutation')
% figure(99), subplot(2,2,2); scatter(Uncorrected_permuted_p,FWE_corrected_p), title('Uncorrected vs. FWE')
% figure(99), subplot(2,2,3); scatter(Uncorrected_permuted_p,FDR_corrected_p), title('Uncorrected vs. FDR')
% figure(99), subplot(2,2,4); scatter(FWE_corrected_p,FDR_corrected_p), title('FWE vs. FDR')

if search_mode==1
    probe_id=[];
    
    entrez_id=unique(entrez_id);
    gene_info = table2struct(readtable(gen_list));
    clear gene_symbol;
    for i=1:size(entrez_id,1)
        probe_id{i}='NaN';
        ind_search_entrezid=find([gene_info.entrez_id] == entrez_id(i));
        gene_symbol{i}=gene_info(ind_search_entrezid(1)).gene_symbol;
    end
end
[path,name,ext]=fileparts(path_result_mat);

str_title=['Results n-way ANOVA with user-specified permutations (n=' num2str(n_rep) ')'];

%%%%% find path from current m file location
parts = strsplit(mfilename('fullpath'), filesep);
jugex_root_path = strjoin(parts(1,1:end-2),filesep);

reference_brain=[jugex_root_path '\visualize_AB_exp_Data\visualization_template\mni_icbm152_t1_masked.nii'];


spm_header_reference_brain=spm_vol(reference_brain);
spm_header_reference_origin=spm_header_reference_brain.mat(1:3,4);

spm_header_map=spm_vol(parameter.pmap{1,1});
spm_header_map_origin=spm_header_map.mat(1:3,4);

% shift because jubrain maps with 256x256x256 dim and mni152 193x229x193
% dim
shift_between_map_and_ref=abs(spm_header_reference_origin)-abs(spm_header_map_origin);


if search_mode==1
    table_header={'Index' 'Entrez_ID' 'Gene_Symbol' 'p_FWE_corrected'};
    
    % edit rev2 print results_file
    res_file=fopen(['data_' name '.txt'],'w');
    %create header line of txt file
    fprintf(res_file,'VOI\tMNI\t');
    fprintf(res_file,'%s\t',gene_symbol{:});
    fprintf(res_file,'\n');
    voi_koords=[area1_koords;area2_koords];
   
    
    origin_correct=repmat(spm_header_reference_origin',size(voi_koords,1),1);
    shift_correct=repmat(shift_between_map_and_ref',size(voi_koords,1),1);
   
    voi_koords=(voi_koords+origin_correct)+shift_correct;
    
    for idx_res=1:size(combined_zscores,1)
        fprintf(res_file,'%s\t',factor_area{idx_res});
        %fprintf(res_file,'%f\t',Uncorrected_permuted_p(idx_res));
        %fprintf(res_file,'%f\t',FWE_corrected_p(idx_res));
        fprintf(res_file,'%d,%d,%d\t',voi_koords(idx_res,1),voi_koords(idx_res,2),voi_koords(idx_res,3));
        fprintf(res_file,'%f\t',combined_zscores(idx_res,:));
        fprintf(res_file,'\n');
    end
    fprintf(res_file,'\n');
    fprintf(res_file,'\n');
    fprintf(res_file,'p_uc\t\t');
    fprintf(res_file,'%f\t',Uncorrected_permuted_p);
    fprintf(res_file,'\n');
    fprintf(res_file,'p_FWE\t\t');
    fprintf(res_file,'%f\t',FWE_corrected_p);
    fprintf(res_file,'\n');
    fclose(res_file);
    % ende edit rev2
    
    
elseif search_mode==2
    table_header={'Index' 'Probe_name' 'Entrez_ID' 'Gene_Symbol' 'p_FWE_corrected'};
    
    % edit rev2 print results_file
    res_file=fopen(['data_' name '.txt'],'w');
    %create header line of txt file
    fprintf(res_file,'VOI\tMNI\t');
    fprintf(res_file,'%d\t',probe_id);
    fprintf(res_file,'\n');
    voi_koords=[area1_koords;area2_koords];
    
        origin_correct=repmat(spm_header_reference_origin',size(voi_koords,1),1);
    shift_correct=repmat(shift_between_map_and_ref',size(voi_koords,1),1);
   
    voi_koords=(voi_koords+origin_correct)+shift_correct;
    for idx_res=1:size(combined_zscores,1)
        fprintf(res_file,'%s\t',factor_area{idx_res});
        %fprintf(res_file,'%f\t',Uncorrected_permuted_p(idx_res));
        %fprintf(res_file,'%f\t',FWE_corrected_p(idx_res));
        fprintf(res_file,'%d,%d,%d\t',voi_koords(idx_res,1),voi_koords(idx_res,2),voi_koords(idx_res,3));
        fprintf(res_file,'%f\t',combined_zscores(idx_res,:));
        fprintf(res_file,'\n');
    end
    fprintf(res_file,'\n');
    fprintf(res_file,'\n');
    fprintf(res_file,'p_uc\t\t');
    fprintf(res_file,'%f\t',Uncorrected_permuted_p);
    fprintf(res_file,'\n');
    fprintf(res_file,'p_FWE\t\t');
    fprintf(res_file,'%f\t',FWE_corrected_p);
    fprintf(res_file,'\n');
    fclose(res_file);
    % ende edit rev2
    
    
end


index_sig_gens=find(FWE_corrected_p<p2use);
cw_output(index_sig_gens,probe_id,entrez_id,gene_symbol,FWE_corrected_p,str_title,table_header,'FWE','005',name,search_mode,n_rep);

%%% plots
do_further_plots=0;
if do_further_plots~=0
    for i=1:size(index_sig_gens,2)
        id=index_sig_gens(i);
        fign=gcf;
        fign=fign.Number+1;
        if isempty(fign)
            fign=1;
        end
        factor_age_numeric=cellfun(@str2num, factor_age(:,1:end));
        figure(fign);
        min_z_scores=min(combined_zscores(:,id));
        max_z_scores=max(combined_zscores(:,id));
        set(gcf,'numbertitle','off','Name',gene_symbol{id});
        figure(fign), subplot(6,2,1);scatter(factor_age_numeric(1:n_samples_area1),area1_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. Age Area1'),lsline(gca);
        figure(fign), subplot(6,2,2);scatter(factor_age_numeric(n_samples_area1+1:n_samples_area1+n_samples_area2),area2_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. Age Area2'),lsline(gca);
        for j=1:size(unique(factor_gender),2)
            u_fg=unique(factor_gender);
            index = find(strcmp(factor_gender, u_fg{j}));
            n_fg(index) = j;
        end
        figure(fign), subplot(6,2,3);scatter(n_fg(1:n_samples_area1),area1_zscores(:,id)),ylim([min_z_scores max_z_scores]),ylim([min_z_scores max_z_scores]), title('z-scores vs. Gender Area1'),lsline(gca);
        figure(fign), subplot(6,2,4);scatter(n_fg(n_samples_area1+1:n_samples_area1+n_samples_area2),area2_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. Gender Area2'),lsline(gca);
        for j=1:size(unique(factor_race),2)
            u_fr=unique(factor_race);
            index = find(strcmp(factor_race, u_fr{j}));
            n_fr(index) = j;
        end
        figure(fign), subplot(6,2,5);scatter(n_fr(1:n_samples_area1),area1_zscores(:,id)),ylim([min_z_scores max_z_scores]),ylim([min_z_scores max_z_scores]), title('z-scores vs. ethnic group Area1'),lsline(gca);
        figure(fign), subplot(6,2,6);scatter(n_fr(n_samples_area1+1:n_samples_area1+n_samples_area2),area2_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. ethnic group Area2'),lsline(gca);
        figure(fign), subplot(6,2,7);scatter(area1_koords(:,2),area1_zscores(:,id)),ylim([min_z_scores max_z_scores]),ylim([min_z_scores max_z_scores]), title('z-scores vs. X-position Area1'),lsline(gca);
        figure(fign), subplot(6,2,8);scatter(area2_koords(:,2),area2_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. X-position Area2'),lsline(gca);
        figure(fign), subplot(6,2,9);scatter(area1_koords(:,1),area1_zscores(:,id)),ylim([min_z_scores max_z_scores]),ylim([min_z_scores max_z_scores]), title('z-scores vs. Y-position Area1'),lsline(gca);
        figure(fign), subplot(6,2,10);scatter(area2_koords(:,1),area2_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. Y-position Area2'),lsline(gca);
        figure(fign), subplot(6,2,11);scatter(area1_koords(:,3),area1_zscores(:,id)),ylim([min_z_scores max_z_scores]),ylim([min_z_scores max_z_scores]), title('z-scores vs. Z-position Area1'),lsline(gca);
        figure(fign), subplot(6,2,12);scatter(area2_koords(:,3),area2_zscores(:,id)),ylim([min_z_scores max_z_scores]), title('z-scores vs. Z-position Area2'),lsline(gca);
        
    end
end
%%%
fprintf('\nNum samples %s: %d \n',area1_name,n_samples_area1);
fprintf('Num samples %s: %d \n',area2_name,n_samples_area2);
fprintf('Num samples      : %d \n',n_samples);
if search_mode==1
    fprintf('Num genes        : %d \n',n_genes);
elseif search_mode==2
    fprintf('Num probes        : %d \n',n_genes);
end
fprintf('Num permutations : %d \n\n',n_rep);
fprintf('Gene list        : %s \n',gene_liste_name);
toc;
diary off
end




function [path_result_mat,gen_list,n_rep,area1_name,area2_name,search_mode]=convert_parameter(parameter)

path_result_mat=parameter.extracted_data_file{1};
gen_list=parameter.gene_list{1};
n_rep=parameter.nrepetitions;
area1_name=parameter.name{1}; % identic+case sensitive to main_r.name
area2_name=parameter.name{2}; % identic+case sensitive to main_r.name
search_mode=parameter.analysis_mode;  %[1= gensymbol mode; 2= probe_id mode]

end
