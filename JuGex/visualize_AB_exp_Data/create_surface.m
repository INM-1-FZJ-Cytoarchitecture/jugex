function create_surface(mode,parameter)
%clear all;
%close all;
tic;
if strcmp(mode,'man')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EDIT VALUES
    reference_brain=load_nii('mni_icbm152_t1_masked.nii');
    path_result_mat='.\output\extracted_data\fp1_L_vs_CC_th_2.mat';
    gen_list='.\gene_list/MDD_Gene_List.csv';
    
    area1_name='ba10p_l_N10_nlin2Stdcolin27'; % identic+case sensitive to main_r.name
    area2_name='4008'; % identic+case sensitive to main_r.name
    gene_index=21;
    cutoff=[180,256;
        50,256;
        95,256];
    %%% DO NOT EDIT CODE BELOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'gui')
    [path_result_mat,gen_list,area1_name,area2_name,gene_index,marker1_index,marker2_index,cutoff,cutoff_mode,cutoff_margin]=convert_parameter(parameter);
    
    [vis_path,~,~] = fileparts(mfilename('fullpath'));
    reference_brain=[vis_path '\visualization_template\mni_icbm152_t1_masked.nii'];

    %reference_brain=load_nii(reference_brain);
    spm_header_reference_brain=spm_vol(reference_brain);
    spm_header_reference_origin=spm_header_reference_brain.mat(1:3,4);
    reference_brain = spm_read_vols(spm_header_reference_brain);
else
    error('No running mode selected. [man:manual; gui:call via config_gui]');
end

%%%% Load coords and zscores=> build vars
load(path_result_mat);
p_value_based_sample_allocation=1;
if p_value_based_sample_allocation==1
    % function pval_filter is saved in utilities and removes duplicate
    % entries of main_r only the TS with the higher pval will be kept
    [main_r, TS]=pval_filter(main_r);
end

% initialize vars
area1_coords=[];
area2_coords=[];
area1_coords_n=[];
area2_coords_n=[];
area1_zscores=[];
area2_zscores=[];
area1_specimen=[];
area2_specimen=[];
area1_area=[];
area2_area=[];
combined_coords=[];

%%%% get exp lvl to corresponding coordinates
[probe_id,gene_symbol,entrez_id] = read_gen_file(gen_list);
% create struct & groups
for i=1:size(main_r,2)
    if strcmp(main_r(i).name,area1_name)
        area1_zscores=[area1_zscores;main_r(i).validated_zscores];
        for j=1:size(main_r(i).validated_zscores,1)
            area1_coords_n{end+1}=main_r(i).data2plot{1,1}(:,j);
            area1_specimen{end+1}=main_r(i).specimen;
            area1_area{end+1}=main_r(i).name;
            combined_coords(:,end+1)=main_r(i).data2plot{1,1}(:,j);
        end
    else
        area2_zscores=[area2_zscores;main_r(i).validated_zscores];
        for j=1:size(main_r(i).validated_zscores,1)
            area2_coords_n{end+1}=main_r(i).data2plot{1,1}(:,j);
            area2_specimen{end+1}=main_r(i).specimen;
            area2_area{end+1}=main_r(i).name;
            combined_coords(:,end+1)=main_r(i).data2plot{1,1}(:,j);
        end
    end
end



spm_header_map=spm_vol(parameter.pmap{1,1});
spm_header_map_origin=spm_header_map.mat(1:3,4);

% shift because jubrain maps with 256x256x256 dim and mni152 193x229x193
% dim
shift_between_map_and_ref=abs(spm_header_reference_origin)-abs(spm_header_map_origin);
%shift_between_map_and_ref=[0 0 0];

% ref_origin=[reference_brain.hdr.hist.qoffset_x,reference_brain.hdr.hist.qoffset_y,reference_brain.hdr.hist.qoffset_z];
% 
% map=load_nii(parameter.pmap{1,1});
% map_origin=[map.hdr.hist.qoffset_x,map.hdr.hist.qoffset_y,map.hdr.hist.qoffset_z];
% 
% shift_between_map_and_ref=abs(map_origin)-abs(ref_origin);

%     a=[shift_between_map_and_ref(1),shift_between_map_and_ref(2),shift_between_map_and_ref(3)]; 
%     b=repmat(a,size(combined_coords,2),1);
%     
%     combined_coords=combined_coords+b';

clear ref_origin map map_origin


if cutoff_mode==2   % wird nur bei margin ausgeführt
    % margin_cutoff=0;
    % if margin_cutoff==1
    a=[shift_between_map_and_ref(1),shift_between_map_and_ref(2),shift_between_map_and_ref(3)];
    b=repmat(a,size(combined_coords,2),1);
    combined_coords=combined_coords+b';
    margin=cutoff_margin;
    
    margin_min_x=min(combined_coords(2,:))-margin(1);
    if margin_min_x<1
       margin_min_x=1;
    end
    cutoff(1,1)=margin_min_x;
    
    margin_max_x=max(combined_coords(2,:))+margin(1);
    if margin_max_x>size(reference_brain,2)
       margin_max_x=size(reference_brain,2);
    end
    cutoff(1,2)=margin_max_x;
    
    margin_min_y=min(combined_coords(1,:))-margin(2);
    if margin_min_y<1
       margin_min_y=1;
    end
    cutoff(2,1)=margin_min_y;
    
    margin_max_y=max(combined_coords(1,:))+margin(2);
    if margin_max_y>size(reference_brain,1)
       margin_max_y=size(reference_brain,1);
    end
    cutoff(2,2)=margin_max_y;
    
    margin_min_z=min(combined_coords(3,:))-margin(3);
    if margin_min_z<1
       margin_min_z=1;
    end
    cutoff(3,1)=margin_min_z;
    margin_max_z=max(combined_coords(3,:))+margin(3);
    if margin_max_z>size(reference_brain,3)
       margin_max_z=size(reference_brain,3);
    end
    cutoff(3,2)=margin_max_z;
    
    %end
end
h_figure=figure;





% slices to cut above
% nr_xslice=cutoff(1,1);
% nr_yslice=cutoff(2,1);
% nr_zslice=cutoff(3,1);
iso=65;
%iso=70;


cropped_colin=reference_brain;
cropped_colin(cutoff(2,1):cutoff(2,2),cutoff(1,1):cutoff(1,2),cutoff(3,1):cutoff(3,2))=0;

% Get the input voxel volume
data=struct;
data.volume=uint8(cropped_colin);

% Calculate min and max value to scale the intensities
data.vmax=double(max(data.volume(:)));
data.vmin=double(min(data.volume(:)));

% Get or Calculate the staritng iso value
data.iso=double(iso);
%data.iso=double(max(data.volume(:)));

% Default input values
data.scaling=1;
data.posx=0.5;
data.posy=0.5;
data.posz=0.5;
data.scales=[1 1 1];
data.color=[.95 .87 .73];

data.sizes=size(data.volume);
data.isorender=false;
data.slicerender=false;
data.h1=[];
data.h2=[];
data.h3=[];
% make downsampled volumes, for fast iso surface rendering.
data.small1=imresize3d(data.volume,0.25,[],'linear');
data.small2=imresize3d(data.volume,0.5,[],'linear');

% create isosurface
data=make_isosurface(data);
% plot isosurface to figure
data=show_isosurface(data);
xlabel('x');
ylabel('y');
zlabel('z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add orthoslice to figure
%do_cutoff=1;
if cutoff_mode==2 || cutoff_mode==3
    
    ref_dims=[size(reference_brain,1),size(reference_brain,2),size(reference_brain,3)];
    x_dim_ref=ref_dims(1);
    y_dim_ref=ref_dims(2);
    z_dim_ref=ref_dims(3);
    
    hold on
    % show non transparant X start slice (back to front)
    hs_x_s=showXSlice(reference_brain, cutoff(1,1));
%     slice(:,:)=reference_brain.img(:,cutoff(1,1),:);
%     slice(:,1:cutoff(3,1)-1)=0;
%     %slice(:,cutoff(3,2)+1:256)=0;
%     slice(:,cutoff(3,2)+1:ref_dims(3))=0;
%     slice(1:cutoff(2,1)-1,:)=0;
%     %slice(cutoff(2,2)+1:256,:)=0;
%     slice(cutoff(2,2)+1:ref_dims(2),:)=0;
    slice=hs_x_s.CData(:,:,1);
    slice(slice>0)=1;
    %if margin_cutoff==1
%     slice=rot90(slice);
%     %end
%     slice=flipud(slice);
    if nnz(slice)>1
        if islogical(slice)
            slice=imfill(slice,'holes');
        else
            slice=imfill(im2bw(slice,graythresh(slice)),'holes');
        end
    end
    set(hs_x_s,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',slice);
        % show non transparant X end slice (back to front)
    hs_x_e=showXSlice(reference_brain, cutoff(1,2));
%     slice(:,:)=reference_brain.img(:,cutoff(1,2),:);
%     slice(:,1:cutoff(3,1)-1)=0;
%     %slice(:,cutoff(3,2)+1:256)=0;
%     slice(:,cutoff(3,2)+1:ref_dims(3))=0;
%     slice(1:cutoff(2,1)-1,:)=0;
%     %slice(cutoff(2,2)+1:256,:)=0;
%     slice(cutoff(2,2)+1:ref_dims(2),:)=0;
slice=hs_x_e.CData(:,:,1);
    slice(slice>0)=1;
    %slice=rot90(slice);slice=flipud(slice);
    
    if nnz(slice)>1
        if islogical(slice)
            slice=imfill(slice,'holes');
        else
            slice=imfill(im2bw(slice,graythresh(slice)),'holes');
        end
    end
    set(hs_x_e,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',slice);
    
        % show non transparant Y start slice (back to front)
    hs_y_s=showYSlice(reference_brain, cutoff(2,1));
%     slice(:,:)=reference_brain.img(cutoff(2,1),:,:);
slice=hs_y_s.CData(:,:,1);
    slice(slice>0)=1;%slice=rot90(slice);slice=flipud(slice);
%     slice(:,1:cutoff(3,1)-1)=0;
%     %slice(:,cutoff(3,2)+1:256)=0;
%     slice(:,cutoff(3,2)+1:ref_dims(3))=0;
%     slice(1:cutoff(1,1)-1,:)=0;
%     %slice(cutoff(1,2)+1:256,:)=0;
%     slice(cutoff(1,2)+1:ref_dims(1),:)=0;
    if nnz(slice)>1
        if islogical(slice)
            slice=imfill(slice,'holes');
        else
            slice=imfill(im2bw(slice,graythresh(slice)),'holes');
        end
    end
    set(hs_y_s,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',slice);
        % show non transparant Y end slice (back to front)

    hs_y_e=showYSlice(reference_brain, cutoff(2,2));
    %slice(:,:)=reference_brain.img(cutoff(2,2),:,:);
    slice=hs_y_e.CData(:,:,1);
    slice(slice>0)=1;%slice=rot90(slice);slice=flipud(slice);
%     slice(:,1:cutoff(3,1)-1)=0;
%     %slice(:,cutoff(3,2)+1:256)=0;
%     slice(:,cutoff(3,2)+1:ref_dims(3))=0;
%     slice(1:cutoff(1,1)-1,:)=0;
%     %slice(cutoff(1,2)+1:256,:)=0;
%     slice(cutoff(1,2)+1:ref_dims(1),:)=0;
    
    
    if nnz(slice)>1
        if islogical(slice)
            slice=imfill(slice,'holes');
        else
            slice=imfill(im2bw(slice,graythresh(slice)),'holes');
        end
    end
    set(hs_y_e,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',slice);
        % show non transparant Z start slice (back to front)
    hs_z_s=showZSlice(reference_brain, cutoff(3,1));
    %slice(:,:)=reference_brain.img(:,:,cutoff(3,1));
   slice=hs_z_s.CData(:,:,1);
    slice(slice>0)=1;%slice=rot90(slice);slice=flipud(slice);
%     slice(:,1:cutoff(1,1)-1)=0;
%     %slice(:,cutoff(1,2)+1:256)=0;
%     slice(:,cutoff(1,2)+1:ref_dims(1))=0;
%     slice(1:cutoff(2,1)-1,:)=0;
%     %slice(cutoff(2,2)+1:256,:)=0;
%     slice(:,cutoff(1,2)+1:ref_dims(1))=0;
    if nnz(slice)>1
        if islogical(slice)
            slice=imfill(slice,'holes');
        else
            slice=imfill(im2bw(slice,graythresh(slice)),'holes');
        end
    end
    set(hs_z_s,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',slice);
    
    % show non transparant Z end slice (back to front)
    hs_z_e=showZSlice(reference_brain, cutoff(3,2));
    slice=hs_z_e.CData(:,:,1);
%     slice(:,:)=reference_brain.img(:,:,cutoff(3,2));
    slice(slice>0)=1;%slice=rot90(slice);slice=flipud(slice);
%     slice(:,1:cutoff(1,1)-1)=0;
%     %slice(:,cutoff(1,2)+1:256)=0;
%     slice(:,cutoff(1,2)+1:ref_dims(1))=0;
%     slice(1:cutoff(2,1)-1,:)=0;
%     %slice(cutoff(2,2)+1:256,:)=0;
%     slice(cutoff(2,2)+1:ref_dims(2),:)=0;
    
    if nnz(slice)>1
        if islogical(slice)
            slice=imfill(slice,'holes');
        else
            slice=imfill(im2bw(slice,graythresh(slice)),'holes');
        end
    end
    set(hs_z_e,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',slice);
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END add orthoslice to figure


clc;

factor_specimen=[area1_specimen area2_specimen];
factor_area=[area1_area area2_area];
combined_zscores=[area1_zscores;area2_zscores];

[ combined_zscores,area1_zscores,area2_zscores,unique_entrez_id  ] = switch2gensymbol_lvl( entrez_id,combined_zscores,size(area1_zscores,1),size(area2_zscores,1) );


area2_zscores_searched_gene=area2_zscores(:,gene_index);
area1_zscores_searched_gene=area1_zscores(:,gene_index);

szc=[area2_zscores_searched_gene;area1_zscores_searched_gene];
min_szc=sort(szc,'ascend');min_szc=min_szc(3);
%min_szc=0.13;
max_szc=sort(szc,'descend');max_szc=max_szc(3);
%max_szc=0.65;

%%%%%%%%%%%%%%%%%%%%  Min Max für zwei abbildung fest und gleich setzen
% min_szc=0.044754;
% max_szc=0.88142;
%%%%%%%%%%%%%%%%%%%%%  ende festsetzen min max

colorline = linspace(min_szc,max_szc,256);
%disp('test');
counter_area1=1;
counter_area2=1;
%load('redgreenmaps.mat');
% load('orange_white_green_cmap.mat');
% colmap=orange_white_green_cmap;
load('magenta_white_green_colormap.mat');
colmap=magenta_white_green_colormap;

a=[shift_between_map_and_ref(1),shift_between_map_and_ref(2),shift_between_map_and_ref(3)]; 

for i=1:size(main_r,2)
    if strcmp(main_r(i).name,area1_name)   % es soll der Marker für area1 geplottet werden
        num_coords= size(main_r(i).data2plot{1,1},2);
        for j=1:num_coords
            trip=main_r(i).data2plot{1,1}(:,j);
            area1_coords=[area1_coords,trip];
            trip=double(trip);
            trip=trip+a';
            value2color=area1_zscores_searched_gene(counter_area1);
            [c index] = min(abs(colorline-value2color));
            cube_color=colmap(index,:);
            %cube_color=[1 0 0];
            plot_marker(marker1_index,trip,cube_color,1);
            counter_area1=counter_area1+1;
            value2color=[];
        end
    else    % es soll der Marker für area2 geplottet werden
        num_coords= size(main_r(i).data2plot{1,1},2);
        for j=1:num_coords
            trip=main_r(i).data2plot{1,1}(:,j);
            area2_coords=[area2_coords,trip];
            trip=double(trip);
            trip=trip+a';
            value2color=area2_zscores_searched_gene(counter_area2);
            [c index] = min(abs(colorline-value2color));
            cube_color=colmap(index,:);
            %cube_color=[1 0 0];
            plot_marker(marker2_index,trip,cube_color,1);
            counter_area2=counter_area2+1;
            value2color=[];
        end
    end
end

% pattern = '.*(\d).mat';  % old code to extract thresh...
% [tok, mat]= regexp(parameter.extracted_data_file,pattern,'tokens','match');
% thresh=cell2mat(tok{1,1}{1,1});

%%% load maps, add maps transparently to img
%a=[shift_between_map_and_ref(2),shift_between_map_and_ref(1),shift_between_map_and_ref(3)]; 
%     spm_header_reference_brain=spm_vol(reference_brain);
%     spm_header_reference_origin=spm_header_reference_brain.mat(1:3,4);
%     reference_brain = spm_read_vols(spm_header_reference_brain);
%a=shift_between_map_and_ref';
% map against ABA macro label, plot one map
%a=[-16 -32 -32]

a=[shift_between_map_and_ref(2),shift_between_map_and_ref(1),shift_between_map_and_ref(3)]; 


voi1_header=spm_vol(parameter.pmap{1,1});
voi1=spm_read_vols(voi1_header);
%threshold anwenden
voi1(voi1<(str2num(parameter.map_threshold{:})/10))=0;
patch_hdl_fp1=plot_area_as_surface(voi1,'red',0.2);
b=repmat(a,size(patch_hdl_fp1.Vertices,1),1);
patch_hdl_fp1.Vertices=patch_hdl_fp1.Vertices+b;



if parameter.search_mode{:}==1  % map againstsecond map, so plot both
%     fp2=load_nii(parameter.pmap{1,2});
%     fp2=fp2.img;
    
    spm_readed_header=spm_vol(parameter.pmap{1,2});
    fp2=spm_read_vols(spm_readed_header);
    
    
    %threshold anwenden
    fp2(fp2<(str2num(parameter.map_threshold{:})/10))=0;

    patch_hdl_fp2=plot_area_as_surface(fp2,'blue',0.2);
    b=repmat(a,size(patch_hdl_fp2.Vertices,1),1);
    patch_hdl_fp2.Vertices=patch_hdl_fp2.Vertices+b;
end

drawnow


colormap(colmap)
% cb = colorbar('vert');
% zlab = get(cb,'ylabel');
% set(zlab,'String','Exp. Level [z scores]');

% auskommentieren für colorbar

h = colorbar('vert');
if min_szc>max_szc
    min_szc_tmp=min_szc;
    max_szc_tmp=max_szc;
    max_szc=min_szc;
    min_szc=max_szc_tmp;
end
set(gca, 'CLim', [min_szc, max_szc])
set(h, 'XTick', [min_szc, max_szc])
set(h,'XTickLabel',{num2str(min_szc) ,num2str(max_szc)},'FontSize',12,'FontWeight','bold', 'XColor', 'w') %# don't add units here...
xlabel(h, 'zscores')                                  %# ...use xlabel to add units


scrsz = get(0,'ScreenSize'); %getting the screensize of the 1 screen
set(gcf,'Position',[1 0 1000 1000],'MenuBar','figur','ToolBar','figur','resize','on','Color',[0 0 0]) % fullscreen
set(gca, 'XColor', 'w', 'YColor', 'w','ZColor', 'w','Color',[0 0 0])
%set(gcf,'Position',[1 0 200 200],'MenuBar','none','ToolBar','none','resize','off') % fullscreen
set(gca,'XDir','reverse');
axis tight

% view für abbildung aus abstrakt
%view(-124.9395,24.9668);

% ba10 figures
%view(-122,14);
%zoom(2.5)
%lightsource=camlight(-122,20);
lightsource=camlight('headlight');

% nina figures
% view(-35,20);
% zoom(2.5)
% lightsource=camlight(10,200);

%%%% Beschriftungen und Colorbar ausschalten
% axis off
%legend('hide')
%colorbar('hide')



[probe_id,gene_symbol,entrez_id] = read_gen_file(gen_list);
%gene_symbol=unique(gene_symbol);
unique_entrez_id=unique(entrez_id);
gene_symbol_according2_unique_entrezid = cell(size(unique_entrez_id,1),1);
for i=1:size(unique_entrez_id,1)
    tmp=find(entrez_id==unique_entrez_id(i));
    gene_symbol_according2_unique_entrezid(i)=unique(gene_symbol(tmp));
end
ugs= gene_symbol_according2_unique_entrezid;
% ugs=unique(gene_symbol);
strn = regexprep(ugs{gene_index},'"','');
set(gcf,'NumberTitle','off');
string=[area1_name '_vs_' area2_name '_Gene_' strn];
set(gcf,'Name',string);
h_curfig=gcf;

createtextbox(gcf,strrep(string,'_','\_'));




% save_str=['.' filesep 'output' filesep 'img' filesep area1_name '_vs_' area2_name '_Gene_' strn '_th_' thresh '_without_Map.png'];
% w = warning ('off','all');
% export_fig(save_str);
% w = warning ('on','all');
%save_str_p=['.' filesep 'output' filesep 'img' filesep area1_name '_vs_' area2_name '_Gene_' strn '_without_Map_print'];
%print(h_curfig,save_str_p,'-dtiffn','-r300','-opengl');



%ugs=unique(gene_symbol);
strn = regexprep(ugs{gene_index},'"','');
set(gcf,'NumberTitle','off');
string=[area1_name '_vs_' area2_name '_Gene_' strn];
set(gcf,'Name',string);
writetable(TS,['.' filesep 'output' filesep 'img' filesep string '_TS_data.txt']);
save_str=['.' filesep 'output' filesep 'img' filesep area1_name '_vs_' area2_name '_Gene_' strn '_th_' parameter.map_threshold{:} '.fig'];
savefig(save_str);
disp(['>>>>File saved to ' save_str]);
save_str=['.' filesep 'output' filesep 'img' filesep area1_name '_vs_' area2_name '_Gene_' strn '_th_' parameter.map_threshold{:} '_with_Map.png'];
% w = warning ('off','all');
% export_fig(save_str);
% w = warning ('on','all');

print(gcf,save_str,'-dpng')


toc;

TS
%create_movie
end

function [path_result_mat,gen_list,area1_name,area2_name,gene_index,marker1_index,marker2_index,cutoff,cutoff_mode,cutoff_margin]=convert_parameter(parameter)

path_result_mat=parameter.extracted_data_file{1};
gen_list=parameter.gene_list{1};
area1_name=parameter.name{1}; % identic+case sensitive to main_r.name
area2_name=parameter.name{2}; % identic+case sensitive to main_r.name
%reference_brain=parameter.reference_brain{1};
gene_index=parameter.gene_index(1);
marker1_index=parameter.marker1_index(1);
marker2_index=parameter.marker2_index(1);
cutoff=parameter.cutoff;
cutoff_mode=parameter.cutoff_mode;
cutoff_margin=parameter.cutoff_margin;
end

function createtextbox(figure1,string)
%CREATETEXTBOX(FIGURE1)
%  FIGURE1:  annotation figure

%  Auto-generated by MATLAB on 12-Feb-2015 10:40:51

% Create textbox
annotation(figure1,'textbox',[0.3046875 0.867 0.3765625 0.0709999999999998],...
    'Color',[1 1 1],...
    'String',string,...
    'HorizontalAlignment','center',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FitBoxToText','off');
end

function plot_marker(marker,trip,cube_color,alpha)
    switch marker
        case 2
            bubbleplot3(trip(2),trip(1),trip(3),2,cube_color,alpha);
        case 3
            %%%%%% orig und funktioniereder cube plot
            %plotcube([2 2 2],[trip(2) trip(1) trip(3)],.7,cube_color);
            %%%%%% aus 7 cubes ein plus zu basteln
            plotcube([1 1 1],[trip(2) trip(1) trip(3)],alpha,.2,cube_color);
           
            plotcube([1 1 1],[trip(2)-1 trip(1) trip(3)],alpha,.2,cube_color);
            plotcube([1 1 1],[trip(2)+1 trip(1) trip(3)],alpha,.2,cube_color);
            
            plotcube([1 1 1],[trip(2) trip(1)-1 trip(3)],alpha,.2,cube_color);
            plotcube([1 1 1],[trip(2) trip(1)+1 trip(3)],alpha,.2,cube_color);
            
            plotcube([1 1 1],[trip(2) trip(1) trip(3)-1],alpha,.2,cube_color);
            plotcube([1 1 1],[trip(2) trip(1) trip(3)+1],alpha,.2,cube_color);
        case 4
            hold on
            scatter3(trip(2),trip(1),trip(3),'o','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',1);
            hold off
        case 5
            hold on
            h_scatter=scatter3(trip(2),trip(1),trip(3),'+','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',5);
            set (h_scatter,'SizeData',200);
            hold off
        case 6
            hold on
            h_scatter=scatter3(trip(2),trip(1),trip(3),'*','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',2);
            set (h_scatter,'SizeData',100);
            hold off
        case 7
            hold on
            scatter3(trip(2),trip(1),trip(3),'.','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',1);
            hold off
        case 8
            hold on
            scatter3(trip(2),trip(1),trip(3),'x','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',1);
            hold off
        case 9
            hold on
            scatter3(trip(2),trip(1),trip(3),'s','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',1);
            hold off
        case 10
            hold on
            scatter3(trip(2),trip(1),trip(3),'d','MarkerFaceColor',cube_color,'MarkerEdgeColor',cube_color,'LineWidth',1);
            hold off            
    end
end