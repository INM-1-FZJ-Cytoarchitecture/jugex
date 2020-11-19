function varargout = Visualization(varargin)
% VISUALIZATION MATLAB code for Visualization.fig
%      VISUALIZATION, by itself, creates a new VISUALIZATION or raises the existing
%      singleton*.
%
%      H = VISUALIZATION returns the handle to a new VISUALIZATION or the handle to
%      the existing singleton*.
%
%      VISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZATION.M with the given input arguments.
%
%      VISUALIZATION('Property','Value',...) creates a new VISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Visualization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Visualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Visualization

% Last Modified by GUIDE v2.5 04-Oct-2016 16:52:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Visualization_OpeningFcn, ...
    'gui_OutputFcn',  @Visualization_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Visualization is made visible.
function Visualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Visualization (see VARARGIN)

%uistack(handles.checkbox_cutoff,'top');

% Choose default command line output for Visualization
handles.output = hObject;

set(gcf,'name','Visualization','numbertitle','off')

set(handles.slider_x,'Visible','Off');
set(handles.slider_y,'Visible','Off');
set(handles.slider_z,'Visible','Off');

set(handles.Sslider,'Visible','Off');
% Update handles structure
guidata(hObject, handles);

uistack(handles.none,'top');





% UIWAIT makes Visualization wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%createColormapDropdown (gcf,gca);


% --- Outputs from this function are returned to the command line.
function varargout = Visualization_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_select_reference_brain_Callback(hObject, eventdata, handles)
% hObject    handle to text_select_reference_brain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_select_reference_brain as text
%        str2double(get(hObject,'String')) returns contents of text_select_reference_brain as a double


% --- Executes during object creation, after setting all properties.
function text_select_reference_brain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_select_reference_brain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

p = mfilename('fullpath');
parts = strsplit(p, filesep);
str='';
for i=1:size(parts,2)-1
    str=[str parts{i} filesep];
end
str=[str 'mask' filesep 'colin27T1_for_mapping.nii'];
myhandles = guihandles(gcf);
set(myhandles.text_select_reference_brain,'string','colin27T1_for_mapping.nii');
set(myhandles.text_select_reference_brain,'TooltipString',str);


% --- Executes on button press in pushbutton_select_reference_brain.
function pushbutton_select_reference_brain_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_reference_brain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.nii'},'Select reference brain (*.nii file)');
f = fullfile(pathname,filename);
parts = strsplit(pathname, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep filename];
myhandles = guihandles(gcf);
set(myhandles.text_select_reference_brain,'string',disp_str);
set(myhandles.text_select_reference_brain,'TooltipString',f);


function edit_extracted_data_file_Callback(hObject, eventdata, handles)
% hObject    handle to text_extracted_data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_extracted_data_file as text
%        str2double(get(hObject,'String')) returns contents of text_extracted_data_file as a double


% --- Executes during object creation, after setting all properties.
function text_extracted_data_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_extracted_data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_change_input_file.
function pushbutton_change_input_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_change_input_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

myhandles = guihandles(gcf);
[filename, pathname] = uigetfile({'*_parameter.mat'},'Select file with extracted data',['.' filesep 'output' filesep 'extracted_data' filesep]);
f = fullfile(pathname,filename);
handles.path_parameter_file=f;
guidata(gcf, myhandles);
guidata(gcf, handles);
load(f);
output_parameter=parameter;
clear parameter
% ab hier liegt output_parameter struct vor  Ausgabe in
% myhandles.text_extracted_data_file,'string', anpassen

full_str=[pathname output_parameter.project_name{1,1} '_th_' num2str(output_parameter.map_threshold{1,1}) '.mat'];

parts = strsplit(pathname, filesep);
DirPart = parts{end};
disp_str=['.'  parts{end-2} filesep parts{end-1} filesep output_parameter.project_name{1,1} '_th_' num2str(output_parameter.map_threshold{1,1}) '.mat'];


% myhandles = guihandles(gcf);
set(myhandles.text_extracted_data_file,'string',disp_str);
set(myhandles.text_extracted_data_file,'TooltipString',full_str);

if output_parameter.search_mode{1}==1
    set(myhandles.text_name_area1,'string',output_parameter.name{1});
    set(myhandles.text_name_area2,'string',output_parameter.name{2});
elseif output_parameter.search_mode{1}==2
    set(myhandles.text_name_area1,'string',output_parameter.name{1});
    ontology = readtable('tmp_ontology.csv');
    ind_structure=find(ontology.id==output_parameter.struct_id{1});
    set(myhandles.text_name_area2,'string',ontology.name{ind_structure});
end

f = output_parameter.gene_list{1,1};
% parts = strsplit(pathname, filesep);
% DirPart = parts{end};
% disp_str=['.' filesep DirPart filesep filename];
% myhandles = guihandles(gcf);
% set(myhandles.text_select_gene_list,'string',disp_str);
% set(myhandles.text_select_gene_list,'TooltipString',f);

[probe_id,gene_symbol,entrez_id] = read_gen_file(f);
%gene_symbol=unique(gene_symbol);
unique_entrez_id=unique(entrez_id);
gene_symbol_according2_unique_entrezid = cell(size(unique_entrez_id,1),1);
for i=1:size(unique_entrez_id,1)
    tmp=find(entrez_id==unique_entrez_id(i));
    gene_symbol_according2_unique_entrezid(i)=unique(gene_symbol(tmp));
end
set(handles.popupmenu_select_gene,'string',gene_symbol_according2_unique_entrezid); %%% ich brauche die unique id bzw den index der id des ausgewählten gene_symbol
guidata(hObject, handles);


function edit_select_gene_list_Callback(hObject, eventdata, handles)
% hObject    handle to text_select_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_select_gene_list as text
%        str2double(get(hObject,'String')) returns contents of text_select_gene_list as a double


% --- Executes during object creation, after setting all properties.
function text_select_gene_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_select_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_select_gene_list.
function pushbutton_select_gene_list_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.csv'},'Select gene list file',['.' filesep 'gene_list' filesep]);
f = fullfile(pathname,filename);
parts = strsplit(pathname, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep filename];
myhandles = guihandles(gcf);
set(myhandles.text_select_gene_list,'string',disp_str);
set(myhandles.text_select_gene_list,'TooltipString',f);

[probe_id,gene_symbol,entrez_id] = read_gen_file(f);
%gene_symbol=unique(gene_symbol);
unique_entrez_id=unique(entrez_id);
gene_symbol_according2_unique_entrezid = cell(size(unique_entrez_id,1),1);
for i=1:size(unique_entrez_id,1)
    tmp=find(entrez_id==unique_entrez_id(i));
    gene_symbol_according2_unique_entrezid(i)=unique(gene_symbol(tmp));
end
set(handles.popupmenu_select_gene,'string',gene_symbol_according2_unique_entrezid); %%% ich brauche die unique id bzw den index der id des ausgewählten gene_symbol
guidata(hObject, handles);




% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_select_gene.
function popupmenu_select_gene_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_gene contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_gene
% contents = cellstr(get(hObject,'String')); %returns popupmenu_select_gene contents as cell array
% contents{get(hObject,'Value')}; %returns selected item from popupmenu_select_gene
%
% myhandles = guihandles(gcf);
% f=get(myhandles.text_select_gene_list,'TooltipString');
%
% [probe_id,gene_symbol,entrez_id] = read_gen_file(f);
% gene_symbol=unique(gene_symbol);
% entrez_id=unique(entrez_id);
%
% find(contents,contents{get(hObject,'Value')});

% --- Executes during object creation, after setting all properties.
function popupmenu_select_gene_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_start_visualization.
function pushbutton_start_visualization_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_visualization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


load(handles.path_parameter_file);
output_parameter=parameter;
clear parameter
myhandles = guihandles(gcf);
output_parameter.name={get(myhandles.text_name_area1,'String'),get(myhandles.text_name_area2,'String')};
%output_parameter.gene_list={get(myhandles.text_select_gene_list,'TooltipString')};
output_parameter.extracted_data_file={get(myhandles.text_extracted_data_file,'TooltipString')};
contents = cellstr(get(myhandles.popupmenu_select_gene,'String')); %returns popupmenu_select_gene contents as cell array
%contents{get(myhandles.popupmenu_select_gene,'Value')}; %returns selected item from popupmenu_select_gene
output_parameter.gene_index=find(ismember(contents,contents{get(myhandles.popupmenu_select_gene,'Value')}));

contents = cellstr(get(myhandles.popupmenu_select_marker_area1,'String')); %returns popupmenu_select_gene contents as cell array
%contents{get(myhandles.popupmenu_select_gene,'Value')}; %returns selected item from popupmenu_select_gene
output_parameter.marker1_index=find(ismember(contents,contents{get(myhandles.popupmenu_select_marker_area1,'Value')}));
contents = cellstr(get(myhandles.popupmenu_select_marker_area2,'String')); %returns popupmenu_select_gene contents as cell array
%contents{get(myhandles.popupmenu_select_gene,'Value')}; %returns selected item from popupmenu_select_gene
output_parameter.marker2_index=find(ismember(contents,contents{get(myhandles.popupmenu_select_marker_area2,'Value')}));


%output_parameter.reference_brain={get(myhandles.text_select_reference_brain,'TooltipString')};
if get(myhandles.radiobutton_none,'Value')==1
    output_parameter.cutoff=[255,256;255,256;255,256];
    output_parameter.cutoff_mode=1;
    output_parameter.cutoff_margin=[];
elseif get(myhandles.radiobutton_margin,'Value')==1
    output_parameter.cutoff=[255,256;255,256;255,256];
    output_parameter.cutoff_margin=[str2num(get(myhandles.x_margin,'String')),str2num(get(myhandles.y_margin,'String')),str2num(get(myhandles.z_margin,'String'))];
    output_parameter.cutoff_mode=2;
elseif get(myhandles.radiobutton_free,'Value')==1
    output_parameter.cutoff=[str2num(get(myhandles.x_start,'String')),str2num(get(myhandles.x_end,'String'));...
        str2num(get(myhandles.y_start,'String')),str2num(get(myhandles.y_end,'String'));...
        str2num(get(myhandles.z_start,'String')),str2num(get(myhandles.z_end,'String'))];
    output_parameter.cutoff_mode=3;
    output_parameter.cutoff_margin=[];
end

create_surface('gui',output_parameter);


% --- Executes on selection change in popupmenu_select_marker_area1.
function popupmenu_select_marker_area1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_marker_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_marker_area1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_marker_area1


% --- Executes during object creation, after setting all properties.
function popupmenu_select_marker_area1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_marker_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_select_marker_area2.
function popupmenu_select_marker_area2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_marker_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_marker_area2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_marker_area2


% --- Executes during object creation, after setting all properties.
function popupmenu_select_marker_area2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_marker_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function uipanel_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    [vis_path,~,~] = fileparts(mfilename('fullpath'));
    % change v1.1 Sebastian
    reference_brain_str=[vis_path '\visualize_AB_exp_Data\visualization_template\mni_icbm152_t1_masked.nii'];
    %reference_brain='D:\matlab_files\matlab_add_ons\Toolboxes\JuGex\tmp_2_del_mask\colin27T1_for_mapping.nii';
    %reference_brain=load_nii(reference_brain);
    spm_readed_header=spm_vol(reference_brain_str);
    reference_brain.img=spm_read_vols(spm_readed_header);
    % end change 1.1 Sebastian
    ref_dims=[size(reference_brain.img,1),size(reference_brain.img,2),size(reference_brain.img,3)];
    x_dim_ref=ref_dims(2);
    y_dim_ref=ref_dims(1);
    z_dim_ref=ref_dims(3);
min_start=1;
%max_start=256;
%update_box='x';
%ssx=superSlider('parent',findobj('Tag','uipanel_cutoff'),'numSlides', 2,'Min',0,'Max',256,'tagName','slider_x','position',[.07 .265 .35 .05],'stepsize',8,'value',[10 240],'discrete',true);
uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','text','String','X:','Position', [1 78 20 15]);
h_start_x=uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','edit','Tag','x_start','String',num2str(min_start),'Position', [250 77 40 20]);
h_end_x=uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','edit','Tag','x_end','String',num2str(x_dim_ref),'Position', [295 77 40 20]);
%para=struct('handle_start',h_start_x,'handle_end',h_end_x,'slider_tag','slider_x');
ssx=superSlider('parent',findobj('Tag','uipanel_cutoff'),'numSlides', 2,'Min',1,'Max',x_dim_ref,'tagName','slider_x','position',[.07 .33 .35 .05],'stepsize',8,'value',[min_start x_dim_ref],'discrete',true','callback',@update_edit_box);
%set(ssx,'callback','{@update_edit_box,parawww}');
% update_box='y';
ssy=superSlider('parent',findobj('Tag','uipanel_cutoff'),'numSlides', 2,'Min',1,'Max',y_dim_ref,'tagName','slider_y','position',[.07 .215 .35 .05],'stepsize',8,'value',[min_start y_dim_ref],'discrete',true,'callback',@update_edit_box);
uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','text','String','Y:','Position', [1 42 20 15]);
h_start_y=uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','edit','Tag','y_start','String',num2str(min_start),'Position', [250 41 40 20]);
h_end_y=uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','edit','Tag','y_end','String',num2str(y_dim_ref),'Position', [295 41 40 20]);
% set(ssy,'callback','{@update_edit_box,h_start_y,h_end_y,ssy}');
% update_box='z';
ssz=superSlider('parent',findobj('Tag','uipanel_cutoff'),'numSlides', 2,'Min',1,'Max',z_dim_ref,'tagName','slider_z','position',[.07 .09 .35 .05],'stepsize',8,'value',[min_start z_dim_ref],'discrete',true,'callback',@update_edit_box);
uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','text','String','Z:','Position', [1 6 20 15]);
h_start_z=uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','edit','Tag','z_start','String',num2str(min_start),'Position', [250 5 40 20]);
h_end_z=uicontrol('parent',findobj('Tag','uipanel_cutoff'),'Style','edit','Tag','z_end','String',num2str(z_dim_ref),'Position', [295 5 40 20]);
% set(ssz,'callback','{@update_edit_box,h_start_z,h_end_z,ssz}');
%disp('...');





function update_edit_box(test,h_object)
slider_tag=get(test,'Tag');
parts = strsplit(slider_tag,'_');
slider_orientation=parts{end};

h_start_edit=findobj('Tag',[slider_orientation '_start']);
h_end_edit=findobj('Tag',[slider_orientation '_end']);

infoMatrix = get(test, 'UserData');

set(h_start_edit,'String',infoMatrix(1,1));
set(h_end_edit,'String',infoMatrix(1,2));


% --- Executes on button press in checkbox_cutoff.
function checkbox_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cutoff

if get(hObject,'Value')==1
    set(handles.cutoff_margin,'Visible','On');
    set(handles.cutoff_margin,'Parent',gcf);
    uistack(handles.cutoff_margin,'top');
else
    set(handles.cutoff_margin,'Visible','Off');
end

function z_margin_Callback(hObject, eventdata, handles)
% hObject    handle to z_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_margin as text
%        str2double(get(hObject,'String')) returns contents of z_margin as a double


% --- Executes during object creation, after setting all properties.
function z_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_margin_Callback(hObject, eventdata, handles)
% hObject    handle to y_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_margin as text
%        str2double(get(hObject,'String')) returns contents of y_margin as a double


% --- Executes during object creation, after setting all properties.
function y_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_margin_Callback(hObject, eventdata, handles)
% hObject    handle to x_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_margin as text
%        str2double(get(hObject,'String')) returns contents of x_margin as a double


% --- Executes during object creation, after setting all properties.
function x_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_none.
function radiobutton_none_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_none
if get(hObject,'Value')==1
    set(handles.slider_x,'Visible','Off');
    set(handles.slider_y,'Visible','Off');
    set(handles.slider_z,'Visible','Off');
    set(handles.Sslider,'Visible','Off');
    set(handles.none,'Visible','On');
    uistack(handles.none,'top');
    set(handles.cutoff_margin,'Visible','Off');
    set(handles.uipanel_cutoff,'Visible','Off');
    %set(handles.radiobutton_none,'Value',0);
    set(handles.radiobutton_margin,'Value',0);
    set(handles.radiobutton_free,'Value',0);
    set(handles.none,'Parent',handles.uipanel1);
    
    %set(handles.cutoff_margin,'Visible','Off');
end


% --- Executes on button press in radiobutton_margin.
function radiobutton_margin_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_margin
if get(hObject,'Value')==1
    set(handles.slider_x,'Visible','Off');
    set(handles.slider_y,'Visible','Off');
    set(handles.slider_z,'Visible','Off');
    set(handles.Sslider,'Visible','Off');
    set(handles.none,'Visible','Off');
    set(handles.uipanel_cutoff,'Visible','Off');
    set(handles.radiobutton_none,'Value',0);
    %set(handles.radiobutton_margin,'Value',0);
    set(handles.radiobutton_free,'Value',0);
    set(handles.cutoff_margin,'Visible','On');
    set(handles.cutoff_margin,'Parent',handles.uipanel1);
    uistack(handles.cutoff_margin,'top');
end



% --- Executes on button press in radiobutton_free.
function radiobutton_free_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_free (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_free


if get(hObject,'Value')==1
    set(handles.slider_x,'Visible','On');
    set(handles.slider_y,'Visible','On');
    set(handles.slider_z,'Visible','On');
    set(handles.Sslider,'Visible','On');
    set(handles.cutoff_margin,'Visible','Off');
    set(handles.none,'Visible','Off');
    set(handles.radiobutton_none,'Value',0);
    set(handles.radiobutton_margin,'Value',0);
    set(handles.uipanel_cutoff,'Visible','On');
    %set(handles.radiobutton_free,'Value',0);
    %uipanel_cutoff_CreateFcn(hObject, eventdata, handles)
    set(handles.uipanel_cutoff,'Parent',handles.uipanel1);
    uistack(handles.uipanel_cutoff,'top');
end


% --- Executes during object creation, after setting all properties.
function none_CreateFcn(hObject, eventdata, handles)
% hObject    handle to none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%set(handles.none,'Visible','On');
%uistack(handles.none,'top');
