function varargout = Configuration(varargin)
% CONFIGURATION MATLAB code for Configuration.fig
%      CONFIGURATION, by itself, creates a new CONFIGURATION or raises the existing
%      singleton*.
%
%      H = CONFIGURATION returns the handle to a new CONFIGURATION or the handle to
%      the existing singleton*.
%
%      CONFIGURATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONFIGURATION.M with the given input arguments.
%
%      CONFIGURATION('Property','Value',...) creates a new CONFIGURATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Configuration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Configuration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Configuration

% Last Modified by GUIDE v2.5 05-Jul-2016 11:43:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Configuration_OpeningFcn, ...
                   'gui_OutputFcn',  @Configuration_OutputFcn, ...
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


% --- Executes just before Configuration is made visible.
function Configuration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Configuration (see VARARGIN)

% Choose default command line output for Configuration
set(gcf,'name','Configuration','numbertitle','off');

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Configuration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Configuration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_start_analysis.
function pushbutton_start_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
myhandles = guihandles(gcf);
output_parameter.pmap={get(myhandles.text_area1,'TooltipString'),get(myhandles.text_area2,'TooltipString')};
output_parameter.name={get(myhandles.edit_name_area1,'String'),get(myhandles.edit_name_area2,'String')};
output_parameter.gene_list={get(myhandles.text_gene_list,'TooltipString')};
output_parameter.output_folder={get(myhandles.text_output_folder,'TooltipString')};
output_parameter.project_name={get(myhandles.edit_project_name,'String')};
output_parameter.map_threshold={get(myhandles.edit_map_threshold,'String')};
if get(myhandles.checkbox_vs_structure,'Value')==1
    search_mode=2;
    ontology = readtable('tmp_ontology.csv');
    struct_id=ontology.id(get(myhandles.popupmenu_select_structure,'Value')-1);
    output_parameter.struct_id={struct_id};
elseif get(myhandles.checkbox_vs_structure,'Value')==0
    search_mode=1;
    struct_id=1; %please select
    output_parameter.struct_id={struct_id};
end
output_parameter.search_mode={search_mode};



strings = get(myhandles.listbox_donors, 'string');
curval = get(myhandles.listbox_donors, 'value');
for i=1:size(curval,2)
    donors{i}=strings{curval(i)};
end
output_parameter.donors=donors;
save('latest_settings.mat','output_parameter');
extract_exp_lvl('gui',output_parameter,2);% last parameter: 0, no output during download, 1, print api urls, 2, print tissueblock details

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to text_output_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_output_folder as text
%        str2double(get(hObject,'String')) returns contents of text_output_folder as a double


% --- Executes during object creation, after setting all properties.
function text_output_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_output_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if isdir('output') ~= 0
myhandles = guihandles(gcf);
set(myhandles.text_output_folder,'string',['.' filesep 'output' filesep]);
findpath = what('output');
set(myhandles.text_output_folder,'TooltipString',findpath.path);
end

% --- Executes on button press in pushbutton_select_output_folder.
function pushbutton_select_output_folder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_output_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[pathname] =uigetdir(['.' filesep 'output' filesep],'Select output folder');
parts = strsplit(pathname, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep];
f = fullfile(pathname);
myhandles = guihandles(gcf);
set(myhandles.text_output_folder,'string',disp_str);
set(myhandles.text_output_folder,'TooltipString',f);


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to text_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_gene_list as text
%        str2double(get(hObject,'String')) returns contents of text_gene_list as a double


% --- Executes during object creation, after setting all properties.
function text_gene_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_gene_list (see GCBO)
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
set(myhandles.text_gene_list,'string',disp_str);
set(myhandles.text_gene_list,'TooltipString',f);


% --- Executes on selection change in listbox_donors.
function listbox_donors_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_donors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_donors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_donors


% --- Executes during object creation, after setting all properties.
function listbox_donors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_donors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% %liest alle Donors für MicroArray ein
% s=xml2struct('http://api.brain-map.org/api/v2/data/Donor/query.xml?include=age&criteria=products[id$eq2]');
% 
% % s.Response.donors.donor{1,i}.name.Text
% % s.Response.donors.donor{1,i}.id.Text
% % size(s.Response.donors.donor,2)
% donors={};
% for i=1:size(s.Response.donors.donor,2)
% donors(i).name=s.Response.donors.donor{1,i}.name.Text;
% end
% 
% set(hObject,'string',{donors.name});
% guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to text_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_area1 as text
%        str2double(get(hObject,'String')) returns contents of text_area1 as a double


% --- Executes during object creation, after setting all properties.
function text_area1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to text_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_area2 as text
%        str2double(get(hObject,'String')) returns contents of text_area2 as a double


% --- Executes during object creation, after setting all properties.
function text_area2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_select_area1.
function pushbutton_select_area1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.hdr';'*.nii'},'Select first map',['.' filesep 'maps' filesep]);
f = fullfile(pathname,filename);
[pathstr,name,ext] = fileparts(filename);
parts = strsplit(pathname, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep filename];
myhandles = guihandles(gcf);
set(myhandles.text_area1,'string',disp_str);
set(myhandles.edit_name_area1,'string',name);
set(myhandles.text_area1,'TooltipString',f);

% --- Executes on button press in pushbutton_select_area2.
function pushbutton_select_area2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.hdr';'*.nii'},'Select second map',['.' filesep 'maps' filesep]);
f = fullfile(pathname,filename);
[pathstr,name,ext] = fileparts(filename);
parts = strsplit(pathname, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep filename];
myhandles = guihandles(gcf);
set(myhandles.text_area2,'string',disp_str);
set(myhandles.edit_name_area2,'string',name);
set(myhandles.text_area2,'TooltipString',f);

% --- Executes on button press in radiobutton_against_GM.
function radiobutton_against_GM_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_against_GM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_against_GM
   
    set(handles.radiobutton_against_Cortex,'Value',0);
    ontology = readtable('tmp_ontology.csv');
    ind_structure=find(ontology.id==4006);
    set(handles.popupmenu_select_structure,'Value',ind_structure+1);

% --- Executes on button press in radiobutton_against_Cortex.
function radiobutton_against_Cortex_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_against_Cortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_against_Cortex
    set(handles.radiobutton_against_GM,'Value',0);
    ontology = readtable('tmp_ontology.csv');
    ind_structure=find(ontology.id==4008);
    set(handles.popupmenu_select_structure,'Value',ind_structure+1);


function edit_name_area1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_name_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_name_area1 as text
%        str2double(get(hObject,'String')) returns contents of edit_name_area1 as a double


% --- Executes during object creation, after setting all properties.
function edit_name_area1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_name_area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_name_area2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_name_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_name_area2 as text
%        str2double(get(hObject,'String')) returns contents of edit_name_area2 as a double


% --- Executes during object creation, after setting all properties.
function edit_name_area2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_name_area2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_reload_latest_settings.
function pushbutton_reload_latest_settings_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reload_latest_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
load('latest_settings.mat')
myhandles = guihandles(gcf);
[pathstr,name,ext] = fileparts(output_parameter.pmap{1});
parts = strsplit(pathstr, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep name ext];
set(myhandles.text_area1,'string',disp_str);
set(myhandles.edit_name_area1,'string',name);
set(myhandles.text_area1,'TooltipString',output_parameter.pmap{1});

[pathstr,name,ext] = fileparts(output_parameter.pmap{2});
parts = strsplit(pathstr, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep  name ext];
set(myhandles.text_area2,'string',disp_str);
set(myhandles.edit_name_area2,'string',name);
set(myhandles.text_area2,'TooltipString',output_parameter.pmap{2});

[pathstr,name,ext] = fileparts(output_parameter.gene_list{1});
parts = strsplit(pathstr, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep  name ext];
set(myhandles.text_gene_list,'string',disp_str);
set(myhandles.text_gene_list,'TooltipString',output_parameter.gene_list{1});

[pathstr,name,ext] = fileparts(output_parameter.output_folder{1});
parts = strsplit(pathstr, filesep);
DirPart = parts{end};
disp_str=['.' filesep DirPart filesep  name ext];
set(myhandles.text_output_folder,'string',disp_str);
set(myhandles.text_output_folder,'TooltipString',output_parameter.output_folder{1});

% [pathstr,name,ext] = fileparts(output_parameter.project_name{1});
% parts = strsplit(pathstr, filesep);
% DirPart = parts{end};
% disp_str=['.' filesep DirPart filesep  name ext];
set(myhandles.edit_project_name,'string',output_parameter.project_name{1});
set(myhandles.edit_project_name,'TooltipString',output_parameter.project_name{1});

set(myhandles.edit_map_threshold,'string',output_parameter.map_threshold{1});

available_donors=get(myhandles.listbox_donors,'string');
   
select_donors=output_parameter.donors;
mem=ismember(available_donors,select_donors);
vals=find(mem);
set(myhandles.listbox_donors,'Value',vals);

if output_parameter.search_mode{1}==2
    set(myhandles.checkbox_vs_structure,'Value',1);
    set(handles.text_area2,'Enable','off');     
    set(handles.edit_name_area2,'Enable','off');
    set(handles.pushbutton_select_area2,'Enable','off');   
    set(handles.popupmenu_select_structure,'Enable','on');
    set(handles.radiobutton_against_GM,'Enable','on');
    set(handles.radiobutton_against_Cortex,'Enable','on');
    set(handles.radiobutton_against_GM,'Value',0);
    set(handles.radiobutton_against_Cortex,'Value',0);
    ontology = readtable('tmp_ontology.csv');
    ind_structure=find(ontology.id==output_parameter.struct_id{1});
    set(handles.popupmenu_select_structure,'Value',ind_structure+1);
end

catch
    h = msgbox('Sorry! Loading latest setting was not possible.');
    
end



function edit_map_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_map_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_map_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_map_threshold as a double


% --- Executes during object creation, after setting all properties.
function edit_map_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_map_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vs_structure.
function checkbox_vs_structure_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vs_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vs_structure
value_checkbox_vs_structure=get(hObject,'Value');
if value_checkbox_vs_structure==1   % checkbox wurde angewählt, name2 ausblenden
    set(handles.text_area2,'Enable','off');     
    set(handles.edit_name_area2,'Enable','off');
    set(handles.pushbutton_select_area2,'Enable','off');   
    set(handles.popupmenu_select_structure,'Enable','on');
    %set(handles.radiobutton_against_GM,'Enable','on');
    %set(handles.radiobutton_against_Cortex,'Enable','on');
elseif value_checkbox_vs_structure==0
    set(handles.text_area2,'Enable','on'); 
    set(handles.edit_name_area2,'Enable','on');
    set(handles.pushbutton_select_area2,'Enable','on');    
    set(handles.popupmenu_select_structure,'Enable','off');
    %set(handles.radiobutton_against_GM,'Enable','off');
    %set(handles.radiobutton_against_Cortex,'Enable','off');
    
end

% --- Executes on selection change in popupmenu_select_structure.
function popupmenu_select_structure_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_structure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_structure
contents = cellstr(get(hObject,'String'));
contents{get(hObject,'Value')} % structure name
ontology = readtable('tmp_ontology.csv');
ontology.id(get(hObject,'Value')-1) %structure id  (minus 1 wegen "please select structure" am anfang
%     set(handles.radiobutton_against_Cortex,'Value',0);
%     set(handles.radiobutton_against_GM,'Value',0);

% --- Executes during object creation, after setting all properties.
function popupmenu_select_structure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% url='http://api.brain-map.org/api/v2/data/query.csv?criteria=model::Structure,rma::criteria,[ontology_id$eq7],rma::options[order$eq%27structures.graph_order%27][num_rows$eqall]';
% urlwrite(url,'tmp_ontology.csv');
% ontology = readtable('tmp_ontology.csv');
% first={'please select structure'};
% names=cell(ontology.name);
% combined=[first;names];
% set(hObject,'String',combined);
% %set(hObject,'Value',cell(ontology.id));
% % test=ontology.id;
% % set(hObject,'Value',test)
