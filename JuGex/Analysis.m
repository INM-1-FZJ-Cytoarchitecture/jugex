function varargout = Analysis(varargin)
% ANALYSIS MATLAB code for Analysis.fig
%      ANALYSIS, by itself, creates a new ANALYSIS or raises the existing
%      singleton*.
%
%      H = ANALYSIS returns the handle to a new ANALYSIS or the handle to
%      the existing singleton*.
%
%      ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSIS.M with the given input arguments.
%
%      ANALYSIS('Property','Value',...) creates a new ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Analysis

% Last Modified by GUIDE v2.5 29-Sep-2016 14:43:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Analysis_OpeningFcn, ...
    'gui_OutputFcn',  @Analysis_OutputFcn, ...
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


% --- Executes just before Analysis is made visible.
function Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Analysis (see VARARGIN)

% Choose default command line output for Analysis
handles.output = hObject;

set(gcf,'name','Analysis','numbertitle','off');

% Update handles structure
guidata(hObject, handles);

if ~isempty(varargin)
    if length(varargin{1})>0
        if strcmp(varargin,'gui')
            load('latest_settings.mat');
            myhandles = guihandles(gcf);
            %         [pathstr,name,ext] = fileparts(output_parameter.gene_list{1});
            %         parts = strsplit(pathstr, filesep);
            %         DirPart = parts{end};
            %         disp_str=['.' filesep DirPart filesep  name ext];
            %         set(myhandles.text_gene_list,'string',disp_str);
            %         set(myhandles.text_gene_list,'TooltipString',output_parameter.gene_list{1});
            
            [pathstr,name,ext] = fileparts(output_parameter.output_folder{1});
            parts = strsplit(pathstr, filesep);
            DirPart = parts{end};
            disp_str=['.' filesep name filesep 'extracted_data' filesep output_parameter.project_name{1} '_th_' output_parameter.map_threshold{1} '.mat'];
            guidata(hObject, handles);
            set(myhandles.text_extracted_data_file,'string',disp_str);
            set(myhandles.text_extracted_data_file,'TooltipString',[output_parameter.output_folder{1} filesep 'extracted_data' filesep output_parameter.project_name{1} '_th_' output_parameter.map_threshold{1} '.mat']);
            if output_parameter.search_mode{1}==1
                set(myhandles.edit_name_area1,'string',output_parameter.name{1});
                set(myhandles.edit_name_area2,'string',output_parameter.name{2});
            elseif output_parameter.search_mode{1}==2
                set(myhandles.edit_name_area1,'string',output_parameter.name{1});
                ontology = readtable('tmp_ontology.csv');
                ind_structure=find(ontology.id==output_parameter.struct_id{1});
                set(myhandles.edit_name_area2,'string',ontology.name{ind_structure});
            end
            
        end
    end
end
% UIWAIT makes Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Analysis_OutputFcn(hObject, eventdata, handles)
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

myhandles = guihandles(gcf);
% output_parameter.name={get(myhandles.edit_name_area1,'String'),get(myhandles.edit_name_area2,'String')};
% output_parameter.gene_list={get(myhandles.text_gene_list,'TooltipString')};
% f = fullfile(pathname,filename);
% myhandles.path_parameter_file=f;
[path file ext]=fileparts(handles.text_extracted_data_file.TooltipString);
para_path=[path filesep file '_parameter' ext];
load(para_path);
output_parameter=parameter;
clear parameter
output_parameter.extracted_data_file={get(myhandles.text_extracted_data_file,'TooltipString')};


if get(myhandles.radiobutton_analyze_gene_level,'Value')==1 && get(myhandles.radiobutton_analyze_probeid_level,'Value')==0
    output_parameter.analysis_mode=1;
elseif get(myhandles.radiobutton_analyze_gene_level,'Value')==0 && get(myhandles.radiobutton_analyze_probeid_level,'Value')==1
    output_parameter.analysis_mode=2;
else
    error('No analysis mode selected');
end
output_parameter.nrepetitions=str2num(get(myhandles.edit_nrepetitions,'String'));


anovan_zscores_factor1_specimen_factor2_area('gui',output_parameter);


function edit_nrepetitions_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nrepetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nrepetitions as text
%        str2double(get(hObject,'String')) returns contents of edit_nrepetitions as a double


% --- Executes during object creation, after setting all properties.
function edit_nrepetitions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nrepetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_extracted_data_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_extracted_data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_gene_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_change_input_file.
function pushbutton_change_input_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_change_input_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myhandles = guihandles(gcf);
[filename, pathname] = uigetfile({'*_parameter.mat'},'Select file with extracted data',['.' filesep 'output' filesep 'extracted_data' filesep]);
f = fullfile(pathname,filename);
myhandles.path_parameter_file=f;
guidata(gcf, myhandles);
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

%
%         [pathstr,name,ext] = fileparts(output_parameter.gene_list{1});
%         parts = strsplit(pathstr, filesep);
%         DirPart = parts{end};
%         disp_str=['.' filesep DirPart filesep  name ext];
%         set(myhandles.text_gene_list,'string',disp_str);
%         set(myhandles.text_gene_list,'TooltipString',output_parameter.gene_list{1});
%
%         [pathstr,name,ext] = fileparts(output_parameter.output_folder{1});
%         parts = strsplit(pathstr, filesep);
%         DirPart = parts{end};
%         disp_str=['.' filesep name filesep 'extracted_data' filesep output_parameter.project_name{1} '_th_' output_parameter.map_threshold{1} '.mat'];
%         set(myhandles.text_extracted_data_file,'string',disp_str);
%         set(myhandles.text_extracted_data_file,'TooltipString',[output_parameter.output_folder{1} filesep 'extracted_data' filesep output_parameter.project_name{1} '_th_' output_parameter.map_threshold{1} '.mat']);
if output_parameter.search_mode{1}==1
    set(myhandles.edit_name_area1,'string',output_parameter.name{1});
    set(myhandles.edit_name_area2,'string',output_parameter.name{2});
elseif output_parameter.search_mode{1}==2
    set(myhandles.edit_name_area1,'string',output_parameter.name{1});
    ontology = readtable('tmp_ontology.csv');
    ind_structure=find(ontology.id==output_parameter.struct_id{1});
    set(myhandles.edit_name_area2,'string',ontology.name{ind_structure});
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
disp_str=['.' DirPart filesep filename];
myhandles = guihandles(gcf);
set(myhandles.text_gene_list,'string',disp_str);
set(myhandles.text_gene_list,'TooltipString',f);


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
