function varargout = NOW_GUI(varargin)
% NOW_GUI MATLAB code for NOW_GUI.fig
%      NOW_GUI, by itself, creates a new NOW_GUI or raises the existing
%      singleton*.
%
%      H = NOW_GUI returns the handle to a new NOW_GUI or the handle to
%      the existing singleton*.
%
%      NOW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NOW_GUI.M with the given input arguments.
%
%      NOW_GUI('Property','Value',...) creates a new NOW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NOW_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NOW_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NOW_GUI

% Last Modified by GUIDE v2.5 20-Sep-2020 20:11:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NOW_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NOW_GUI_OutputFcn, ...
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


% --- Executes just before NOW_GUI is made visible.
function NOW_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NOW_GUI (see VARARGIN)


[absolutePathToNOW, ~, ~] = fileparts(which('NOW_GUI'));
addpath(absolutePathToNOW)

% set(handles.mainWindow, 'units', 'normalized', 'position', [0.05 0.15 0.5 0.6])

% Update handles structure
settings.name = 'NOW1';
problem = optimizationProblem(settings);

handles.problem = problem;
handles.queue = [];
% Choose default command line output for NOW_GUI
handles.output = [];

set(handles.discretizationStepsTextBox, 'String', num2str(problem.N))
set(handles.targetTensorTable, 'Data', problem.targetTensor)
set(handles.encodingTensorTable, 'Data', problem.targetTensor)
set(handles.slewRateTextBox, 'String', num2str(problem.sMax))
set(handles.gMaxTextBox, 'String', num2str(problem.gMax))
set(handles.maxNormRadioButton, 'Value', problem.useMaxNorm);
set(handles.EuclideanNormRadioButton, 'Value', ~problem.useMaxNorm);

% Timings panel
totalTimeRequested = problem.durationFirstPartRequested + problem.durationZeroGradientRequested + problem.durationSecondPartRequested;
set(handles.totalTimeRequestedText, 'String', num2str(totalTimeRequested))
set(handles.durationFirstPartTextBox, 'String', num2str(problem.durationFirstPartRequested))
set(handles.durationZeroGradientTextBox, 'String', num2str(problem.durationZeroGradientRequested))
set(handles.durationSecondPartTextBox, 'String', num2str(problem.durationSecondPartRequested))
set(handles.durationFirstPartActualText, 'String', sprintf('%0.2f',problem.durationFirstPartActual))
set(handles.durationZeroGradientActualText, 'String', sprintf('%0.2f',problem.durationZeroGradientActual))
set(handles.durationSecondPartActualText, 'String', sprintf('%0.2f',problem.durationSecondPartActual))
set(handles.totalTimeActualText, 'String', sprintf('%0.2f',problem.totalTimeActual))

set(handles.enforceSymmetryCheckBox,'Value', problem.enforceSymmetry)
%set(handles.redoIfFailedCheckBox,'Value', problem.redoIfFailed)
set(handles.doMaxwellCheckBox,'Value', problem.doMaxwellComp)
set(handles.heatDissipationTextBox, 'String', num2str(problem.eta))
set(handles.nameTextBox, 'String', problem.name)

% Set gradient plot to be default
set(handles.plotButtonGroup,  'SelectedObject', handles.gRadioButton)

guidata(hObject, handles)

% UIWAIT makes NOW_GUI wait for user response (see UIRESUME)
% uiwait(handles.mainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = NOW_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function updateTimings(hObject, handles)
handles.problem.durationFirstPartRequested = str2double(get(handles.durationFirstPartTextBox, 'String'));
handles.problem.durationZeroGradientRequested = str2double(get(handles.durationZeroGradientTextBox, 'String'));
handles.problem.durationSecondPartRequested = str2double(get(handles.durationSecondPartTextBox, 'String'));
totalTime = handles.problem.durationFirstPartRequested + handles.problem.durationZeroGradientRequested + handles.problem.durationSecondPartRequested;

set(handles.totalTimeRequestedText, 'String', num2str(totalTime))

handles.problem = optimizationProblem(handles.problem);

set(handles.durationFirstPartActualText, 'String', sprintf('%0.2f',handles.problem.durationFirstPartActual))
set(handles.durationZeroGradientActualText, 'String', sprintf('%0.2f',handles.problem.durationZeroGradientActual))
set(handles.durationSecondPartActualText, 'String', sprintf('%0.2f',handles.problem.durationSecondPartActual))
set(handles.totalTimeActualText, 'String', sprintf('%0.2f',handles.problem.totalTimeActual))

% if handles.problem.durationFirstPartActual ~= handles.problem.durationSecondPartActual
%     handles.problem.enforceSymmetry = false;
%     handles.problem = optimizationProblem(handles.problem);
%     set(handles.enforceSymmetryCheckBox, 'Value', false)
%     set(handles.enforceSymmetryCheckBox, 'Enable', 'off')
% else
%     set(handles.enforceSymmetryCheckBox, 'Enable', 'on')
% end

guidata(hObject, handles)

function gMaxTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to gMaxTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gMaxTextBox as text
%        str2double(get(hObject,'String')) returns contents of gMaxTextBox as a double
gMax = str2double(get(hObject,'String'));
handles.problem.gMax = gMax;
handles.problem = optimizationProblem(handles.problem);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function gMaxTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gMaxTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slewRateTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to slewRateTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slewRateTextBox as text
%        str2double(get(hObject,'String')) returns contents of slewRateTextBox as a double
slewRate = str2double(get(hObject,'String'));
handles.problem.sMax = slewRate;
handles.problem = optimizationProblem(handles.problem);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slewRateTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slewRateTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function durationFirstPartTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to durationFirstPartTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of durationFirstPartTextBox as text
%        str2double(get(hObject,'String')) returns contents of durationFirstPartTextBox as a double
updateTimings(hObject, handles)


% --- Executes during object creation, after setting all properties.
function durationFirstPartTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to durationFirstPartTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function settingsPanel_ButtonDownFcn(hObject,eventdata, handles)
% Intentionally empty

function durationZeroGradientTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to durationZeroGradientTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of durationZeroGradientTextBox as text
%        str2double(get(hObject,'String')) returns contents of durationZeroGradientTextBox as a double
updateTimings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function durationZeroGradientTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to durationZeroGradientTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function durationSecondPartTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to durationSecondPartTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of durationSecondPartTextBox as text
%        str2double(get(hObject,'String')) returns contents of durationSecondPartTextBox as a double
updateTimings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function durationSecondPartTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to durationSecondPartTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function totalTimeRequestedText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalTimeRequestedText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function heatDissipationTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to heatDissipationTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heatDissipationTextBox as text
%        str2double(get(hObject,'String')) returns contents of heatDissipationTextBox as a double
eta = str2double(get(hObject,'String'));
if eta > 1
    eta = 1;
    set(hObject,'String',num2str(eta))
elseif eta < 0
    eta = 0;
    set(hObject,'String',num2str(eta))
end
handles.problem.eta = eta;
handles.problem = optimizationProblem(handles.problem);
guidata(hObject, handles)
    

% --- Executes during object creation, after setting all properties.
function heatDissipationTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heatDissipationTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in enforceSymmetryCheckBox.
function enforceSymmetryCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to enforceSymmetryCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of enforceSymmetryCheckBox
enforceSymmetry = get(hObject,'Value');
handles.problem.enforceSymmetry = enforceSymmetry;
handles.problem = optimizationProblem(handles.problem);
guidata(hObject, handles)
updateTimings(hObject, handles)


% --- Executes on button press in maxNormRadioButton.
function maxNormRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to maxNormRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxNormRadioButton
useMaxNorm = get(hObject,'Value');
handles.problem.useMaxNorm = useMaxNorm;
handles.problem = optimizationProblem(handles.problem);
set(handles.maxNormRadioButton, 'Value', useMaxNorm);
set(handles.EuclideanNormRadioButton, 'Value', ~useMaxNorm);

handles.problem.N = str2double(get(handles.discretizationStepsTextBox,'String'));
handles.problem = optimizationProblem(handles.problem);
guidata(hObject, handles);
updateTimings(hObject, handles);


% --- Executes on button press in runOptimizationPushButton.
function runOptimizationPushButton_Callback(hObject, eventdata, handles)
% hObject    handle to runOptimizationPushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[result, problem] = optimize(handles.problem);
optimization.result = result;
optimization.problem = problem;
handles.output = [handles.output, optimization];
guidata(hObject, handles);
addNewResult = true;
updateResults(hObject, handles, addNewResult);
updateName(hObject, handles);

function updateName(hObject,handles)
newDefaultName = ['NOW' num2str(length(handles.output)+length(handles.queue)+1)];
set(handles.nameTextBox, 'String', newDefaultName)
handles.problem.name = newDefaultName;
guidata(hObject, handles);

% --- Executes on button press in EuclideanNormRadioButton.
function EuclideanNormRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to EuclideanNormRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EuclideanNormRadioButton
useMaxNorm = ~get(hObject,'Value');
handles.problem.useMaxNorm = useMaxNorm;
set(handles.maxNormRadioButton, 'Value', useMaxNorm);
set(handles.EuclideanNormRadioButton, 'Value', ~useMaxNorm);

handles.problem.N = str2double(get(handles.discretizationStepsTextBox,'String'));
guidata(hObject, handles);
updateTimings(hObject, handles);



% --- Executes when entered data in editable cell(s) in targetTensorTable.
function targetTensorTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to targetTensorTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.problem.targetTensor(eventdata.Indices(1),eventdata.Indices(2)) = str2double(eventdata.EditData);
guidata(hObject, handles)


function plotResult(hObject, handles)
obj = get(handles.plotButtonGroup,'SelectedObject');
switch get(obj,'Tag') % Get Tag of selected object.
    case 'qRadioButton'
        plotq(hObject, handles);
    case 'gRadioButton'
        plotGradient(hObject, handles);
    case 'slewRadioButton'
        plotSlew(hObject, handles);
    case 'maxwellRadioButton'
        plotMaxwell(hObject, handles);
    otherwise
        error('Plot command not recognized!')
end


function plotGradient(hObject, handles)
if isfield(handles.output, 'result')
    index = get(handles.outputNameDropDown,'Value');
    dt = handles.output(index).problem.dt;
    T = handles.output(index).problem.totalTimeActual;
    gMax = handles.output(index).problem.gMax;
    plot(handles.plotAxes,0:dt:T,handles.output(index).result.g)
    axis([0 T -1.1*gMax 1.1*gMax])
    xlabel(handles.plotAxes,'Time [ms]')
    ylabel(handles.plotAxes,'Gradient amplitude [mT/m]')
end

function plotq(hObject, handles)
if isfield(handles.output, 'result')
    index = get(handles.outputNameDropDown,'Value');
    dt = handles.output(index).problem.dt;
    T = handles.output(index).problem.totalTimeActual;
    plot(handles.plotAxes,dt/2:dt:(T-dt/2),handles.output(index).result.q*1e-6)
    xlabel(handles.plotAxes,'Time [ms]')
    ylabel(handles.plotAxes,'q [(\mu m)^{-1}]')
end

function plotMaxwell(hObject, handles)
if isfield(handles.output, 'result')
    index = get(handles.outputNameDropDown,'Value');
    [k, m, g2t] = now_maxwell_coeff(handles.output(index));
    dt = handles.output(index).problem.dt;
    T = handles.output(index).problem.totalTimeActual;
    plot(dt/2:dt:(T+dt/2), g2t')
    legend('xx', 'yy', 'zz', 'xy', 'xz', 'yz')
    text(.1, .1, 'k = ', 'units', 'normalized', 'hor', 'right')
    text(.1, .1, num2str(k, ' %0.1e'), 'units', 'normalized', 'hor', 'left')
    text(.1, .9, 'm = ', 'units', 'normalized', 'hor', 'right')
    text(.1, .9, [num2str(m*1e9, ' %0.1e') '  [mT/m s]'], 'units', 'normalized', 'hor', 'left')
    title('G(t)^T\cdotG(t)')
    xlabel(handles.plotAxes,'Time [ms]')
    ylabel(handles.plotAxes,'Maxwell')
end

function plotSlew(hObject, handles)
if  isfield(handles.output, 'result')
    index = get(handles.outputNameDropDown,'Value');
    dt = handles.output(index).problem.dt;
    T = handles.output(index).problem.totalTimeActual;
    sMax = handles.output(index).problem.sMax;
    plot(handles.plotAxes,-dt/2:dt:(T+dt/2),handles.output(index).result.slew)
    axis([-dt/2 (T+dt/2) -1.1*sMax 1.1*sMax])
    xlabel(handles.plotAxes,'Time [ms]')
    ylabel(handles.plotAxes,'Slew rate [T/m/s]')
end

% --- Executes when selected object is changed in plotButtonGroup.
function plotButtonGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in plotButtonGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
%     case 'qRadioButton'
%         plotq(hObject, handles);
%     case 'gRadioButton'
%         plotGradient(hObject, handles);
%     case 'slewRadioButton'
%         plotSlew(hObject, handles);
% end


% --- Executes on button press in queuePushButton.
function queuePushButton_Callback(hObject, eventdata, handles)
% hObject    handle to queuePushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'queue')
    handles.queue = handles.problem;
else
    handles.queue = [handles.queue, handles.problem];
end
guidata(hObject, handles);
updateName(hObject,handles)


% --- Executes on button press in runQueuePushButton.
function runQueuePushButton_Callback(hObject, eventdata, handles)
% hObject    handle to runQueuePushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = arrayfun(@(i) optimize(handles.queue(i)), 1:length(handles.queue), 'UniformOutput',false);
optimization = struct('result', tmp, 'problem', num2cell(handles.queue));
handles.output = [handles.output, optimization];
handles.queue = [];
guidata(hObject, handles);
addNewResult = true;
updateResults(hObject, handles, addNewResult)
set(handles.outputNameDropDown,'Value', length(handles.output));
updateName(hObject, handles)


% --- Executes on button press in viewQueuePushButton.
function viewQueuePushButton_Callback(hObject, eventdata, handles)
% hObject    handle to viewQueuePushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savePushButton.
function savePushButton_Callback(hObject, eventdata, handles)
% hObject    handle to savePushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

index = get(handles.outputNameDropDown,'Value');
names = get(handles.outputNameDropDown,'String');
selectedName = names{index};
[file,path] = uiputfile([selectedName '.txt'],'Save file name');
createSaveFile(handles.output(index).result, file, path)

function createSaveFile(result, out_name, out_dir)

if nargin < 2
    out_name = 'NOW';
end

if nargin < 3
    out_dir = pwd;
end

formatspec = '%8.5f %8.5f %8.5f\r\n';
g_norm = max(abs(result.g(:)));

[~, out_name, ext] = fileparts(out_name);

out_mat_whole = result.g' / g_norm;

outfile_whole = fopen([out_dir filesep out_name ext], 'w');
fprintf(outfile_whole, '%8.0i\r\n', size(out_mat_whole, 2));
fprintf(outfile_whole, formatspec, out_mat_whole);
fclose (outfile_whole);

z_ind = result.zind;

% If the waveform is split into two parts, it is also stored as two
% individual files with endings "_A" and "_B".
if ~isempty(z_ind)
    end_a     = z_ind(1) + 1;
    beg_b     = z_ind(end) + 1;
    out_mat_a = result.g(1:end_a  , :)' / g_norm;
    out_mat_b = result.g(beg_b:end, :)' / g_norm;
    
    outfile_a = fopen([out_dir filesep out_name '_A' ext], 'w');
    fprintf(outfile_a, '%8.0i\r\n', size(out_mat_a, 2));
    fprintf(outfile_a, formatspec, out_mat_a);
    fclose (outfile_a);
    
    % The sign of the second part is inverted, so that the waveform
    % conforms to the gradient system. For example, a Stejskal-Tanner
    % waveform should be [0 1 1 0 ... 0 1 1 0] along one axis.
    outfile_b = fopen([out_dir filesep out_name '_B' ext], 'w');
    fprintf(outfile_b, '%8.0i\r\n', size(out_mat_b, 2));
    fprintf(outfile_b, formatspec, -out_mat_b); 
    fclose (outfile_b);
end



% --- Executes on selection change in outputNameDropDown.
function outputNameDropDown_Callback(hObject, eventdata, handles)
% hObject    handle to outputNameDropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns outputNameDropDown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from outputNameDropDown
addNewResult = false;
updateResults(hObject, handles, addNewResult);

% --- Executes during object creation, after setting all properties.
function outputNameDropDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputNameDropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateResults(hObject, handles, addNewResult)
set(handles.savePushButton, 'Enable', 'On');

names = cell(1, length(handles.output));
for i = 1:length(names)
    names{i} = handles.output(i).problem.name;
end
set(handles.outputNameDropDown, 'String', names);
if addNewResult
    set(handles.outputNameDropDown,'Value',length(handles.output))
    assignin('base', 'optimization', handles.output)
end
index = get(handles.outputNameDropDown,'Value');
handles.problem = handles.output(index).problem;
guidata(hObject, handles);
plotResult(hObject, handles);
setAllFieldsAsProblem(hObject, handles, index)

function setAllFieldsAsProblem(hObject, handles, index)
set(handles.bResultText, 'String', [sprintf('%0.2f',handles.output(index).result.b) ' s/mm^2'])
set(handles.encodingTensorTable, 'Data', handles.output(index).result.B/trace(handles.output(index).result.B)*trace(handles.output(index).problem.targetTensor))

set(handles.targetTensorTable, 'Data', handles.output(index).problem.targetTensor)
set(handles.slewRateTextBox, 'String', num2str(handles.output(index).problem.sMax))
set(handles.gMaxTextBox, 'String', num2str(handles.output(index).problem.gMax))
set(handles.maxNormRadioButton, 'Value', handles.output(index).problem.useMaxNorm);
set(handles.EuclideanNormRadioButton, 'Value', ~handles.output(index).problem.useMaxNorm);

order = max(handles.output(index).problem.motionCompensation.order);
if isempty(order)
    order = 0;
end
  
switch order
    case 0
        set(handles.motionCompensationOrder0, 'Value', true)
        set(handles.motionCompensationOrder1, 'Value', false)
        set(handles.motionCompensationOrder2, 'Value', false)  
    case 1
        set(handles.motionCompensationOrder0, 'Value', false)
        set(handles.motionCompensationOrder1, 'Value', true)
        set(handles.motionCompensationOrder2, 'Value', false)  
    case 2
        set(handles.motionCompensationOrder0, 'Value', false)
        set(handles.motionCompensationOrder1, 'Value', false)
        set(handles.motionCompensationOrder2, 'Value', true)  
end

% Timings panel
totalTimeRequested = handles.output(index).problem.durationFirstPartRequested + handles.output(index).problem.durationZeroGradientRequested + handles.output(index).problem.durationSecondPartRequested;
set(handles.totalTimeRequestedText, 'String', num2str(totalTimeRequested))
set(handles.durationFirstPartTextBox, 'String', num2str(handles.output(index).problem.durationFirstPartRequested))
set(handles.durationZeroGradientTextBox, 'String', num2str(handles.output(index).problem.durationZeroGradientRequested))
set(handles.durationSecondPartTextBox, 'String', num2str(handles.output(index).problem.durationSecondPartRequested))
set(handles.durationFirstPartActualText, 'String', sprintf('%0.2f',handles.output(index).problem.durationFirstPartActual))
set(handles.durationZeroGradientActualText, 'String', sprintf('%0.2f',handles.output(index).problem.durationZeroGradientActual))
set(handles.durationSecondPartActualText, 'String', sprintf('%0.2f',handles.output(index).problem.durationSecondPartActual))
set(handles.totalTimeActualText, 'String', sprintf('%0.2f',handles.output(index).problem.totalTimeActual))

set(handles.enforceSymmetryCheckBox,'Value', handles.output(index).problem.enforceSymmetry)
%set(handles.redoIfFailedCheckBox,'Value', handles.output(index).problem.redoIfFailed)
set(handles.doMaxwellCheckBox,'Value', handles.output(index).problem.doMaxwellComp)
set(handles.heatDissipationTextBox, 'String', num2str(handles.output(index).problem.eta))
set(handles.nameTextBox, 'String', handles.output(index).problem.name)


function nameTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to nameTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nameTextBox as text
%        str2double(get(hObject,'String')) returns contents of nameTextBox as a double
handles.problem.name = get(hObject,'String');
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function nameTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nameTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in redoIfFailedCheckBox.
function redoIfFailedCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to redoIfFailedCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of redoIfFailedCheckBox
redoIfFailed = get(hObject,'Value');
handles.problem.redoIfFailed = redoIfFailed;
guidata(hObject, handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over qRadioButton.
function qRadioButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to qRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% plotq(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over gRadioButton.
function gRadioButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to gRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% plotGradient(hObject, handles);


% --- Executes on button press in slewRadioButton.
function slewRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to slewRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slewRadioButton
plotSlew(hObject, handles);


% --- Executes on button press in gRadioButton.
function gRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to gRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gRadioButton
plotGradient(hObject, handles);


% --- Executes on button press in qRadioButton.
function qRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to qRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of qRadioButton
plotq(hObject, handles)


% --- Executes on button press in maxwellRadioButton.
function Maxwellradiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to maxwellRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxwellRadioButton
plotMaxwell(hObject, handles)


% --- Executes on button press in doMaxwellCheckBox.
function doMaxwellCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to doMaxwellCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doMaxwellCheckBox
doMaxwellComp = get(hObject,'Value');
handles.problem.doMaxwellComp = doMaxwellComp;
guidata(hObject, handles)



function discretizationStepsTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to discretizationStepsTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of discretizationStepsTextBox as text
%        str2double(get(hObject,'String')) returns contents of discretizationStepsTextBox as a double
handles.problem.N = str2double(get(hObject,'String'));
guidata(hObject, handles);
updateTimings(hObject, handles);


% --- Executes during object creation, after setting all properties.
function discretizationStepsTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to discretizationStepsTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in motionCompensationOrder0.
function motionCompensationOrder0_Callback(hObject, eventdata, handles)
% hObject    handle to motionCompensationOrder0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of motionCompensationOrder0
order = []; % order zero is always ensured by the spin echo condition
handles.problem.motionCompensation.order = order;
handles.problem = optimizationProblem(handles.problem);
set(handles.motionCompensationOrder0, 'Value', true);
set(handles.motionCompensationOrder1, 'Value', false);
set(handles.motionCompensationOrder2, 'Value', false);

guidata(hObject, handles);


% --- Executes on button press in motionCompensationOrder1.
function motionCompensationOrder1_Callback(hObject, eventdata, handles)
% hObject    handle to motionCompensationOrder1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of motionCompensationOrder1
order = 1; % order zero is always ensured by the spin echo condition
handles.problem.motionCompensation.order = order;
handles.problem.motionCompensation.linear = true;
handles.problem = optimizationProblem(handles.problem);
set(handles.motionCompensationOrder0, 'Value', false);
set(handles.motionCompensationOrder1, 'Value', true);
set(handles.motionCompensationOrder2, 'Value', false);

guidata(hObject, handles);

% --- Executes on button press in motionCompensationOrder2.
function motionCompensationOrder2_Callback(hObject, eventdata, handles)
% hObject    handle to motionCompensationOrder2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of motionCompensationOrder2
order = [1, 2]; % order zero is always ensured by the spin echo condition
handles.problem.motionCompensation.order = order;
handles.problem.motionCompensation.linear = [true, true];
handles.problem = optimizationProblem(handles.problem);
set(handles.motionCompensationOrder0, 'Value', false);
set(handles.motionCompensationOrder1, 'Value', false);
set(handles.motionCompensationOrder2, 'Value', true);

guidata(hObject, handles);
