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

% Last Modified by GUIDE v2.5 21-Sep-2016 13:00:51

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



% set(handles.mainWindow, 'units', 'normalized', 'position', [0.05 0.15 0.5 0.6])

% Update handles structure
settings.name = 'NOW1';
problem = optimizationProblem(settings);

handles.problem = problem;
handles.queue = [];
% Choose default command line output for NOW_GUI
handles.output = [];

set(handles.targetTensorTable, 'Data', problem.targetTensor)
set(handles.encodingTensorTable, 'Data', problem.targetTensor)
set(handles.slewRateTextBox, 'String', num2str(problem.sMax))
set(handles.gMaxTextBox, 'String', num2str(problem.gMax))
set(handles.maxNormRadioButton, 'Value', problem.useMaxNorm);
set(handles.EuclideanNormRadioButton, 'Value', ~problem.useMaxNorm);

precomputedDiscretizations = getPrecomputedDiscretizations(hObject, handles);
set(handles.discretizationStepsDropDown, 'String', precomputedDiscretizations)
set(handles.discretizationStepsDropDown, 'Value', getIndexOfSameN(handles.problem, precomputedDiscretizations))

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
set(handles.redoIfFailedCheckBox,'Value', problem.redoIfFailed)
set(handles.heatDissipationTextBox, 'String', num2str(problem.eta))
set(handles.nameTextBox, 'String', problem.name)


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

function precomputedDiscretizations = getPrecomputedDiscretizations(hObject, handles)
    if handles.problem.useMaxNorm
        fileEnding = 'MaxNorm.m';      
    else
        fileEnding = '2Norm.m';
    end
    D = dir(['private/*' fileEnding]);
    precomputedDiscretizations = cell(length(D),1);
    for i = 1:length(D)
        n1 = length('nonlcon');
        n2 = length('points');
        n3 = length(fileEnding);
        
        precomputedDiscretizations{i} = D(i).name((n1+1):(end-(n2+n3)));
    end
    precomputedDiscretizations = cellfun(@(s) num2str(s), num2cell(sort(str2double(precomputedDiscretizations))),'UniformOutput',false);


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

if handles.problem.durationFirstPartActual ~= handles.problem.durationSecondPartActual
    handles.problem.enforceSymmetry = false;
    handles.problem = optimizationProblem(handles.problem);
    set(handles.enforceSymmetryCheckBox, 'Value', false)
    set(handles.enforceSymmetryCheckBox, 'Enable', 'off')
else
    set(handles.enforceSymmetryCheckBox, 'Enable', 'on')
end

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
precomputedDiscretizations = getPrecomputedDiscretizations(hObject, handles);
set(handles.discretizationStepsDropDown, 'String', precomputedDiscretizations)
set(handles.discretizationStepsDropDown, 'Value', min(get(handles.discretizationStepsDropDown,'Value'), length(precomputedDiscretizations)));
contents = cellstr(get(handles.discretizationStepsDropDown,'String'));
handles.problem.N = str2double(contents{get(handles.discretizationStepsDropDown,'Value')});
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
precomputedDiscretizations = getPrecomputedDiscretizations(hObject, handles);
set(handles.discretizationStepsDropDown, 'String', precomputedDiscretizations)
set(handles.discretizationStepsDropDown, 'Value', min(get(handles.discretizationStepsDropDown,'Value'), length(precomputedDiscretizations)));
contents = cellstr(get(handles.discretizationStepsDropDown,'String'));
handles.problem.N = str2double(contents{get(handles.discretizationStepsDropDown,'Value')});
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


% --- Executes on selection change in discretizationStepsDropDown.
function discretizationStepsDropDown_Callback(hObject, eventdata, handles)
% hObject    handle to discretizationStepsDropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns discretizationStepsDropDown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from discretizationStepsDropDown
contents = cellstr(get(hObject,'String'));
handles.problem.N = str2double(contents{get(hObject,'Value')});
guidata(hObject, handles)
updateTimings(hObject, handles);


% --- Executes during object creation, after setting all properties.
function discretizationStepsDropDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to discretizationStepsDropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addDiscretizationStepPushButton.
function addDiscretizationStepPushButton_Callback(hObject, eventdata, handles)
% hObject    handle to addDiscretizationStepPushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.problem.useMaxNorm == true
    normStr = 'Max norm';
else
    normStr = 'Euclidean norm';
end
prompt = ['Number of discretization steps (' normStr '):'];
titleString = '(SLOW) Add discretization step';
N = str2double(inputdlg(prompt, titleString));
if ~isempty(N)
    createConstraintGradientFunction(N,handles.problem.useMaxNorm)
    precomputedDiscretizations = getPrecomputedDiscretizations(hObject, handles);
    set(handles.discretizationStepsDropDown, 'String', precomputedDiscretizations)
    guidata(hObject, handles)
end

function plotResult(hObject, handles)
obj = get(handles.plotButtonGroup,'SelectedObject');
switch get(obj,'Tag') % Get Tag of selected object.
    case 'qRadioButton'
        plotq(hObject, handles);
    case 'gRadioButton'
        plotGradient(hObject, handles);
    case 'slewRadioButton'
        plotSlew(hObject, handles);
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

outfile_1st = fopen([out_dir filesep out_name], 'w');

formatspec = '%8.5f %8.5f %8.5f\r\n';

out_mat_1st = result.g';

fprintf(outfile_1st, '%8.0i\r\n', size(out_mat_1st, 2));
fprintf(outfile_1st, formatspec, out_mat_1st);


fclose(outfile_1st);



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
set(handles.encodingTensorTable, 'Data', handles.output(index).problem.targetTensor)
set(handles.slewRateTextBox, 'String', num2str(handles.output(index).problem.sMax))
set(handles.gMaxTextBox, 'String', num2str(handles.output(index).problem.gMax))
set(handles.maxNormRadioButton, 'Value', handles.output(index).problem.useMaxNorm);
set(handles.EuclideanNormRadioButton, 'Value', ~handles.output(index).problem.useMaxNorm);

precomputedDiscretizations = getPrecomputedDiscretizations(hObject, handles.output(index));
set(handles.discretizationStepsDropDown, 'String', precomputedDiscretizations)
set(handles.discretizationStepsDropDown, 'Value', getIndexOfSameN(handles.output(index).problem, precomputedDiscretizations))

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
set(handles.redoIfFailedCheckBox,'Value', handles.output(index).problem.redoIfFailed)
set(handles.heatDissipationTextBox, 'String', num2str(handles.output(index).problem.eta))
set(handles.nameTextBox, 'String', handles.output(index).problem.name)

function index = getIndexOfSameN(problem, precomputedDiscretizations)
index = find(problem.N == str2double(precomputedDiscretizations), 1);
% for i = 1:length(precomputedDiscretizations)
%     if problem.N == str2double(precomputedDiscretizations{i})
%         set(handles.discretizationStepsDropDown, 'Value', i)
%         break
%     end
% end

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

% function index = findIndexOfName(hObject, handles, name)
% for i = 1:length(handles.output)
%     if strcmp(handles.output(i).name,
% end


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
