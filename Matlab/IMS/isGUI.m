%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Basic Image Source Model
%Adapted from that taught by Adam Hill
%Built upon to include dual Sources 
%Soon to finish implemented re-additive bass after pre-filtering 
%Then may add binaural crosstalk/cancellation 
%
%Simon Durbridge
%LSSDO
%MSc Audio Engineering 2016
%
%
%Run isGUI to start the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varargout = isGUI(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @isGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @isGUI_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPENING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isGUI_OpeningFcn(hObject, eventdata, handles, varargin)

%Set Absorption Coefs
handles.marble = 0.005;
handles.wood = 0.1;
handles.heavyFabric = 0.45;
handles.lightFabric = 0.14;
handles.glass = 0.03;
handles.plaster = 0.05;
handles.tiles = 0.03;
handles.carpet = 0.37;

%Set Reflection Coefs
handles.alphaxpos = 1 - handles.plaster;
handles.alphaxneg = 1 - handles.plaster;
handles.alphaypos = 1 - handles.heavyFabric;
handles.alphayneg = 1 - handles.plaster;
handles.alphazpos = 1 - handles.tiles;
handles.alphazneg = 1 - handles.carpet;

handles.filtord = 2;
handles.lpfco = 10000;
handles.hppford = 4;
handles.hppfco = 250;
%Do Filtering for HF loss due to distance

handles.Fs = 48000; %sample rate
handles.output = hObject;
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function varargout = isGUI_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT TEXT BOXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RSX_Callback(hObject, eventdata, handles)
input = str2num(get(handles.RSX,'String'));
if input < 1
    input = 1;
end
set(handles.RSX, 'Value', input);
set(handles.RSX, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function RSX_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RSY_Callback(hObject, eventdata, handles)
input = str2num(get(handles.RSY,'String'));
if input < 1
    input = 1;
end
set(handles.RSY, 'Value', input);
set(handles.RSY, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function RSY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RSZ_Callback(hObject, eventdata, handles)
input = str2num(get(handles.RSZ,'String'));
if input < 1
    input = 1;
end
set(handles.RSZ, 'Value', input);
set(handles.RSZ, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function RSZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LLX_Callback(hObject, eventdata, handles)
input = str2num(get(handles.LLX,'String'));
if input < 1
    input = 1;
end
set(handles.LLX, 'Value', input);
set(handles.LLX, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function LLX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LLY_Callback(hObject, eventdata, handles)
input = str2num(get(handles.LLY,'String'));
if input < 1
    input = 1;
end
set(handles.LLY, 'Value', input);
set(handles.LLY, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function LLY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LLZ_Callback(hObject, eventdata, handles)
input = str2num(get(handles.LLZ,'String'));
if input < 1
    input = 1;
end
set(handles.LLZ, 'Value', input);
set(handles.LLZ, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function LLZ_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SLX1_Callback(hObject, eventdata, handles)
input = str2num(get(handles.SLX1,'String'));
if input < 1
    input = 1;
end
set(handles.SLX1, 'Value', input);
set(handles.SLX1, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function SLX1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SLY1_Callback(hObject, eventdata, handles)
input = str2num(get(handles.SLY1,'String'));
if input < 1
    input = 1;
end
set(handles.SLY1, 'Value', input);
set(handles.SLY1, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function SLY1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SLZ1_Callback(hObject, eventdata, handles)
input = str2num(get(handles.SLZ1,'String'));
if input < 1
    input = 1;
end
set(handles.SLZ1, 'Value', input);
set(handles.SLZ1, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function SLZ1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SLX2_Callback(hObject, eventdata, handles)
input = str2num(get(handles.SLX2,'String'));
if input < 1
    input = 1;
end
set(handles.SLX2, 'Value', input);
set(handles.SLX2, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function SLX2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SLY2_Callback(hObject, eventdata, handles)
input = str2num(get(handles.SLY2,'String'));
if input < 1
    input = 1;
end
set(handles.SLY2, 'Value', input);
set(handles.SLY2, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function SLY2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SLZ2_Callback(hObject, eventdata, handles)
input = str2num(get(handles.SLZ2,'String'));
if input < 1
    input = 1;
end
set(handles.SLZ2, 'Value', input);
set(handles.SLZ2, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function SLZ2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MRO_Callback(hObject, eventdata, handles)
input = str2num(get(handles.MRO,'String'));
if input < 1
    input = 1;
end
set(handles.MRO, 'Value', input);
set(handles.MRO, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
function MRO_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HFOR_Callback(hObject, eventdata, handles)
input = str2num(get(handles.HFOR,'String'));
if input < 1
    input = 1;
end
set(handles.HFOR, 'Value', input);
set(handles.HFOR, 'String', num2str(input));
handles.filtord = get(handles.HFOR, 'Value');
guidata(hObject, handles);
function HFOR_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HFCO_Callback(hObject, eventdata, handles)
input = str2num(get(handles.HFCO,'String'));
if input < 50
    input = 50;
end
set(handles.HFCO, 'Value', input);
set(handles.HFCO, 'String', num2str(input));
handles.lpfco = get(handles.HFCO, 'Value');
guidata(hObject, handles);
function HFCO_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HPPFOrd_Callback(hObject, eventdata, handles)
input = str2num(get(handles.HPPFOrd,'String'));
if input < 1
    input = 1;
end
set(handles.HPPFOrd, 'Value', input);
set(handles.HPPFOrd, 'String', num2str(input));
handles.hppford = get(handles.HPPFOrd, 'Value');
guidata(hObject, handles);
function HPPFOrd_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HPPFco_Callback(hObject, eventdata, handles)
input = str2num(get(handles.HPPFco,'String'));
if input < 20
    input = 20;
end
set(handles.HPPFco, 'Value', input);
set(handles.HPPFco, 'String', num2str(input));
handles.hppfco = get(handles.HPPFco, 'Value');
guidata(hObject, handles);
function HPPFco_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MATERIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XPMAT_Callback(hObject, eventdata, handles)
choice = get(handles.XPMAT, 'Value');
switch choice
    case 1 % marble
        temp = 1 - handles.marble;
    case 2 % wood
        temp = 1 - handles.wood;
    case 3 % heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 % light fabric
        temp = 1 - handles.lightFabric;
    case 5 % glass
        temp = 1 - handles.glass;
    case 6 % plaster
        temp = 1 - handles.plaster;
    case 7 % tiles
        temp = 1 - handles.tiles;
    case 8 % carpet
        temp = 1 - handles.carpet;
end
handles.alphaxpos = temp;
guidata(hObject, handles);
function XPMAT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function XNMAT_Callback(hObject, eventdata, handles)
choice = get(handles.XNMAT, 'Value');
switch choice
    case 1 % marble
        temp = 1 - handles.marble;
    case 2 % wood
        temp = 1 - handles.wood;
    case 3 % heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 % light fabric
        temp = 1 - handles.lightFabric;
    case 5 % glass
        temp = 1 - handles.glass;
    case 6 % plaster
        temp = 1 - handles.plaster;
    case 7 % tiles
        temp = 1 - handles.tiles;
    case 8 % carpet
        temp = 1 - handles.carpet;
end
handles.alphaxneg = temp;
guidata(hObject, handles);
function XNMAT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function YPMAT_Callback(hObject, eventdata, handles)
choice = get(handles.YPMAT, 'Value');
switch choice
    case 1 % marble
        temp = 1 - handles.marble;
    case 2 % wood
        temp = 1 - handles.wood;
    case 3 % heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 % light fabric
        temp = 1 - handles.lightFabric;
    case 5 % glass
        temp = 1 - handles.glass;
    case 6 % plaster
        temp = 1 - handles.plaster;
    case 7 % tiles
        temp = 1 - handles.tiles;
    case 8 % carpet
        temp = 1 - handles.carpet;
end
handles.alphaypos = temp;
guidata(hObject, handles);
function YPMAT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function YNMAT_Callback(hObject, eventdata, handles)
choice = get(handles.YNMAT, 'Value');
switch choice
    case 1 % marble
        temp = 1 - handles.marble;
    case 2 % wood
        temp = 1 - handles.wood;
    case 3 % heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 % light fabric
        temp = 1 - handles.lightFabric;
    case 5 % glass
        temp = 1 - handles.glass;
    case 6 % plaster
        temp = 1 - handles.plaster;
    case 7 % tiles
        temp = 1 - handles.tiles;
    case 8 % carpet
        temp = 1 - handles.carpet;
end
handles.alphayneg = temp;
guidata(hObject, handles);
function YNMAT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ZPMAT_Callback(hObject, eventdata, handles)
choice = get(handles.ZPMAT, 'Value');
switch choice
    case 1 % marble
        temp = 1 - handles.marble;
    case 2 % wood
        temp = 1 - handles.wood;
    case 3 % heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 % light fabric
        temp = 1 - handles.lightFabric;
    case 5 % glass
        temp = 1 - handles.glass;
    case 6 % plaster
        temp = 1 - handles.plaster;
    case 7 % tiles
        temp = 1 - handles.tiles;
    case 8 % carpet
        temp = 1 - handles.carpet;
end
handles.alphazpos = temp;
guidata(hObject, handles);
function ZPMAT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ZNMAT_Callback(hObject, eventdata, handles)
choice = get(handles.ZNMAT, 'Value');
switch choice
    case 1 % marble
        temp = 1 - handles.marble;
    case 2 % wood
        temp = 1 - handles.wood;
    case 3 % heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 % light fabric
        temp = 1 - handles.lightFabric;
    case 5 % glass
        temp = 1 - handles.glass;
    case 6 % plaster
        temp = 1 - handles.plaster;
    case 7 % tiles
        temp = 1 - handles.tiles;
    case 8 % carpet
        temp = 1 - handles.carpet;
end
handles.alphazneg = temp;
guidata(hObject, handles);
function ZNMAT_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BUTTONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DoCalc_Callback(hObject, eventdata, handles)
%Do Calculation
handles = IMS(handles);
%Tempsave IRs
tempirl = handles.irl;
tempirr = handles.irr;
%Extend shortest IR for plotting if unequal lengths
if length(tempirl) > length(tempirr)
   tempirr((length(tempirr)+1):length(tempirl)) = zeros(1,(length(tempirl) - length(tempirr)));
elseif length(tempirl) < length(tempirr)
   tempirl((length(tempirl)+1):length(tempirr)) = zeros(1,(length(tempirr) - length(tempirl))); 
end
%Create time vector and plot
set(handles.text27, 'Visible', 'on');
axes(handles.axes2), set(handles.axes2, 'Visible', 'On');
cla(handles.axes2), reset(handles.axes2), hold off;
tl = linspace(0,length(tempirl)/handles.Fs, length(tempirl));
tr = linspace(0,length(tempirr)/handles.Fs, length(tempirr));
plot(tl, tempirl, 'b', tr, tempirr, 'r'), axis tight, grid on;
xlabel('Time (s)'), ylabel('amplitude');
guidata(hObject, handles);

function PLIR_Callback(hObject, eventdata, handles)

%combine impulse responses for stereo
irtl = handles.irl;
irtr = handles.irr;
if length(irtl) > length(irtr)
irtl = irtl(1:length(irtr));
else
irtr = irtr(1:length(irtl)); 
end
impresplayer(:,1) = irtl;
impresplayer(:,2) = irtr;

impresplayer = audioplayer(impresplayer, handles.Fs);
playblocking(impresplayer);
guidata(hObject, handles);

function SVIR_Callback(hObject, eventdata, handles)
irtl = handles.irl;
irtr = handles.irr;
if length(irtl) > length(irtr)
irtl = irtl(1:length(irtr));
else
irtr = irtr(1:length(irtl)); 
end
impres(:,1) = irtl;
impres(:,2) = irtr;
audiowrite('IMG_SRC_RESP.wav', impres, handles.Fs);

function LoadAudio_Callback(hObject, eventdata, handles)
[handles.audio audioFs] = getFile(' the audio sample');
if ~isequal(audioFs, handles.Fs)
    handles.audio = resample(handles.audio, handles.Fs, audioFs);
end

guidata(hObject, handles);

function PlayDry_Callback(hObject, eventdata, handles)
audiotemp = handles.audio;
dryplayer = audioplayer(audiotemp, handles.Fs);
playblocking(dryplayer);
guidata(hObject, handles);

function PlayWet_Callback(hObject, eventdata, handles)
audiotemp = handles.audio;
audiotl = handles.audio(1,:);
audiotr = handles.audio(2,:);
%Calculate Prefilter
[B2, A2] = butter(handles.hppford, handles.hppfco/(handles.Fs/2),'high');
[B3, A3] = butter(handles.hppford, handles.hppfco/(handles.Fs/2),'low');
%Do Filtering for LF rejection due to basic modeling method
audiotl = filter(B2, A2, audiotl);
audiotr = filter(B2, A2, audiotr);
audiolpl = filter(B3, A3, audiotl);
audiolpr = filter(B3, A3, audiotr); 
%Calculate Filter 
[B1, A1] = butter(handles.filtord, handles.lpfco/(handles.Fs/2),'low');
%Do Filtering for HF loss due to distance
irlt = filter(B1, A1, handles.irl);
irrt = filter(B1, A1, handles.irr);

%Do Conv
audioL = conv(audiotl,irlt);
audioR = conv(audiotr,irrt);



% audioL = audioL + audiolpl;
% audioR = audioR + audiolpr;

%combine tracks
if length(audioL) > length(audioR)
audioL = audioL(1:length(audioR));
else
audioR = audioR(1:length(audioL)); 
end
audioout(:,1) = audioL;
audioout(:,2) = audioR;

wetplayer = audioplayer(audioout, handles.Fs);
playblocking(wetplayer);
guidata(hObject, handles);

function SaveWet_Callback(hObject, eventdata, handles)
audiotemp = handles.audio;
audiotl = handles.audio(1,:);
audiotr = handles.audio(2,:);
%Calculate Filter 
[B1, A1] = butter(handles.filtord, handles.lpfco/(handles.Fs/2),'low');
%Do Filtering for HF loss due to distance
irlt = filter(B1, A1, handles.irl);
irrt = filter(B1, A1, handles.irr);

%Do Conv
audioL = conv(audiotl,irlt);
audioR = conv(audiotr,irrt);

%combine tracks
if length(audioL) > length(audioR)
audioL = audioL(1:length(audioR));
else
audioR = audioR(1:length(audioL)); 
end
audioout(:,1) = audioL;
audioout(:,2) = audioR;

audiowrite('IMG_SRC_processedaudio2.wav', audioout, handles.Fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOTTING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotter_Callback(hObject, eventdata, handles)
axes(handles.axes1), set(handles.axes1, 'Visible', 'On');
cla(handles.axes1),reset(handles.axes1), hold off;
plot(get(handles.SLX1, 'Value'), get(handles.SLY1, 'Value'), ...
    'k', 'Marker', 'o', 'MarkerSize', 10), hold on;
plot(get(handles.SLX2, 'Value'), get(handles.SLY2, 'Value'), ...
    'k', 'Marker', 's', 'MarkerSize', 10);
plot(get(handles.LLX, 'Value'), get(handles.LLY, 'Value'), ...
    'k', 'Marker', 'x', 'MarkerSize', 10), hold off;
grid on, axis([0 get(handles.LLX, 'Value') 0 get(handles.LLY, 'Value')]);
axis([0 get(handles.RSX, 'Value') 0 get(handles.RSY, 'Value')]);
xlabel('Left/Right'), ylabel('Front/Back');
title(' o = source, x = listener');
