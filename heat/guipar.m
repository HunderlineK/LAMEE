function fig = guipar()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load guipar

h0 = figure('BackingStore','off', ...
	'Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','C:\Program\Matlab53\m-filer\lab\guipar.m', ...
	'MenuBar','none', ...
	'Name','Simulator Parameters', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[530 162 369 201], ...
	'Tag','Fig2', ...
	'ToolBar','none');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat1, ...
	'ListboxTop',0, ...
	'Position',[288 244.5 32.25 13.5], ...
	'String','2000', ...
	'Style','edit', ...
	'Tag','k');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[231 239.25 54 21.75], ...
	'String','Heat Transfer Coefficient', ...
	'Style','text', ...
	'Tag','StaticText6');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[10 10 38.25 19.5], ...
	'String','Holdup volume', ...
	'Style','text', ...
	'Tag','StaticText7');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[10 10 34.5 21], ...
	'String','No of nodes', ...
	'Style','text', ...
	'Tag','StaticText8');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[10 10 32.25 13.5], ...
	'String','20', ...
	'Style','edit', ...
	'Tag','n');
if nargout > 0, fig = h0; end
