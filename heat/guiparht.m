function fig = guiparht()
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

load guiparht

h0 = figure('Units','points', ...
	'BackingStore','off', ...
	'Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','C:\Program\Matlab53\m-filer\lab\guiparht.m', ...
	'MenuBar','none', ...
	'Name','Heat Transfer Parameters', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[269.25 89.25 375 262.5], ...
	'Tag','Fig1', ...
	'ToolBar','none');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','setthispar;', ...
	'ListboxTop',0, ...
	'Position',[60 192.75 45 15], ...
	'Style','edit', ...
	'Tag','coldp');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','setthispar;', ...
	'ListboxTop',0, ...
	'Position',[59.25 165.75 45 15], ...
	'Style','edit', ...
	'Tag','hotp');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','setthispar;', ...
	'ListboxTop',0, ...
	'Position',[57 141 45 15], ...
	'Style','edit', ...
	'Tag','k');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[222.75 33 45 15], ...
	'Style','edit', ...
	'Tag','Vh');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[9 71.25 50.25 22.5], ...
	'String','Volume per area', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[6.75 189.75 45 15], ...
	'String','Cold Passes', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[6 163.5 45 15], ...
	'String','Hot Passes', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[11.25 138.75 45 15], ...
	'String','k', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[176.25 33 45 15], ...
	'String','V Hot', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','setthispar;', ...
	'ListboxTop',0, ...
	'Position',[58.5 75 45 15], ...
	'Style','edit', ...
	'Tag','vpera');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[175.5 13.5 45 15], ...
	'String','V Cold', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[224.25 14.25 45 15], ...
	'Style','edit', ...
	'Tag','Vc');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[10.5 117 45 15], ...
	'String','Plate Area', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','setthispar;', ...
	'ListboxTop',0, ...
	'Position',[55.5 119.25 45 15], ...
	'Style','edit', ...
	'Tag','platea');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','setthispar;', ...
	'ListboxTop',0, ...
	'Position',[55.5 99 45 15], ...
	'Style','edit', ...
	'Tag','nplates');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[9.75 96.75 45 15], ...
	'String','Plates', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[175.5 52.5 45 29.25], ...
	'String','Heat Transfer Area', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[222 57.75 45 15], ...
	'Style','edit', ...
	'Tag','A');
if nargout > 0, fig = h0; end
