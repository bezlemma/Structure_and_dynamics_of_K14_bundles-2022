%% Setup
global umPerPixel;
global yBin; %Bin size for y-data
global xRgt;
global yBotEnd;
global TPH_TOP;
global TPH_BOT;
sensitivity = 30; graph = 0;
umPerPixel = 116/1024;
pth_sdt = 'C:\Data\';
% Set up frames list and time list

%Data info
yBin =40; trim = 30; 
FrameStart = 20;  FrameEnd = 50;
frames = (FrameStart:FrameEnd) - (FrameStart -1);

%% Load Time Stamp information
%time_list_path = [pth_sdt '\time_list'];
%or if the framerate is constant, you can just do something like:
%time_list = frames*3;


%% Top Bleach
tph_name = 'TopBleach\';
[TPH_TOP, data1D] = func_TPH_read(pth_sdt, tph_name, frames,FrameStart);
xRgt = length(data1D(:,1,1)); yBotEnd = length(TPH_TOP(:,1,1));
[SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid] =...
    func_guesspeaks(data1D,TPH_TOP,frames,sensitivity,graph,trim);
[lftloc_TOP,rgtloc_TOP,lftamp_TOP,rgtamp_TOP,lftwid_TOP,rgtwid_TOP,lftamp_err_TOP,rgtamp_err_TOP,~] =  func_doubleline_fit( ...
    data1D,frames,TPH_TOP,graph,SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid);

%% Bot Bleach
tph_name = 'BotBleach\';
[TPH_BOT, data1D] = func_TPH_read(pth_sdt, tph_name, frames,FrameStart);
xRgt = length(data1D(:,1,1)); yBotEnd = length(TPH_BOT(:,1,1));
[SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid] =...
    func_guesspeaks(data1D,TPH_BOT,frames,sensitivity,graph,trim);
[lftloc_BOT,rgtloc_BOT,lftamp_BOT,rgtamp_BOT,lftwid_BOT,rgtwid_BOT,lftamp_err_BOT,rgtamp_err_BOT,~] =  func_doubleline_fit( ...
    data1D,frames,TPH_BOT,graph,SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid);

%% Calculate other things
PAD = length(TPH_BOT(1,:,1));
BOT_MIDLINE = (lftloc_BOT+rgtloc_BOT)/2;
TOP_MIDLINE = (lftloc_TOP+rgtloc_TOP)/2 + PAD;

%% Movie
TPH = [TPH_BOT TPH_TOP];
low = min(min(min(TPH))) + 1;
high = max(max(max(TPH))) - 1;

yBotEnd = length(TPH_BOT(:,1,1));

ylist = (1/2):(yBotEnd-yBin);
for fr = frames(end):-1:frames(1)
    figure();
    imshow(TPH(:,:,fr),[low high])
    hold on;
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'visible','off')
   
    %Plot BOT
    plot(lftloc_BOT(fr,1:10:end),  ylist(1:10:end)+yBin/2,'k.', 'MarkerSize',12 );
    plot(rgtloc_BOT(fr,1:10:end),  ylist(1:10:end)+yBin/2,'k.', 'MarkerSize',12 );
    plot(BOT_MIDLINE(fr,:),  ylist+yBin/2,'b-', 'LineWidth', 4 );
        
    %Plot TOP
    plot(PAD+lftloc_TOP(fr,1:10:end),  ylist(1:10:end)+yBin/2,'k.', 'MarkerSize',12 );
    plot(PAD+rgtloc_TOP(fr,1:10:end),  ylist(1:10:end)+yBin/2,'k.', 'MarkerSize',12 );
    plot(TOP_MIDLINE(fr,:),  ylist+yBin/2,'b-', 'LineWidth', 4 );
    
    exportpath = [pth_sdt '\Movie' num2str(FrameStart) '.tiff'];
    export_fig C:\Data\test.tif -append;
    
    close all;
end
