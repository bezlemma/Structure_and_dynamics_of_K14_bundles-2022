%% Setup
global umPerPixel;
global yBin; %Bin size for y-data
global xRgt;
global yTopEnd;
global TPH_TOP;
sensitivity = 30; graph = 0;
umPerPixel = 1/8.8;
% Set up frames list and time list
pth_sdt = 'C:\Data\';

%1 PEG
yBin =20; trim = 30; 
FrameStart = 1;  FrameEnd = 78;
frames = (FrameStart:FrameEnd) - (FrameStart -1);

%% Load Time Stamp information
%time_list_path = [pth_sdt '\time_list'];
%or if the framerate is constant, you can just do something like:
%time_list = frames*3;


%% Top Bleach
tph_name = '1PEG_Split2\';
[TPH_TOP, data1D] = func_TPH_read(pth_sdt, tph_name, frames,FrameStart);
xRgt = length(data1D(:,1,1)); 
yTopEnd = length(TPH_TOP(:,1,1));
[SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid] = func_guesspeaks(data1D,TPH_TOP,frames,sensitivity,graph,trim);
[lftloc_TOP,rgtloc_TOP,lftamp_TOP,rgtamp_TOP,lftwid_TOP,rgtwid_TOP,lftamp_err_TOP,rgtamp_err_TOP,~] =  func_doubleline_fit_forward( ...
    data1D,frames,TPH_TOP,graph,SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid);

%% Calculate other things
PAD = length(TPH_TOP(1,:,1));

%% Movie
TPH = [TPH_TOP];
low = min(min(min(TPH))) + 1;
high = max(max(max(TPH))) - 1;

yTopEnd = length(TPH_TOP(:,1,1));

ylist = (1/2):(yTopEnd-yBin);
for fr = frames(1):frames(end)
    figure();
    imshow(TPH(:,:,fr),[low high])
    hold on;
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'visible','off')
        
    %Plot TOP
    plot(lftloc_TOP(fr,1:10:end),  ylist(1:10:end)+yBin/2,'b.', 'MarkerSize',12 );
    plot(rgtloc_TOP(fr,1:10:end),  ylist(1:10:end)+yBin/2,'b.', 'MarkerSize',12 );
    
    exportpath = [pth_sdt '\Movie' num2str(FrameStart) '.tiff'];
    export_fig C:\Data\test.tif -append;
  %pause  
    close all;
end
