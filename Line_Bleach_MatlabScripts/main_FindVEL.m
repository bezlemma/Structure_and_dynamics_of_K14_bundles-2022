%Size per Pixel
global umPerPixel;  
global yBin; %Bin size for y-data
yBin = 50;
global xRgt;
global yBotEnd;
global TPH;

%% Particular data information
 trim = 30;
 sensitivity = 30;
 FrameStart = 20; 
 FrameEnd = 50;
 pth_sdt = 'C:\Data\';
 tph_name = '\BleachData\';
 umPerPixel = 116/1024;   


%% Set up frames list with sane numbering
frames = (FrameStart:FrameEnd) - (FrameStart -1);

%% Read in TPH data and flatten into 1D strips
[TPH, data1D] = func_TPH_read(pth_sdt, tph_name, frames,FrameStart);
xRgt = length(data1D(:,1,1));
yBotEnd = length(TPH(:,1,1));

%% Determine a guesss set of Peaks
graph = 1;
[SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid] =...
    func_guesspeaks(data1D,TPH,frames,sensitivity,graph,trim);

%% Fit double gaussians to 1D strips
graph = 0;
[lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,lftamp_err,rgtamp_err,SkipList] =  func_doubleline_fit( ...
    data1D,frames,TPH,graph,SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid);

%% Load Time Stamp information
%time_list_path = [pth_sdt '\time_list'];
%or if the framerate is constant, you can just do something like:
%time_list = frames*3;

%% Calculate velocities
graph = 0;
[V_list_L,V_list_R, vErrL, vErrR, CoM] = func_vel_fit_list_CoM(frames, lftloc,rgtloc,SkipList,time_list,graph);
%V_list_sum = (V_list_R - V_list_L);
V_list_sum = (abs(V_list_R) + abs(V_list_L));

%% Plot Velocity over time...
for i=5:length(frames)
    [V_list_L,V_list_R, vErrL, vErrR] = func_vel_fit_list_CoM(frames(1:i), lftloc,rgtloc,SkipList,time_list,graph);
    V_list_sum = (abs(V_list_R) + abs(V_list_L));  
    V_Scatter(i,:) = V_list_sum;
    Vmean(i) = nanmean(V_list_sum);
    Vstderr(i) =  nanstd(V_list_sum) / sqrt(length(V_list_sum));
end

figure1 = figure;
axes1 = axes('Parent',figure1);
hold on
xlabel('Time, (s)','FontSize',20,'FontName','Montserrat','interpreter','latex');
set(axes1,'FontName','Montserrat','FontSize',14,'XGrid','on');
ylabel('Mean Velocity (nm/s)','FontSize',20,'FontName','Montserrat');
errorbar(time_list(1:end),Vmean(1:end),Vstderr(1:end),'r-','MarkerSize',25, 'LineWidth',2)
plot(time_list(5:end),V_Scatter(1:end,:),'.')
clear V_Scatter

%% Movie that shows position, w/ color for (from blue to red) velocity.
low = min(min(min(TPH))) + 1;
high = max(max(max(TPH))) - 1;

ylist = (1/2):(yBotEnd-yBin);
for fr = frames(end):-1:frames(1)
figure();
imshow(TPH(:,:,fr),[low high]) 
hold on;
set(gca,'position',[0 0 1 1],'units','normalized')
set(gca,'visible','off')
  plot(lftloc(fr,1:10:end),  ylist(1:10:end)+yBin/2,'k.', 'MarkerSize',12 );
  plot(rgtloc(fr,1:10:end),  ylist(1:10:end)+yBin/2,'k.', 'MarkerSize',12 );
  plot( (lftloc(fr,:)+rgtloc(fr,:))/2,  ylist+yBin/2,'b-', 'LineWidth', 4 );
  exportpath = [pth_sdt '\Movie' num2str(FrameStart) '.tiff'];
  export_fig C:\Data.tif -append;

close all;
end