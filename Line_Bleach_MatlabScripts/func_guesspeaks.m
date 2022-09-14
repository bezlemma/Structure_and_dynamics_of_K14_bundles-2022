function [SkipList,lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,avgwid]  = ...
    func_guesspeaks(data1D,TPH,frames,sensitivity,graph,trim)

global yBin
%Height of data
global yBotEnd
global xRgt

%% ---------------- Fitting parameters
trimStart = trim;
trimEnd = xRgt-trim;
clear rgtloc;
clear lftloc;

%Intialize/Reset Variables
loc   = NaN*zeros(length(frames),yBotEnd-yBin,2);
lftloc   = NaN*zeros(length(frames),yBotEnd-yBin);
rgtloc   = NaN*zeros(length(frames),yBotEnd-yBin);
wid   = NaN*zeros(length(frames),yBotEnd-yBin,2);
lftwid   = NaN*zeros(length(frames),yBotEnd-yBin);
rgtwid   = NaN*zeros(length(frames),yBotEnd-yBin);
amp  = NaN*zeros(length(frames),yBotEnd-yBin,2);
lftamp   = NaN*zeros(length(frames),yBotEnd-yBin);
rgtamp   = NaN*zeros(length(frames),yBotEnd-yBin);
SkipList =  zeros(yBotEnd-yBin,length(frames));
ylist = (1/2):(yBotEnd-yBin);
yIntList = 1:(yBotEnd-yBin) ;


%% Make an intial set of guesses by flattening the data to a single strip
for fr = (frames(end)):-1:frames(1)
    
    FlatFrameData = mean(data1D(trimStart:trimEnd,1:end,fr),2) ;
    
    %subtract of a linear bit from this
     FlatFrameData =  detrend( FlatFrameData);
    %FrameData = data1D(trimStart:trimEnd,yy,fr) ;
    
    [amp_pks, locs_pks, wid_pks, prom_pks] = ...
    findpeaks(smooth(FlatFrameData),'MinPeakProminence',max(FlatFrameData)/(sensitivity), ...
    'MinPeakHeight',mean(FlatFrameData));
    [prom_pks, a_order] = sort(prom_pks,'descend');
    locs_pks =  locs_pks(a_order,:);
    
    if length(amp_pks) < 2
        disp([ 'Aborted at ' num2str(fr)])
        [amp_pks, locs_pks, wid_pks, prom_pks] = ...
    findpeaks(FlatFrameData,'MinPeakProminence',max(FlatFrameData)/(sensitivity*2), ...
    'MinPeakHeight',mean(FlatFrameData));
        [prom_pks, a_order] = sort(prom_pks,'descend');
        locs_pks =  locs_pks(a_order,:);
        
        if length(amp_pks) < 2
            disp([ 'Aborted at ' num2str(fr)])
            [amp_pks, locs_pks, wid_pks, prom_pks] = ...
    findpeaks(FlatFrameData,'MinPeakProminence',max(FlatFrameData)/(sensitivity*4), ...
    'MinPeakHeight',mean(FlatFrameData));
            [prom_pks, a_order] = sort(prom_pks,'descend');
            locs_pks =  locs_pks(a_order,:);
            if length(amp_pks) < 2
                disp([ 'Aborted at ' num2str(fr)])
                fuck
            end
        end
        
    end
    
    
    
    for yy = yIntList
            loc(fr,yy,2) = locs_pks(2)+trimStart-1;
            loc(fr,yy,1) = locs_pks(1)+trimStart-1;
    end
end



%% Sort out our data into a left column and a right column

LeftMaster  = (loc(:,:,1) > loc(:,:,2)) .* (loc(:,:,1) > 0);
RightMaster = (loc(:,:,1) < loc(:,:,2)) .* (loc(:,:,2) > 0);

lftloc =    LeftMaster.*loc(:,:,2) + RightMaster.*loc(:,:,1);
rgtloc =   LeftMaster.*loc(:,:,1) + RightMaster.*loc(:,:,2);

%% Possible Graphing Step
if graph == 1
    
    for fr = frames(end):-1:frames(1)
        figure('Name',['Frame: ' num2str(fr)],'NumberTitle','off');
        imagesc(TPH(:,:,fr))
        colormap('gray')
        hold on;
        plot(rgtloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0)+yBin/2,'r.','MarkerSize',15)
        plot(lftloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0)+yBin/2,'b.','MarkerSize',15)
        
        title(['Raw Guess ' num2str(fr) ' Fit'],'FontSize',18,'interpreter','latex');
        ylabel('Bleach Axis [um]','FontSize',18,'interpreter','latex');
        xlabel('Ordered Axis[um]','FontSize',18,'interpreter','latex');
        set(gca,'fontsize',18)
        
    end
end



%Fit slopes to the lines
for fr = frames
    %3rd order
    %BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),3);
    %BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),3);
    %predictR(fr,:) = BleachSlopeR(4) + BleachSlopeR(3)*yIntList + BleachSlopeR(2)*yIntList.^2 + BleachSlopeR(1)*yIntList.^3;
    %predictL(fr,:) = BleachSlopeL(4) + BleachSlopeL(3)*yIntList + BleachSlopeL(2)*yIntList.^2 + BleachSlopeL(1)*yIntList.^3;
    
    %1st order, for straight data
    BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),1);
    BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),1);
    predictR(fr,:) = BleachSlopeR(2) + BleachSlopeR(1)*yIntList;
    predictL(fr,:) = BleachSlopeL(2) + BleachSlopeL(1)*yIntList;
    
    %0th order
%     BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),0);
%     BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),0);
%     predictR(fr,:) = BleachSlopeR(1)*yIntList./yIntList;
%     predictL(fr,:) = BleachSlopeL(1)*yIntList./yIntList;
    
    
end


%% Make a second set of guesses with trim calculated from the first set
% 
for fr = (frames(end)):-1:frames(1)
    for yy = yIntList
        SkipList(yy,fr) = 0;
        lftloc(fr,yy) = NaN;
        rgtloc(fr,yy) = NaN;
        
        buffer = 10;
        
        trimStart = max(floor(predictL(fr,yy) - buffer),1);
        idxValidL = min(trimStart:ceil(predictL(fr,yy) + buffer),length(data1D(:,1,1)));
        idxValidL = idxValidL(idxValidL > 0);
        if idxValidL > buffer
        FrameData = data1D(idxValidL,yy,fr) ;
        else
        FrameData = data1D(trimStart:floor(xRgt/2),yy,fr) ;
        end
        
        %Find the peaks
        [amp_pks, locs_pks, wid_pks, prom_pks] = ...
            findpeaks(smooth(FrameData),'MinPeakProminence',max(FrameData)/sensitivity);
       
        if ~isempty(amp_pks)
            [prom_pks, a_order] = sort(prom_pks,'descend');
            locs_pks =  locs_pks(a_order,:);
            wid_pks =  wid_pks(a_order,:);
            amp_pks =  amp_pks(a_order,:);
            
            lftwid(fr,yy) = wid_pks(1);
            lftamp(fr,yy) = amp_pks(1);
            lftloc(fr,yy) = locs_pks(1)+trimStart-1; %Check about removing the 1 to see if this is wrong
            SkipList(yy,fr) = 0;
        end
              
        trimStart = max(floor(predictR(fr,yy) - buffer),1);
        idxValidR = min(trimStart:ceil(predictR(fr,yy) + buffer),length(data1D(:,1,1)));
        %Find the peaks
        if idxValidR > buffer
        FrameData = data1D(idxValidL,yy,fr) ;
        else
        FrameData = data1D(floor(xRgt/2):xRgt,yy,fr) ;
        end
        
        [amp_pks, locs_pks, wid_pks, prom_pks] = ...
            findpeaks(smooth(FrameData),'MinPeakProminence',max(FrameData)/sensitivity);
   
        if ~isempty(amp_pks)
        [prom_pks, a_order] = sort(prom_pks,'descend');
        locs_pks =  locs_pks(a_order,:);
        wid_pks =  wid_pks(a_order,:);
        amp_pks =  amp_pks(a_order,:);
        
        rgtwid(fr,yy) =  wid_pks(1);
        rgtamp(fr,yy) =  amp_pks(1);
        rgtloc(fr,yy) = locs_pks(1)+trimStart-1; %Check about removing the 1 to see if this is wrong
        SkipList(yy,fr) = 0;
        end
        
    end
end
% 
 avgwid = nanmean(nanmean(lftwid+rgtwid))/2;


%% Check to see if velocities are sane

%Fit a velocity to the motion
for yy = yIntList
    LinFitsR = fitlm( frames , rgtloc(frames,yy) );
    vFitsR(yy) = LinFitsR.Coefficients.Estimate(2);
    LinFitsL = fitlm( frames , lftloc(frames,yy) );
    vFitsL(yy) = LinFitsL.Coefficients.Estimate(2);
end

%% Grab the insane velocities and reset the points
for fr = frames
    
    %3rd order, for curvy data
    %BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),3);
    %BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),3);
    %predictR(fr,:) = BleachSlopeR(4) + BleachSlopeR(3)*yIntList + BleachSlopeR(2)*yIntList.^2 + BleachSlopeR(1)*yIntList.^3;
    %predictL(fr,:) = BleachSlopeL(4) + BleachSlopeL(3)*yIntList + BleachSlopeL(2)*yIntList.^2 + BleachSlopeL(1)*yIntList.^3;

    %1st order, for straight data
    BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),1);
    BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),1);
    predictR(fr,:) = BleachSlopeR(2) + BleachSlopeR(1)*yIntList;
    predictL(fr,:) = BleachSlopeL(2) + BleachSlopeL(1)*yIntList;

    %0th order
%     BleachSlopeL = polyfit( ylist(~isnan(lftloc(fr,:))),lftloc(fr,~isnan(lftloc(fr,:))),0);
%     BleachSlopeR = polyfit( ylist(~isnan(rgtloc(fr,:))),rgtloc(fr,~isnan(rgtloc(fr,:))),0);
%     predictR(fr,:) = BleachSlopeR(1)*yIntList./yIntList;
%     predictL(fr,:) = BleachSlopeL(1)*yIntList./yIntList;
%     
    correctListL = unique([yIntList( abs(vFitsL) > 5*mean(abs(vFitsL)))  yIntList(isnan(lftloc(fr,:)) )]);
    correctListR = unique([yIntList( abs(vFitsL) > 5*mean(abs(vFitsL)))  yIntList(isnan(lftloc(fr,:)) )]);
    
    for yy = correctListL 
        lftloc(fr,yy) = predictL(fr,yy);
        lftwid(fr,yy) = avgwid;
        lftloc(fr,yy)
        lftamp(fr,yy) = data1D(ceil(lftloc(fr,yy)),yy,fr);
    end
    for yy = correctListR
        rgtloc(fr,yy) = predictR(fr,yy);
        rgtwid(fr,yy) = avgwid;
        rgtamp(fr,yy) = data1D(floor(rgtloc(fr,yy)),yy,fr);
    end
   
    %Smooth these initial position guesses via a hample filter
    
    
   % rgtloc(fr,:) = smooth(hampel(rgtloc(fr,:)));
 %   lftloc(fr,:) = smooth(hampel(lftloc(fr,:)));
    
end




%% Possible Graphing Step
if graph == 1
    
    for fr = frames(end):-1:frames(1)
        figure('Name',['Frame: ' num2str(fr)],'NumberTitle','off');
        imagesc(TPH(:,:,fr))
        colormap('gray')
        hold on;
        plot(rgtloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0)+yBin/2,'r.','MarkerSize',15)
        plot(lftloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0)+yBin/2,'b.','MarkerSize',15)
        plot(predictR(fr,:),ylist,'g-','MarkerSize',15)
        plot(predictL(fr,:),ylist,'g-','MarkerSize',15)
        
        title(['Frame ' num2str(fr) ' Fit'],'FontSize',18,'interpreter','latex');
        ylabel('Bleach Axis [um]','FontSize',18,'interpreter','latex');
        xlabel('Ordered Axis[um]','FontSize',18,'interpreter','latex');
        set(gca,'fontsize',18)
        
    end
end


disp([ 'Using ' num2str(yBotEnd-sum(SkipList)) ' out of '  num2str(yBotEnd) ' rows'])



end