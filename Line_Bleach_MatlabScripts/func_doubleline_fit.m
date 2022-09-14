function [lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,lftamp_err,rgtamp_err,SkipList] = ...
    func_doubleline_fit(data1D,frames,TPH,graph,SkipList,lftloc,rgtloc,...
    lftamp,rgtamp,lftwid,rgtwid)

global yBin
global yBotEnd
global xRgt

%% ---------------- Fitting parameters
trimStart = 1; %There is off-by-1 mistmatch here between trim and idxValid!
trimEnd = xRgt;
ylist = (1/2):(yBotEnd-yBin);
yIntList =  1:(yBotEnd-yBin);
NARROW = 20; %this can be 3.5 at smallest, need as many pixels as we have free variables to fit
graphgraph = 0; %if you want to look at individual fits while debugging, set this to 1
     
%% ---------------- Fit gaussians to the following frames
for fr = frames(end):-1:frames(1)
    disp(['Now analayzing frame: ' num2str(fr)] )
    
    for yy = yIntList
        
        if SkipList(yy,fr) == 1 
            continue
        end
        
        %Narrow the window to NARROW pixels centered around each
        %gaussian (a minimum of 8 points is needed to fit, we have NARROW*2*2)
        removeList2 = NaN*zeros(length(1:trimEnd),1);
        for ii = (trimStart+3):(trimEnd-3)
            if  (ii < (lftloc(fr,yy)+ NARROW)) &&  (ii > (lftloc(fr,yy) - NARROW))
                removeList2(ii) = 1;
            end
            if  (ii < (rgtloc(fr,yy)+ NARROW)) &&  (ii > (rgtloc(fr,yy) - NARROW))
                removeList2(ii) = 1;
            end
        end
        
        idxValid2 = ~isnan(removeList2);
            
            %By default say the errors are 0, they will get real errors if
            %a real fit occurs
            lftamp_err(fr,yy) = 0;
            rgtamp_err(fr,yy) = 0;
            
        if ~isnan(lftloc(fr,yy)) && ~isnan(rgtloc(fr,yy))
                [lftloc(fr,yy), rgtloc(fr,yy),lftamp(fr,yy), rgtamp(fr,yy),lftwid(fr,yy),rgtwid(fr,yy),~,~] = ...
            func_doublegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  lftamp(fr,yy),rgtamp(fr,yy),...
            lftloc(fr,yy),rgtloc(fr,yy), ...
            lftwid(fr,yy),rgtwid(fr,yy),...  %Guesses at Widths 
            (lftloc(fr,yy)-NARROW/2),(rgtloc(fr,yy)-NARROW/2),... %Min Locations
            (lftloc(fr,yy)+NARROW/2),(rgtloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        elseif ~isnan(lftloc(fr,yy))
             [lftloc(fr,yy),lftamp(fr,yy),lftwid(fr,yy)] = ...
            func_singlegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  lftamp(fr,yy),lftloc(fr,yy),...
            lftwid(fr,yy),(lftloc(fr,yy)-NARROW/2),... 
            (lftloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        elseif ~isnan(rgtloc(fr,yy))
              [rgtloc(fr,yy),rgtamp(fr,yy),rgtwid(fr,yy)] = ...
            func_singlegauss_fit(data1D(trimStart:trimEnd,yy,fr),...
            idxValid2,  rgtamp(fr,yy),rgtloc(fr,yy),...
            rgtwid(fr,yy),(rgtloc(fr,yy)-NARROW/2),... 
            (rgtloc(fr,yy)+NARROW/2),graphgraph); %Max Locations
        else
            the code should never reach here;
        end             
            
    end   %End Y Loop
end  %End Fr Loop


%% Possible Graphing Step
if graph == 1
    for fr = frames(end):-1:frames(1)
        figure('Name',['Frame: ' num2str(fr)],'NumberTitle','off');
        imagesc(TPH(:,:,fr))
        colormap('gray')
        hold on;
        plot(rgtloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0)+yBin/2,'r.','MarkerSize',15)
        plot(lftloc(fr,SkipList(:,fr) == 0),ylist(SkipList(:,fr) == 0)+yBin/2,'b.','MarkerSize',15)
        title(['Frame ' num2str(fr) ' Fit'],'FontSize',18,'interpreter','latex');
        ylabel('Bleach Axis [um]','FontSize',18,'interpreter','latex');
        xlabel('Ordered Axis[um]','FontSize',18,'interpreter','latex');
        set(gca,'fontsize',18)
        
    end
end

end