%This function creates a list of velocity fits based off the polyfit
%function, while throwing out outlier position data

%This version, compares position to the Center of Mass (CoM) of the whole
%system, by taking the average of all positions. 

function [vFitsL, vFitsR, vErrL, vErrR, CoM] = func_vel_fit_list_CoM( frames, lftloc, rgtloc,SkipList,time_list,graph)
global yBin
global yBotEnd
global TPH
global umPerPixel

yIntList = 1:(yBotEnd-yBin) ;
%Uses R-squared > 0.50 as a threshold for keeping a linear fit.
%See https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
%for details
minFit = 0/100; %(value from 0 to 1 that tells how close a linear fit neads to be to be accepted)

%Calculate Center's of Mass for all time frames
for i = 1:length(frames)
    CoM(i) = ( nanmean(lftloc(i,:)) + nanmean(rgtloc(i,:))   ) / 2 ;
end


%Fit a velocity to the motion
for yy = 1:length(rgtloc(1,:))
    %Remove outliers
    k = yIntList( isoutlier(rgtloc(frames,yy))  );
    rgtloc(k,yy) = NaN;
    k = yIntList( isoutlier(lftloc(frames,yy)) );
    lftloc(k,yy) = NaN;
    
    %Fit
    LinFitsR = fitlm( time_list(frames)' , (rgtloc(frames,yy) - CoM(frames)'     )*umPerPixel );
    LinFitsL = fitlm( time_list(frames)' , (lftloc(frames,yy) - CoM(frames)'     )*umPerPixel );  
    
    %Only save the data if the fit was good
    if  LinFitsR.Rsquared.Adjusted > minFit %&& LinFitsR.Coefficients.Estimate(2) > 0
        vFitsR(yy) = LinFitsR.Coefficients.Estimate(2);
        vErrR(yy) = LinFitsR.RMSE;%LinFitsR.Coefficients.SE(2);
    end
    if LinFitsL.Rsquared.Adjusted > minFit %&& LinFitsL.Coefficients.Estimate(2) < 0
        vFitsL(yy) = LinFitsL.Coefficients.Estimate(2);
        vErrL(yy) = LinFitsL.RMSE;%LinFitsL.Coefficients.SE(2);
    end 
end


%Save the data only if one fit is good, using the mean of the other. 
for yy = 1:length(rgtloc(1,:))
    LinFitsR = fitlm( time_list(frames)'   ,  (rgtloc(frames,yy) - CoM(frames)'     )*umPerPixel );
    LinFitsL = fitlm( time_list(frames)'  ,  (lftloc(frames,yy) - CoM(frames)'     )*umPerPixel );  
    
    if LinFitsR.Rsquared.Adjusted > minFit && LinFitsL.Rsquared.Adjusted > minFit
            %do nothing
    elseif LinFitsR.Rsquared.Adjusted > minFit && LinFitsL.Rsquared.Adjusted < minFit
        vFitsR(yy) = LinFitsR.Coefficients.Estimate(2);
        vFitsL(yy) = nanmean(vFitsL);
        vErrL(yy) = 2*nanmean(vErrL); %this is wrong, but it only happens ~1% of the time
    elseif LinFitsR.Rsquared.Adjusted < minFit && LinFitsL.Rsquared.Adjusted > minFit
        vFitsL(yy) = LinFitsL.Coefficients.Estimate(2);
        vFitsR(yy) = nanmean(vFitsR);
         vErrR(yy) = 2*nanmean(vErrL);
    else
        vFitsL(yy) = NaN;
        vFitsR(yy) = NaN;
    end
    
end


%% Possible Graphing Step to show what lines were kept
if graph == 1
ylist = (1/2):(yBotEnd-yBin);
for fr = frames(end):-1:frames(1)
    figure('Name',['Frame: ' num2str(fr)],'NumberTitle','off');
    imagesc(TPH(:,:,fr))
    colormap('gray')
    hold on;
    
    plot(CoM(fr),length(ylist)/2,'kx','MarkerSize',40)
    
    plot(rgtloc(fr,isnan(vFitsR)),ylist(isnan(vFitsR)),'k.','MarkerSize',15)
    plot(lftloc(fr,isnan(vFitsL)),ylist(isnan(vFitsL)),'k.','MarkerSize',15)
    
    plot(rgtloc(fr,~isnan(vFitsR)),ylist(~isnan(vFitsR)),'r.','MarkerSize',15)
    plot(lftloc(fr,~isnan(vFitsL)),ylist(~isnan(vFitsL)),'b.','MarkerSize',15)
    
    plot( lftloc(fr,~isnan(vFitsL+vFitsR)),  ylist(~isnan(vFitsL+vFitsR)),'g.','MarkerSize',12)
    plot(rgtloc(fr,~isnan(vFitsL+vFitsR)), ylist(~isnan(vFitsL+vFitsR)),'g.','MarkerSize',12)
    
    title(['Frame ' num2str(fr) ' Fit'],'FontSize',18,'interpreter','latex');
    ylabel('Bleach Axis [um]','FontSize',18,'interpreter','latex');
    xlabel('Ordered Axis[um]','FontSize',18,'interpreter','latex');
    legend('Center Of Mass','Right Line','Left Line');
    set(gca,'fontsize',18)
    
end
end
end

