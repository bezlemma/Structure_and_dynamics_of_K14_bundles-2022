function [BAH_list_sum, BAH_list_sum_err] = func_int_fit_list( frames, pickedframe, lftloc,rgtloc,lftamp, rgtamp,lftamp_err, rgtamp_err,lftwid,rgtwid,SkipList,graph)
global yBin
global yBotEnd

%This function simply creates a list of polariteies based off
%comparing amplitudes of bleach lines


%Uses R-squared > 0.50 as a threshold for keeping a linear fit.
%See https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
%for details
yIntList = 1:(yBotEnd-yBin) ;

Lmeany = nanmean(lftwid(:));
Rmeany = nanmean(rgtwid(:));
Lstd = nanstd(lftwid(:)) / 2;
Rstd = nanstd(rgtwid(:)) / 2;
lftwid(lftwid > (Lmeany + Lstd) ) = Lmeany;
rgtwid(rgtwid > (Rmeany + Rstd) ) = Rmeany;
lftwid(lftwid < (Lmeany - Lstd) ) = Lmeany;
rgtwid(rgtwid < (Rmeany - Rstd) ) = Rmeany;

intLamp = lftamp .* lftwid;
intRamp = rgtamp .* rgtwid;

for yy = 1:length(intLamp(1,:))

    fr = pickedframe;
    maxxx =  max(intLamp(fr,yy),intRamp(fr,yy)) + intRamp(fr,yy) + intLamp(fr,yy) -  intRamp(fr,yy) - intLamp(fr,yy); 
    minnn =  min(intLamp(fr,yy),intRamp(fr,yy)) + intRamp(fr,yy) + intLamp(fr,yy) -  intRamp(fr,yy) - intLamp(fr,yy); 
    BAH_list_sum(yy) = (maxxx - minnn) / (maxxx + minnn);
    
    %now do all time points to see the errors
    for fr = 1:length(frames)
    maxxx =  max(intLamp(fr,yy),intRamp(fr,yy)) + intRamp(fr,yy) + intLamp(fr,yy) -  intRamp(fr,yy) - intLamp(fr,yy); 
    minnn =  min(intLamp(fr,yy),intRamp(fr,yy)) + intRamp(fr,yy) + intLamp(fr,yy) -  intRamp(fr,yy) - intLamp(fr,yy); 
    BAH_list_sum_frames(fr,yy) = (maxxx - minnn) / (maxxx + minnn);
    end
    
  BAH_list_sum(BAH_list_sum == 0) = NaN;
    BAH_list_sum(BAH_list_sum == 1) = NaN;  
  
    %now fit and see what that looks like
    LinFits = fitlm(   frames'  , BAH_list_sum_frames(frames,yy) );  
    
    %Only save the data if the fit was good
  
    BAH_list_sum_err(yy) = LinFits.RMSE;%LinFitsR.Coefficients.SE(2);
    
end
  
  
  
  
%

end