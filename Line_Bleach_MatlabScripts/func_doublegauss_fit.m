function [lftloc,rgtloc,lftamp,rgtamp,lftwid,rgtwid,lftamp_err,rgtamp_err] = ...
    func_doublegauss_fit(slice1D,idxValid,lftamp, ...
                         rgtamp,lftloc,rgtloc,lftwid,rgtwid,...
                         minlftloc,minrgtloc,maxlftloc,maxrgtloc,graphgraph)

global xRgt

bckgrndNoise = 5;
trimStart = 1;
trimEnd = xRgt;

xList = (trimStart:trimEnd)';

minExpectedWidth =  5; %The min half-width we expect for a bleach line
maxExpectedWidth = 50; %The max half-width we expect for a bleach line

%fitting parameters go in alphabetical order so
% a1 a2 b1 b2 c1 c2 d
f = fit(xList(idxValid),slice1D(idxValid), ...
    'a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + d1 + e1*x',...
    'start',[lftamp,rgtamp,lftloc,rgtloc,lftwid,rgtwid,bckgrndNoise,0], ...
    'Lower',[lftamp/4,rgtamp/4, ...
    minlftloc,minrgtloc, ...
    minExpectedWidth,minExpectedWidth, -lftamp,-lftamp], ...
    'Upper',[lftamp*4,rgtamp*4,...
    maxlftloc,maxrgtloc,...
      maxExpectedWidth, maxExpectedWidth, lftamp,lftamp] );

%Save the data
lftamp = f.a1; rgtamp = f.a2;
lftloc = f.b1; rgtloc = f.b2;
lftwid = f.c1; rgtwid = f.c2;
ConfInt_95 = confint(f,0.95);

%RMSE = SE = (upper limit – lower limit) / 3.92.
stupidRMSE_left = (   ConfInt_95(2,1) -ConfInt_95(1,1)   )  / 3.92;
stupidRMSE_right = (   ConfInt_95(2,2) -ConfInt_95(1,2)   )  / 3.92;
lftamp_err = stupidRMSE_left;
rgtamp_err = stupidRMSE_right;
lin = f.d1; slope = f.e1;
%conf(:,:) = confint(f)
 
% %Plotting Data to Sanity Check while Debugging
if graphgraph
figure;
plot((1:length(slice1D)),slice1D,'LineWidth',2);
hold on;
plot(xList, lftamp.*exp(-((xList-lftloc)./lftwid).^2) + rgtamp.*exp(-((xList-rgtloc)./rgtwid).^2) + lin + slope.*xList  );
%title(['Y-Pixel: ' num2str ' Fit']);
ylabel('Intensity');
xlabel('Position [um]');
legend('Data','Fit');
% 
end
end