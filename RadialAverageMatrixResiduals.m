%% This script generates Figures 2D and 2E from Kollegger et al. 
clc
close all
load('C:\Users\lorenzotruej\Dropbox\MergedFoulder_RESEARCH\MontclairStateUniversity\Research\Madeline\Manuscript drafts\Code\Anew.mat')
cx=80;%Inlet x coordinate
cy=80;%Inlet y coordinate

Anew(isnan(Anew))=0;%make the NaN values zero
% e=612;
e=1340/5;
w0=50; %radius for 0.25
w=100;%radius for 0.5
w2=200;%1m
w3=244;%1.2m

m=491;
%% Time vector
dt=1;
t=680:dt:1170; %hours
tSL=1:dt:m; %sea level increases with time
nt=length(t);
%% Distance
dx=5;%length of each cell is 5mm
Dmax=612*dx;%total length in mm is number of cells * length of each cell
D=dx:dx:Dmax; %cells*cell width % start at 5mm and end at the total length in mm

resallpos=zeros(e,m);
resallneg=zeros(e,m);

%%Sea Level
Z=0.25*tSL;%increases 0.25mm/hr
Ah=12.25;%MSLR/2
Bh=0.257;%TSLR
Zcp=zeros(1,nt);
Zrstore=zeros(1,m);
 
%% Elevation Vectors
heights0=zeros(1,m); %0.25
heights=zeros(1,m); %0.5
heights2=zeros(1,m);%0.76
heights3=zeros(1,m);%1.2

residual=zeros(1,e-w+1);

one=1;

for q=1:m
q
V=Anew(:,:,q);
[a,b] = size(V);
[X, Y] = meshgrid( (1:a)-cx, (1:b)-cy);
R = sqrt(X.^2 + Y.^2);%how to define radius

Zr=195+Ah*sin(q*Bh);%sea level oscillation without backround rise
Zrstore(q)=Zr; %store Zr

Zc=Zr + Z(q); %sea level oscillations with backround rise
Zcp(q)=Zc; %store Zc
     for i = w0 % radius of the circle
        mask0 = (w0-one<R & R<w0+one); % smooth 1 px around the radius
        values0 = V(mask0); % without smooth
        total0=sum(values0); %add all the heights at the distance
        N0 = nnz(values0); %returns the number of nonzero elements in  values
        average0=total0/N0; %calculate average height at this timestep
     end
heights0(q)=average0;%store average height at each timestep

     for i = w % radius of the circle
        mask = (w-one<R & R<w+one); % smooth 1 px around the radius
        values = V(mask); % without smooth
        total=sum(values); %add all the heights at the distance
        N = nnz(values); %returns the number of nonzero elements in  values
        average=total/N; %calculate average height at this timestep
     end
heights(q)=average;%store average height at each timestep

    for i = w2
        mask2 = (w2-one<R & R<w2+one); % smooth 1 px around the radius
        values2 = V(mask2); % without smooth
        total2=sum(values2); %add all the heights at the distance
        N2 = nnz(values2); %returns the number of nonzero elements in  values
        average2=total2/N2; %calculate average height at this timestep
     end
heights2(q)=average2;%store average height at each timestep

    for i = w3
        mask3 = (w3-one<R & R<w3+one); % smooth 1 px around the radius
        values3 = V(mask3); % without smooth
        total3=sum(values3); %add all the heights at the distance
        N3 = nnz(values3); %returns the number of nonzero elements in  values
        average3=total3/N3; %calculate average height at this timestep
     end
heights3(q)=average3;%store average height at each timestep

 for i= 1:1:e %every point from the first elevation location to the end
        maskLAG = (i-1<R & R<i+1); % smooth 1 px around the distance
        valuesLAG = V(maskLAG); % without smooth
        totalLAG=sum(valuesLAG); %add all the heights at the distance
        NLAG = nnz(valuesLAG); %returns the number of nonzero elements in  values
        averageLAG=totalLAG/NLAG; %calculate average height at this timestep
        residual=averageLAG(:)-Z(q);%residual is average height at the distance -SLR
        residual(isnan(residual))=0;
        RESIDUAL(i)=residual; %store residual
 end
RESIDUALALL(:,q)=RESIDUAL(:);%save all residuals for i= w:1:e
end

residual0=heights0-Z-heights0(1);
residual1=heights-Z-heights(1);
residual2=heights2-Z-heights2(1);
residual3=heights3-Z-heights3(1);
residualSL=Zrstore-195;% sealevel oscillation without rise and minus initial elevation

%% Figure 2D
%% Residuals
figure (1)
subplot (2,1,1)
plot(t,heights-Z-heights(1),'k')
hold on
plot(t,heights2-Z-heights2(1),'r')
hold on
plot(t,heights3-Z-heights3(1),'m')
hold on
yline(0)
xlabel('Time(Hr)')
xlim([680 1170])
ylabel('Elevation(mm)');
legend('0.5m','1m','1.22m');%'Sea Level Cycles');
title('Hour 680-1170 Residuals')

subplot (2,1,2)
plot(t,residualSL,'c')
hold on
yline(0)
xlabel('Time(Hr)')
xlim([680 1170])
ylabel('Elevation(mm)');
legend('Sea Level Cycles');
title('Sea Level Residual')

%% Figure 2E
%% Timelag across the fluvial surface
for l=1:e
 x = RESIDUALALL(l,:)-RESIDUALALL(l,1);
 y = Zrstore-Zrstore(1);
 [xc,lags] = xcorr(x,y); % xcorr returns 2*length of the longest vector-1 lags by default
 xc = xc(length(x):length(x)+24);
 lags = lags(length(x):length(x)+24);
 [val,index] = max(xc);
valstore(l)=xc(index);
timelag(l)=lags(index);
end

figure (2)
Dnew=D(1:268);
scatter(Dnew,timelag,5, 'k','filled')
% plot(Zrstore 'c')
ylim([0 13])
xlim([100 1340])
title(['Timelag Across the Fluvial Surface'])
xlabel('Distance (mm)')
ylabel('Timelag (hr)')
box 'on'

%% Appendix Figures
%% ACTUAL TIMESERIES
% subplot (2,1,1)
% plot(t,heights0-Z-heights0(1),'g')
% hold on
% plot(t,heights,'k')
% hold on
% plot(t,heights2,'r')
% hold on
% plot(t,heights3,'m')
% hold on
% plot(t,Zcp,'c')
% % yline(0)
% xlabel('Time(Hr)')
% xlim([680 1170])
% ylabel('Elevation(mm)');
% legend('0.5m','0.76m','1.22m','Sea Level Cycles');
% title('HMSP Residuals')

% subplot (2,1,2)
% plot(t,residualSL,'c')
% hold on
% yline(0)
% xlabel('Time(Hr)')
% xlim([680 1170])
% ylabel('Elevation(mm)');
% legend('Sea Level Cycles');
% title('Residual Sea Level')

%% Dampening
%%  amplitude...dampening
% fs = 1;% sample frequency
% n = length(residual1h);
% fbins = [(0:1/n:1-1/n)*fs];    %frequency bin vector for plotting (x axis)
% subfbins = [0,0.2];
% calval = n/2;                  %for two-sided ffts, calval should just be N for one-sided
% 
% [fftdatres1] = fft(residual1h);
% fftmagres1 = abs(fftdatres1)/calval;
% 
% [fftdatres2] = fft(residual2h);
% fftmagres2 = abs(fftdatres2)/calval;
% 
% [fftdatres3] = fft(residual3h);
% fftmagres3 = abs(fftdatres3)/calval;
% 
% [fftdatSL] = fft(residualSLh);
% fftmagSL = abs(fftdatSL)/calval;
% semilogx(fbins, fftmagSL,'c')%to visualize expected magnitude
% hold on
% semilogx(fbins, fftmagres1,'k')%to visualize expected magnitude
% hold on
% semilogx(fbins, fftmagres2,'r')%to visualize expected magnitude
% hold on
% semilogx(fbins, fftmagres3,'m')%to visualize expected magnitude
% xlabel('frequency')
% xlim([0.03,0.07])
% ylim([0 13])
% ylabel('magnitide (mm)')
% title('HMSP Amplitudes')
% 
% 
% subax = axes('Position',[0.475 0.5 0.4 0.4]); %left bottom width height
% plot(subax, fbins, fftmagSL,'c');
% hold on
% plot(subax, fbins, fftmagres1,'k');
% hold on
% plot(subax, fbins, fftmagres2,'r');
% hold on
% plot(subax, fbins, fftmagres3,'m');
% xlim([0.035,0.045])
% ylim([0 2])
% legend('Sea Level', 'Upstream','Middle','Downstream');
% plot(fbins, fftmagSL,'c')%to visualize expected magnitude
% hold on
% plot(fbins, fftmagres1,'k')%to visualize expected magnitude
% hold on
% plot(fbins, fftmagres2,'r')%to visualize expected magnitude
% hold on
% plot(fbins, fftmagres3,'m')%to visualize expected magnitude
% xlabel('frequency')
% ylabel('magnitide (mm)')
% title('HMSP Amplitudes')
% 
% subax = axes('Position',[0.3 0.4 0.4 0.4]);
% plot(subax, fbins, fftmagSL,'c');
% hold on
% plot(subax, fbins, fftmagres1,'k');
% hold on
% plot(subax, fbins, fftmagres2,'r');
% hold on
% plot(subax, fbins, fftmagres3,'m');
% xlim([0.02,0.06])
% ylim([0,13])


%% Extra, sedimentation video
% v = VideoWriter('C:\Users\kolleggerm1\Desktop\movies\test.avi');
%    open(v);
% for q = 1:m    
%     if q==1
%     plot(D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q),'w')
%     elseif q>1
%        if Zcp(q)>Zcp(q-1) %sea level rise
%             shade(D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q),D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q-1),'Color','none','FillType',[1,2;2,1],'FillColor', [0 1 0; 1 0 0])
%        elseif Zcp(q)<Zcp(q-1) %sea level fall
%         shade(D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q),D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q-1),'Color','none','FillType',[1,2;2,1],'FillColor', [0 1 0; 1 0 0])      
%        end
%     end
%     hold on
%      
%     ylim([160 350])%change if needed
%     xlim([cut*dx 1340])
%     title(['SL Fall Hour:' num2str(q)])
%     xlabel('Distance (mm)')
%     ylabel('Elevation (mm)')
% %     legend('Deposition', 'Erosion')
% %     hold off    
%     drawnow
%      frame = getframe(gcf);
%     writeVideo(v,frame);
% end
%      close(v);
%    movie2gif(frame(1:end-1), 'C:\Users\kolleggerm1\Desktop\movies\SLF.gif', 'LoopCount', inf, 'DelayTime', 0.05)
% esum=sum(estore,2);
% dsum=sum(dstore,2);
% figure(2)
% plot(D(cut:1:end), esum+dsum,'k')
% hold on
% plot(D(cut:1:end),PROFILEALL(cut:1:end,491),'r')

