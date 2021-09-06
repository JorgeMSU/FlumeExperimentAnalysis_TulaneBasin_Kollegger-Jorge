%% This script generates Figure 4A from Kollegger et al. 
clc
close all
load('C:\Users\lorenzotruej\Dropbox\MergedFoulder_RESEARCH\MontclairStateUniversity\Research\Madeline\Manuscript drafts\Code\LMLP.mat');
cx=80;%Inlet x coordinate
cy=80;%Inlet y coordinate

LMLP(isnan(LMLP))=0;%make the NaN values zero
% e=612;
e=1340/5;
w=100;%radius for 0.5m
w2=152;%0.76m
w3=244;%1.2m

m=491;
%% Time vector
dt=1;
t=50:dt:540; %hours
tSL=1:dt:m; %sea level increases with time
nt=length(t);
%% Distance
dx=5;%length of each cell is 5mm
Dmax=e*dx;%total length in mm is number of cells * length of each cell
D=dx:dx:Dmax; %cells*cell width % start at 5mm and end at the total length in mm

%%Sea Level
Z=0.25*tSL;%increases 0.25mm/hr
Al=3.0625;%mm
Bl=0.0643;
Zcp=zeros(1,nt);
 
%% Elevation Vectors
heightsl=zeros(1,m); %0.5
heights2l=zeros(1,m);%0.76
heights3l=zeros(1,m);%1.2

residuall=zeros(1,e-w+1);

one=1;

for q=1:m
q
Vl=LMLP(:,:,q);
[a,b] = size(Vl);
[X, Y] = meshgrid( (1:a)-cx, (1:b)-cy);
R = sqrt(X.^2 + Y.^2);%how to define radius

Zr=37.5+Al*sin(q*Bl);
Zrstore(q)=Zr;
Zc=Zr + Z(q);

Zcp(q)=Zc;

     for i = w % radius of the circle
        mask = (w-one<R & R<w+one); % smooth 1 px around the radius
        values = Vl(mask); % without smooth
        total=sum(values); %add all the heights at the distance
        N = nnz(values); %returns the number of nonzero elements in  values
        average=total/N; %calculate average height at this timestep
     end
heightsl(q)=average;%store average height at each timestep

    for i = w2
        mask2 = (w2-one<R & R<w2+one); % smooth 1 px around the radius
        values2 = Vl(mask2); % without smooth
        total2=sum(values2); %add all the heights at the distance
        N2 = nnz(values2); %returns the number of nonzero elements in  values
        average2=total2/N2; %calculate average height at this timestep
     end
heights2l(q)=average2;%store average height at each timestep

    for i = w3
        mask3l = (w3-one<R & R<w3+one); % smooth 1 px around the radius
        values3l = Vl(mask3l); % without smooth
        total3l=sum(values3l); %add all the heights at the distance
        N3l = nnz(values3l); %returns the number of nonzero elements in  values
        average3l=total3l/N3l; %calculate average height at this timestep
     end
heights3l(q)=average3l;%store average height at each timestep
 for i= 1:1:e %every point from the first elevation location to the end
        maskLAGl = (i-1<R & R<i+1); % smooth 1 px around the distance
        valuesLAGl = Vl(maskLAGl); % without smooth
        totalLAGl=sum(valuesLAGl); %add all the heights at the distance
        NLAGl = nnz(valuesLAGl); %returns the number of nonzero elements in  values
        averageLAGl=totalLAGl/NLAGl; %calculate average height at this timestep
        residuall=averageLAGl(:)-Z(q);%residual is average height at the distance -SLR
        residuall(isnan(residuall))=0;
        RESIDUALLMLP(i)=residuall; %store residual
 end
RESIDUALALLLMLP(:,q)=RESIDUALLMLP(:);%save all residuals for i= w:1:e
end

residual1=heightsl-Z-heightsl(1);
residual2=heights2l-Z-heights2l(1);
residual3=heights3l-Z-heights3l(1);
residualSL=Zrstore-37.5;% sealevel oscillation without rise and minus initial elevation


%% Figure 4Apart2 (B)
% %% Residual change over time test
% % upstream
% resrise1=(residual1(24)-residual1(1))/49;
% resfall1=(residual1(74)-residual1(24))/49;
% resrise2=(residual1(122)-residual1(74))/49;
% resfall2=(residual1(172)-residual1(122))/49;
% resrise3=(residual1(220)-residual1(172))/49;
% resfall3=(residual1(270)-residual1(120))/49;
% resrise4=(residual1(318)-residual1(270))/49;
% resfall4=(residual1(368)-residual1(318))/49;
% resrise5=(residual1(414)-residual1(368))/49;
% resfall5=(residual1(465)-residual1(414))/49;
% resrise6=(residual1(491)-residual1(465))/49;
% 
% reschangeup=[resrise1,resfall1,resrise2,resfall2,resrise3,resfall3,resrise4,resfall4,resrise5,resfall5,resrise6];
% xax=[1,2,3,4,5,6,7,8,9,10,11];
% % subplot(2,1,1)
% scatter(xax,reschangeup,'k')
% hold on
% ylim([-0.2 0.15])
% title(['Residual Change'])
% ylabel('Residual Elevation Change')
% xlabel ('Time(Hour)')
% 
% %downstream
% resrise1dwn=(residual3(24)-residual3(1))/49;
% resfall1dwn=(residual3(74)-residual3(24))/49;
% resrise2dwn=(residual3(122)-residual3(74))/49;
% resfall2dwn=(residual3(172)-residual3(122))/49;
% resrise3dwn=(residual3(220)-residual3(172))/49;
% resfall3dwn=(residual3(270)-residual3(120))/49;
% resrise4dwn=(residual3(318)-residual3(270))/49;
% resfall4dwn=(residual3(368)-residual3(318))/49;
% resrise5dwn=(residual3(414)-residual3(368))/49;
% resfall5dwn=(residual3(465)-residual3(414))/49;
% resrise6dwn=(residual3(491)-residual3(465))/49;
% 
% reschangedown=[resrise1dwn,resfall1dwn,resrise2dwn,resfall2dwn,resrise3dwn,resfall3dwn,resrise4dwn,resfall4dwn,resrise5dwn,resfall5dwn,resrise6dwn];
% xax=[1,2,3,4,5,6,7,8,9,10,11];
% % subplot(2,1,2)
% scatter(xax,reschangedown,'m')

%% Figure 4A
%% Residual elevation changes
subplot (2,1,1)
plot(t,residual1,'k')
hold on
plot(t,residual2,'r')
hold on
plot(t,residual3,'m')
xlabel('Time(Hr)')
xlim([50 540])
ylabel('Elevation(mm)');
legend('0.5m','0.76m','1.22m');%'Sea Level Cycles');
title('LMLP Residuals')

subplot (2,1,2)
plot(t,residualSL,'c')
xlabel('Time(Hr)')
xlim([50 540])
ylabel('Elevation(mm)');
legend('Sea Level Cycles');
title('Residual Sea Level')

%% test Timelag for one location
% x=RESIDUALALLLMLP(w,:)-RESIDUALALLLMLP(w,1); %upstream
% y=Zrstore-Zrstore(1);%sea level
%  [xc,lags] = xcorr(x,y); % xcorr returns 2*length of the longest vector-1 lags by default
%  xc = xc(length(x):length(x)+98);
%  lags = lags(length(x):length(x)+98);
%  [val,index] = max(xc);
% % valstore=xc(index);
% timelagupstream=lags(index);
% 
% %% Timelag Analysis
% % Timelag across the fluvial surface
% for l=1:e
%  x = RESIDUALALLLMLP(l,:)-RESIDUALALLLMLP(l,1);
%  y = Zrstore-Zrstore(1);
%  [xc,lags] = xcorr(x,y); % xcorr returns 2*length of the longest vector-1 lags by default
%  xc = xc(length(x):length(x)+98);
%  lags = lags(length(x):length(x)+98);
%  [val,index] = max(xc);
% valstore(l)=xc(index);
% timelag(l)=lags(index);
% end
% 
% Dnew=D(1:e);
% scatter(Dnew,timelag,5, 'k','filled')
% % ylim([0 13])
% xlim([100 1340])
% title(['Timelag Across the Fluvial Surface'])
% xlabel('Distance (mm)')
% ylabel('Timelag (hr)')
% box 'on'

% Moving Average
% highstand=zeros(e,m);
% lowstand=zeros(e,m);
% for mm=1:e
% movemeanmatrix(mm,:)= movmean(RESIDUALALLLMLP(mm,:),[10 10]);
% end
% 
% for em=1:e
% for q=1:m
% if RESIDUALALLLMLP(em,q)>movemeanmatrix(em,q)
%     hs=RESIDUALALLLMLP(em,q);
%     ls=0;
% elseif RESIDUALALLLMLP(em,q)<=movemeanmatrix(em,q)
%     hs=0;
%     ls=RESIDUALALLLMLP(em,q);
% end
% highstand(em,q)=hs(:);
% lowstand(em,q)=ls(:);
% end
% end
% Zrpos=Zrstore-Zrstore(1);
% Zrpos(Zrpos<0)=0;
% 
% Zrneg=Zrstore-Zrstore(1);
% Zrneg(Zrneg>0)=0;
% 
% for l=1:e
%  xpos = highstand(l,:);
%  ypos = Zrpos-Zrpos(1);
%  [xc,lags] = xcorr(xpos,ypos); % xcorr returns 2*length of the longest vector-1 lags by default
%  xc = xc(length(xpos):length(xpos)+24);
%  lags = lags(length(xpos):length(xpos)+24);
%  [val,index] = max(xc);
% timelaghighstand(l)=lags(index);
% end
% 
% for l=1:e
%  xpos = lowstand(l,:);
%  ypos = Zrneg-Zrneg(1);
%  [xc,lags] = xcorr(xpos,ypos); % xcorr returns 2*length of the longest vector-1 lags by default
%  xc = xc(length(xpos):length(xpos)+24);
%  lags = lags(length(xpos):length(xpos)+24);
%  [val,index] = max(xc);
% timelaglowstand(l)=lags(index);
% end

%%  amplitude...dampening
% fs = 1;% sample frequency
% n = length(residual1);
% fbins = [(0:1/n:1-1/n)*fs];    %frequency bin vector for plotting (x axis)
% calval = n/2;                  %for two-sided ffts, calval should just be N for one-sided
% 
% [fftdatres1] = fft(residual1);
% fftmagres1 = abs(fftdatres1)/calval;
% 
% [fftdatres2] = fft(residual2);
% fftmagres2 = abs(fftdatres2)/calval;
% 
% [fftdatres3] = fft(residual3);
% fftmagres3 = abs(fftdatres3)/calval;
% 
% [fftdatSL] = fft(residualSL);
% fftmagSL = abs(fftdatSL)/calval;
% 
% plot(fbins, fftmagSL,'c')%to visualize expected magnitude
% hold on
% plot(fbins, fftmagres1,'k')%to visualize expected magnitude
% hold on
% plot(fbins, fftmagres2,'r')%to visualize expected magnitude
% hold on
% plot(fbins, fftmagres3,'m')%to visualize expected magnitude
% xlabel('frequency')
% ylabel('magnitide (mm)')
% xlim([0.005 0.02])
% ylim([0 4])
% title('LMLP Amplitudes')
% legend('Sea Level','Upstream','Middle', 'Downstream');
% 
% subax = axes('Position',[0.475 0.5 0.4 0.4]); %left bottom width height
% plot(subax, fbins, fftmagSL,'c');
% hold on
% plot(subax, fbins, fftmagres1,'k');
% hold on
% plot(subax, fbins, fftmagres2,'r');
% hold on
% plot(subax, fbins, fftmagres3,'m');
% xlim([0.005 0.015])
% ylim([0 2])
% legend('Sea Level', 'Upstream','Middle','Downstream');
