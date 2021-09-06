%% This script generates Figures 4B, 4C and 4D from Kollegger et al. 
clc
close all
load('C:\Users\lorenzotruej\Dropbox\MergedFoulder_RESEARCH\MontclairStateUniversity\Research\Madeline\Manuscript drafts\Code\LMLP.mat')
cx=80;%Inlet x coordinate 
cy=80;%Inlet y coordinate

LMLP(isnan(LMLP))=0;%make the NaN values zero

w=612; %number of cells in the profile vector (long enough to include the shorline)
m=491; %number of hours of the experiment

%% Time vector
dt=1;%timestep
tmax=540;%hour 1170 is the end of the experiment
t=50:dt:tmax; %hours
tSL=1:dt:m; %sea level increases with time
nt=length(t);%number of hours

%% Sea Level
Z=0.25*tSL;%Background sea-level rise rate of 0.25mm/hr
Amp=3.0625;% LMLP - Amplitude of sea-level oscillations in "mm"
B=0.0643;%LMLP Frequency of sea-level oscillations = 2*pi/period (Period is 98 hours)
Zcp=zeros(1,nt);%sea level oscillations from initial height at current hour with backround rise
Zr=zeros(1,nt); %sea level residuals (after substracting backround sea-level rise)
ZcP=zeros(w,m);%matrix of all Zcp

%% Shoreline
sc=zeros(1,m);%x coordinate for shoreline
YB=zeros(1,m); %y coordinate for shoreline
cut=20;%20=100mm 10=50mm % 

%% Distance
dx=5;%length of each cell is 5mm
Dmax=w*dx;%total length in mm is number of cells * length of each cell
D=dx:dx:Dmax; %cells*cell width % start at 5mm and end at the total length in mm

%% Elevation Vectors
heights=zeros(1,m); %0.5m radius

%% Profile
PROFILEALL = zeros(w,m); %matrix with each average profile for each hour
PROFILE=zeros(1,w);%average profile of current hour

%% Volume
erosion=zeros(1,m);
deposition=zeros(1,m);

for q=1:m %%Do this for every hour in the 3D matrix
q
V=LMLP(:,:,q);%pulls out the matrix (from Anew) for the current hour
[a,b] = size(V);%gets the size of the matrix
[X, Y] = meshgrid( (1:a)-cx, (1:b)-cy);
R = sqrt(X.^2 + Y.^2);%how to define radius %distance/width of cell
Zr(q)=Amp*sin(q*B);
Zc=37.5+ Zr(q)  + Z(q); %sea level starts at 37.5 in the LMLP and oscillates ************* Need to justify 37.5 *************
Zcp(q)=Zc; %store sea level at each hour
ZCP=Zcp(q)*ones(1,length(D));%horizontal sea level
ZcP(:,q)=ZCP; %saves the horizotal sea level at each hour

    %%calculate average at each distance for ONE hour
    for i = 1:1:w %do this for each distance until w
        mask = (i-1<R & R<i+1); % smooth 1 px around the radius
        values = V(mask); % without smooth
        total=sum(values); %add all the heights at the distance
        N = nnz(values); %returns the number of nonzero elements in  values
        average=total/N; %calculate average height at this timestep
        profile = average(:);%average profile at that distance
        PROFILE(i)=profile; %saves all the averages of distances together
    end

    PROFILEALL(:,q)=PROFILE(:); %save for every hour
    
    %%Find shoreline
     [x,y]=curveintersect(D,PROFILEALL(:,q),D,ZcP(:,q));% wherever the profile intersects with the sea level is the shoreline
      if x>0
      sc(q)=max(x);%sc is shoreline location
      YB(q)=max(y);
      end
      %

%       %% DERIVATIVES
   scrc(q)=round(sc(q)/dx); %rounds the location of the shoreline
   %FIRST DERIVATIVE
   AP1(q)=sum(diff(PROFILE((cut:1:scrc(q)))))/sc(q);% calculates the average first derivative of each profile between cut and the shoreline

   %AREA UNDER THE PROFILE CURVE
    areareal(q)=trapz(PROFILE((cut:1:(scrc(q)))))*dx;%calculates area under the real profile

    %AREA UNDER THE IDEALIZED CURVE
     Xideal(:,q)=[D(cut) D(cut) D(scrc(q)) D(scrc(q))]; %x values for each idealized curve
     Yideal(:,q)=[0 PROFILE(cut) PROFILE(scrc(q)) 0];%y values for each idealized curve
     polygonideal(:,q)=polyshape(Xideal(:,q),Yideal(:,q)); %makes the polygon of each idealized curve
     areapoly(q)=area(polygonideal(:,q));%calculates the area of the polygon of each idealized curve
      
     %Difference in area % +/- indicates shape of fluvial surface
     difference(q)=(areapoly(q)-areareal(q))/sc(q);% takes the difference in area between the idealized and actual profile
      
%% Volume Eroded (Not used in Kolleger et al.)      
% if q>1
%      if sc(q)<sc(q-1)%if shoreline is retreating% so we dont add the box on the end
%        previousprofile=PROFILEALL((cut:1:scrc(q)),q-1);
%        nowprofile=PROFILEALL(cut:1:scrc(q),q);
%        previous=trapz(previousprofile)*dx;%area of previous timestep
%        now=trapz(nowprofile)*dx; %area of current timestep
%        volumechange=((now-previous))/sc(q);%change in volume is the difference in the area of previousand current time step - the backround erosion
%          
%      elseif sc(q)>=sc(q-1)%shoreline is prograding
%        previousprofile=PROFILEALL(cut:1:scrc(q-1),q-1);
%        nowprofile=PROFILEALL(cut:1:scrc(q-1),q);
%        previous=areareal(q-1);%area of previous timestep
%        now=trapz(nowprofile)*dx;%area of current timestep
%        volumechange=((now-previous))/sc(q); %change in volume is the difference in the area of previousand current time step - the backround erosion        
%      end
%% Volume Change Location (Not used in Kolleger et al.) 
%     newXp=0:.01:1; 
%     xaxis=cut:1:scrc(q);
%      profiledifference=(PROFILEALL(xaxis,q)-PROFILEALL(xaxis,q-1));%area difference beteen current and previous profile
%      d=profiledifference;%setting deposition to profile difference
%      d(d<0)=0;%remove negative values
%      e=profiledifference; %set deposition to profile difference
%      e(e>=0)=0;%removing positive values
%      
%      if Zr(q)>Zr(q-1)&& Zr(q)>0 %for sealevel rise PHASE 1
%          Drise=d;%find deposition for the current hour
%          newX=linspace(0,1,length(Drise));
%          Erise=e;%find erosion for the current hour
%          newDrise(:,q) = interp1(newX, Drise, newXp);
%          newErise(:,q) = interp1(newX, Erise, newXp);
%          phase1d(:,q)=newDrise(:,q);
%          phase1e(:,q)=newErise(:,q);
% %          DR(:,j)=Dr;%save depostion location for each hour 
% %          ER(:,j)=Er;%save erosion location for each hour 
%      elseif Zr(q)>Zr(q-1)&& Zr(q)<=0 %for sealevel rise PHASE 4
%          Drise=d;%find deposition for the current hour
%          newX=linspace(0,1,length(Drise));
%          Erise=e;%find erosion for the current hour
%          newDrise(:,q) = interp1(newX, Drise, newXp);
%          newErise(:,q) = interp1(newX, Erise, newXp);
%          phase4d(:,q)=newDrise(:,q);
%          phase4e(:,q)=newErise(:,q);
%      elseif Zr(q)<=Zr(q-1)&& Zr(q)>0 %for sea level fall PHASE 2
%          Dfall=d;%find deposition for the current hour
%          newX=linspace(0,1,length(Dfall));
%          Efall=e;%find erosion for the current hour
%          newDfall(:,q) = interp1(newX, Dfall, newXp);
%          newEfall(:,q) = interp1(newX, Efall, newXp);
%          phase2d(:,q)=newDfall(:,q);
%          phase2e(:,q)=newEfall(:,q);
% %          DF(:,j)=Df;%save depostion location for each hour 
% %          EF(:,j)=Ef;%save erosion location for each hour 
%      elseif Zr(q)<=Zr(q-1)&& Zr(q)<=0 %for sea level fall PHASE 3
%          Dfall=d;%find deposition for the current hour
%          newX=linspace(0,1,length(Dfall));
%          Efall=e;%find erosion for the current hour
%          newDfall(:,q) = interp1(newX, Dfall, newXp);
%          newEfall(:,q) = interp1(newX, Efall, newXp);
%          phase3d(:,q)=newDfall(:,q);
%          phase3e(:,q)=newEfall(:,q);    
%      end
   
%% Volume Change (Not used in Kollegger et al.)    
%      erosion=profiledifference(profiledifference<0)/sc(q);% takes the negative values from profiledifference and saves as erosion
%      deposition=profiledifference(profiledifference>=0)/sc(q);  %takes the positive values from profiledifference and saves as deposition
%      EROSIONstore(q)=sum(erosion)*dx;%same as volume change erode
%      DEPOSITIONstore(q)=sum(deposition)*dx;%same as volume change deposit
%      VOLUMECHANGE(q)=volumechange;%save volume change at each hour
% end
end
 averageslope=mean(AP1);
%% Figure 4 C and D
%% Derivatives
figure
subplot(3,1,1) %% Figure 4c
plot(t,AP1,'k')
xlim([50 540])
xlabel('Time(Hour)')
title(['First Order Derivative'])

subplot(3,1,2) %% Figure 4d
plot(t,difference,'b')
hold on
% yline(averageCurvature)
% plot(t,RELIEF,'b')
xlim([50 540])
xlabel('Time(Hour)')
% ylabel('mm')
title(['Difference'])

subplot(3,1,3)
plot(t,Zcp-Z,'c')
xlim([50 540])
xlabel('Time(Hour)')
ylabel('Elevation(mm)')
title(['Sea Level'])

%% Figure 4 B
%% LMLP derivative average slope per SL phase 
rise1=(AP1(24)-AP1(1))/49;
fall1=(AP1(74)-AP1(24))/49;
rise2=(AP1(122)-AP1(74))/49;
fall2=(AP1(172)-AP1(122))/49;
rise3=(AP1(220)-AP1(172))/49;
fall3=(AP1(270)-AP1(120))/49;
rise4=(AP1(318)-AP1(270))/49;
fall4=(AP1(368)-AP1(318))/49;
rise5=(AP1(414)-AP1(368))/49;
fall5=(AP1(465)-AP1(414))/49;
rise6=(AP1(491)-AP1(465))/49;
AVERAGEslope=(AP1(491)-AP1(1))/49;

AAslope=[rise1,fall1,rise2,fall2,rise3,fall3,rise4,fall4,rise5,fall5,rise6];
xax=[1,2,3,4,5,6,7,8,9,10,11];
scatter(xax,AAslope)

yline(AVERAGEslope)
title(['Derivative of the Average Slope'])
ylabel('Average Slope')
xlabel ('Time(hr)')



%% Figure 4 F (Not used in Kolleger et al.)
%% change of the curvature (Not used in Kolleger et al.)
% Drise1=(difference(24)-difference(1))/49;
% Dfall1=(difference(74)-difference(24))/49;
% Drise2=(difference(122)-difference(74))/49;
% Dfall2=(difference(172)-difference(122))/49;
% Drise3=(difference(220)-difference(172))/49;
% Dfall3=(difference(270)-difference(120))/49;
% Drise4=(difference(318)-difference(270))/49;
% Dfall4=(difference(368)-difference(318))/49;
% Drise5=(difference(414)-difference(368))/49;
% Dfall5=(difference(465)-difference(414))/49;
% Drise6=(difference(491)-difference(465))/49;
% 
% Dcurvature=[Drise1,Dfall1,Drise2,Dfall2,Drise3,Dfall3,Drise4,Dfall4,Drise5,Dfall5,Drise6];
% xax=[1,2,3,4,5,6,7,8,9,10,11];
% scatter(xax,Dcurvature)
% title(['Derivative of the Curvature'])
% ylabel('Curvature')
% xlabel ('Time(hr)')

%% Appendix (Not used in Kolleger et al.)
%% Phases (Not used in Kolleger et al.)
% phase1=(sum(phase1d,2))+(sum(phase1e,2));
% phase2=(sum(phase2d,2))+(sum(phase2e,2));
% phase3=(sum(phase3d,2))+(sum(phase3e,2));
% phase4=(sum(phase4d,2))+(sum(phase4e,2));

%% Sedimentation (Not used in Kolleger et al.)
% subplot(2,1,1) 
% plot(newXp, phase1+phase4,'r')
% title(['Net Sedimentation during Sea Level Rise'])
% ylabel('mm')
% 
% subplot(2,1,2)
% plot(newXp, phase2+phase3,'b')
% title(['Net Sedimentation during Sea Level Fall'])
% ylabel('mm')
