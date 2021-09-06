%% This script generates Figures 3C and 3D from Kollegger et al. 
clc
close all
load('C:\Users\lorenzotruej\Dropbox\MergedFoulder_RESEARCH\MontclairStateUniversity\Research\Madeline\Manuscript drafts\Code\Anew.mat')% Use "LMLP" when appropriate
cx=80;%Inlet x coordinate
cy=80;%Inlet y coordinate

Anew(isnan(Anew))=0;%make the NaN values zero

w=612; %length of matrix AKA number of horizontal cells
m=491; %number of hours

%% Time vector
dt=1;%timestep
tmax=1170;%hour 1170 is the end of the experiment
t=680:dt:tmax; %hours
tSL=1:dt:m; %sea level increases with time
nt=length(t);%number of hours

%% Sea Level
Z=0.25*tSL;%increases 0.25mm/hr
Amp=12.25;%MSLR/2
B=0.257;%TSLR
Zcp=zeros(1,nt);%sea level oscillations from initial height at current hour with backround rise
Zr=zeros(1,nt); %sea level oscillations without backround rise
ZcP=zeros(w,m);%matrix of all Zcp

%% Shoreline
sc=zeros(1,m);%x coordinate for shoreline
YB=zeros(1,m); %y coordinate for shoreline
cut=20;%20=100mm 10=50mm %clears up inconsistencies at the beginning of the profile

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
V=Anew(:,:,q);%pulls out the matrix (from Anew) for the current hour
[a,b] = size(V);%gets the size of the matrix
[X, Y] = meshgrid( (1:a)-cx, (1:b)-cy);
R = sqrt(X.^2 + Y.^2);%how to define radius %distance/width of cell
Zr(q)=Amp*sin(q*B);
Zc=195+ Zr(q)  + Z(q); %sea level starts at 195mm and oscillates 
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
        PROFILE(i)= profile; %saves all the averages of distances together
    end

    PROFILEALL(:,q)= PROFILE(:); %save for every hour
    
    %%Find shoreline
     [x,y]=curveintersect(D,PROFILEALL(:,q),D,ZcP(:,q));% wherever the profile intersects with the sea level is the shoreline
      if x>0
      sc(q)=max(x);%sc is shoreline location
      YB(q)=max(y);
      end
      sc(66)=sc(65);%DEM scan only went to 1.3 m and these values were beyond that so guessing value is close to previous value
      sc(67)=sc(65); sc(68)=sc(65);  YB(66)=YB(65);   YB(67)=YB(65);   YB(68)=YB(65);
      
newXp=0:.01:1;

      %% DERIVATIVES
   scrc(q)=round(sc(q)/dx); %rounds the location of the shoreline
   %FIRST DERIVATIVE
   averageprofile1(q)=sum(diff(PROFILE((cut:1:scrc(q)))))/sc(q);% calculates the average first derivative of each profile between cut and the shoreline

   %AREA UNDER THE PROFILE CURVE
    areareal(q)=trapz(PROFILE((cut:1:(scrc(q)))))*dx;%calculates area under the real profile

    %AREA UNDER THE IDEALIZED CURVE
     Xideal(:,q)=[D(cut) D(cut) D(scrc(q)) D(scrc(q))]; %x values for each idealized curve
     Yideal(:,q)=[0 PROFILE(cut) PROFILE(scrc(q)) 0];%y values for each idealized curve
     polygonideal(:,q)=polyshape(Xideal(:,q),Yideal(:,q)); %makes the polygon of each idealized curve
     areapoly(q)=area(polygonideal(:,q));%calculates the area of the polygon of each idealized curve
      
     %Difference in area % +/- indicates shape of fluvial surface
     difference(q)=(areapoly(q)-areareal(q))/sc(q);% takes the difference in area between the idealized and actual profile
      
     
if q>1
    %% Volume Eroded 
%      if sc(q)<sc(q-1)%if shoreline is retreating% so we dont add the box on the end
%        previousprofile=PROFILEALL((cut:1:scrc(q)),q-1);
%        nowprofile=PROFILEALL(cut:1:scrc(q),q);
%        previous=trapz(previousprofile)*(dx);%area of previous timestep
%        now=trapz(nowprofile)*dx; %area of current timestep
%        volumechange=((now-previous)-BRE)/sc(q);%change in volume is the difference in the area of previousand current time step - the backround erosion
%          
%      elseif sc(q)>=sc(q-1)%shoreline is prograding
%        previousprofile=PROFILEALL(cut:1:scrc(q-1),q-1);
%        nowprofile=PROFILEALL(cut:1:scrc(q-1),q);
%        previous=areareal(q-1);%area of previous timestep
%        now=trapz(nowprofile)*dx;%area of current timestep
%        volumechange=((now-previous)-BRE)/sc(q); %change in volume is the difference in the area of previousand current time step - the backround erosion        
%      end
%% Volume Change Location
     xaxis=cut:1:scrc(q);
     profiledifference=(PROFILEALL(xaxis,q)-PROFILEALL(xaxis,q-1));%area difference between current and previous profile
     
     ased(q)=sum(profiledifference)/scrc(q);
     
     d=profiledifference;%setting deposition to profile difference
     d(d<0)=0;%remove negative values
     e=profiledifference; %set deposition to profile difference
     e(e>=0)=0;%removing positive values
     
     if Zr(q)>Zr(q-1)&& Zr(q)>0 %for sealevel rise PHASE 1
         Drise=d;%find deposition for the current hour
         newX=linspace(0,1,length(Drise));
         Erise=e;%find erosion for the current hour
         newDrise(:,q) = interp1(newX, Drise, newXp);
         newErise(:,q) = interp1(newX, Erise, newXp);
         phase1d(:,q)=newDrise(:,q);
         phase1e(:,q)=newErise(:,q);
%          DR(:,j)=Dr;%save depostion location for each hour 
%          ER(:,j)=Er;%save erosion location for each hour 
     elseif Zr(q)>Zr(q-1)&& Zr(q)<=0 %for sealevel rise PHASE 4
         Drise=d;%find deposition for the current hour
         newX=linspace(0,1,length(Drise));
         Erise=e;%find erosion for the current hour
         newDrise(:,q) = interp1(newX, Drise, newXp);
         newErise(:,q) = interp1(newX, Erise, newXp);
         phase4d(:,q)=newDrise(:,q);
         phase4e(:,q)=newErise(:,q);
     elseif Zr(q)<=Zr(q-1)&& Zr(q)>0 %for sea level fall PHASE 2
         Dfall=d;%find deposition for the current hour
         newX=linspace(0,1,length(Dfall));
         Efall=e;%find erosion for the current hour
         newDfall(:,q) = interp1(newX, Dfall, newXp);
         newEfall(:,q) = interp1(newX, Efall, newXp);
         phase2d(:,q)=newDfall(:,q);
         phase2e(:,q)=newEfall(:,q);
%          DF(:,j)=Df;%save depostion location for each hour 
%          EF(:,j)=Ef;%save erosion location for each hour 
     elseif Zr(q)<=Zr(q-1)&& Zr(q)<=0 %for sea level fall PHASE 3
         Dfall=d;%find deposition for the current hour
         newX=linspace(0,1,length(Dfall));
         Efall=e;%find erosion for the current hour
         newDfall(:,q) = interp1(newX, Dfall, newXp);
         newEfall(:,q) = interp1(newX, Efall, newXp);
         phase3d(:,q)=newDfall(:,q);
         phase3e(:,q)=newEfall(:,q);    
     end
   
%% Volume Change     
     erosion2=profiledifference(profiledifference<0)/sc(q);% takes the negative values from profiledifference and saves as erosion
     deposition2=profiledifference(profiledifference>=0)/sc(q);  %takes the positive values from profiledifference and saves as deposition
     EROSIONstore(q)=sum(erosion2)*dx;%same as volume change erode
     DEPOSITIONstore(q)=sum(deposition2)*dx;%same as volume change deposit
%      VOLUMECHANGE(q)=volumechange;%save volume change at each hour
%      sedimentationrate(q)=volumechange/scrc(q);
end
 
end

plot(ased)
hold on
plot(Zr/15)

%% Figure 3 C and D
%% Derivatives
figure(1)
subplot(3,1,1)
plot(t,-averageprofile1,'k')
xlim([680 1170])
xlabel('Time(Hour)')
title(['First Order Derivative'])

subplot(3,1,2)
plot(t,difference,'b')
hold on
% yline(averageCurvature)
% plot(t,RELIEF,'b')
xlim([680 1170])
xlabel('Time(Hour)')
ylabel('mm')
title(['Departure from the Linear'])

subplot(3,1,3)
plot(t,Zcp-Z,'c')
xlim([680 1170])
xlabel('Time(Hour)')
ylabel('Elevation(mm)')
title(['Sea Level'])

%% Appendix
%% Phases
phase1=(sum(phase1d,2))+(sum(phase1e,2));
phase2=(sum(phase2d,2))+(sum(phase2e,2));
phase3=(sum(phase3d,2))+(sum(phase3e,2));
phase4=(sum(phase4d,2))+(sum(phase4e,2));

%% Appendix 3
% subplot(2,1,1)
% plot(newXp, phase1+phase4,'r')
% title(['Net Sedimentation During Sea-Level Rise'])
% ylabel('mm')
% 
% subplot(2,1,2)
%  plot(newXp, phase2+phase3,'b')
% title(['Net Sedimentation During Sea-Level Fall'])
% ylabel('mm')
% 
% % Movie
% % FALL RISE DEP EROS 
% v = VideoWriter('C:\Users\kolleggerm1\Desktop\movies\SeaLevelFallHMSP.avi');
%    open(v);
% for q = 1:m    
%     if q==1
%     plot(D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q),'w')
%     elseif q>1
%        if Zcp(q)>Zcp(q-1) %sea level rise
%             shade(D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q),D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q-1),'Color','none','FillType',[1,2;2,1],'FillColor', [0 1 0; 1 0 0])
%        if Zcp(q)<Zcp(q-1) %sea level fall
%         shade(D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q),D(cut:1:scrc(q)),PROFILEALL(cut:1:scrc(q),q-1),'Color','none','FillType',[1,2;2,1],'FillColor', [0 1 0; 1 0 0])      
%        end
%     end
%     hold on
%      
%     ylim([160 350])%change if needed
%     xlim([cut*dx 1340])
%     title(['SL Fall  Hour:' num2str(q)])
%     xlabel('Distance (mm)')
%     ylabel('Elevation (mm)')
%     legend('Deposition', 'Erosion')
%     hold off    
%     drawnow
%      frame = getframe(gcf);
%     writeVideo(v,frame);
% end
%      close(v);
%    movie2gif(frame(1:end-1), 'C:\Users\kolleggerm1\Desktop\movies\SLFHMSP.gif', 'LoopCount', inf, 'DelayTime', 0.05)
% 


