%%    PART A: design & visualizing
%___________  step-response of an IIR-filter
clear; [B,A]=butter(3,0.1);  % 3rd order LP Butterworth filter

step=[zeros(1,30) ones(1,30)]; % step-function
y1=filter(B,A,step);  % step-response
y2=filtfilt(B,A,step);  % step-response --- zero-phase distortion
figure(1),clf,subplot(3,1,1),plot(step,'.-'),grid,ylim([-0.2 1.2])
subplot(3,1,2),plot([1:60],y1,'r.-',[1:60],step,'k-'),grid,ylim([-0.2 1.2]),title('casaul')
subplot(3,1,3),plot([1:60],y2,'g.-',[1:60],step,'k-'),grid,ylim([-0.2 1.2]),title('bidirectional')

%___ Filter Visualization Tool
fvtool(B,A)


%% Part B  :   BRAINWAVES in alpha band (8-12)Hz
%%
clear,close all
load EEG_data

[B,A]= butter(3,[8 13]/(Fs/2));  %%% Translating the band limits [8-13] to normalized frequencies
fdata=filtfilt(B,A,data')';

Sensor_i=10;
figure(1),clf,subplot(1,2,1),strips(data(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('raw')
              subplot(1,2,2),strips(fdata(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('band-limited')

%%____ reminder about the Hilbert transform
 signal=fdata(3,1:1000); tt=1:1000;  
 Hsignal=hilbert(signal)  % notice the complex numbers

figure(2),clf,subplot(2,1,2),plot(tt,signal,tt,5*angle(Hsignal)+13,'r')
subplot(2,1,1),plot(tt,signal,tt,abs(Hsignal),'r')  % the signal envelope
  

%_____ Analysis of Signal-Dependency
  R=corrcoef(fdata(:,:)') ; 
figure(2),clf,subplot(1,2,1),imagesc(R),colormap(jet),colorbar,title('R(xi,xj)'),axis square   % negative correlations are missed
  subplot(1,2,2),imagesc(abs(R)),colormap(jet),colorbar,title('abs(R(xi,xj))'),axis square % both positive and negative correlations are considered
  
%____ ''order the channels''
figure(2), clf, subplot(1,2,1), plot(topoplot_data(:,1),topoplot_data(:,2),'ko','markerfacecolor','r'),axis square
                for i=1:19; text(topoplot_data(i,1),topoplot_data(i,2)+0.05,num2str(i)), end 
lsensors=[1 3 11 5 13 7 15 9]; rsensors=[2 4 12 6 14 8 16 10]; sensor_list=[lsensors,rsensors,17,18,19]
subplot(1,2,2),plot(topoplot_data(lsensors,1),topoplot_data(lsensors,2),'bo',topoplot_data(rsensors,1),topoplot_data(rsensors,2),'ro',topoplot_data(sensor_list(17:19),1),topoplot_data(sensor_list(17:19),2),'black*')
   plot(topoplot_data(lsensors,1),topoplot_data(lsensors,2),'bo-',topoplot_data(rsensors,1),topoplot_data(rsensors,2),'ro-',topoplot_data(sensor_list(17:19),1),topoplot_data(sensor_list(17:19),2),'black*-')

   %_____________  order the dependency-matrix
   Rordered=R(sensor_list',sensor_list');  
   RR=abs(Rordered); XY=topoplot_data(sensor_list',:); 
   figure(3), clf, subplot(1,2,1),imagesc(abs(R)),colorbar,title('abs(R(xi,xj))'),axis square
                    subplot(1,2,2),imagesc(RR), colormap('jet'),colorbar,title('abs( R(x[i],x[j]) )'),axis square
   
   %____ looking at the sensor-graph
   figure(4),clf, gplot(ones(19,19),XY,'-o')  % fully connected : no info about the dependencies
 
   %___ keeping the most important connections
   threshold =0.7;
   figure(4),clf, subplot(1,3,1),gplot(RR> threshold,XY,'-ok'),axis square, title('strongly connected') 
   subplot(1,3,2),plot(XY(:,1),XY(:,2),'ko'),hold, gplot(Rordered> threshold,XY,'-bo'),axis square,title('correlated') % correlated part
   subplot(1,3,3),plot(XY(:,1),XY(:,2),'ko'),hold, gplot(Rordered<- threshold,XY,'-ro'),axis square,title('anti-correlated') % anticorrelated part
  
  %____ within-hemisphere connectivity
  llist= [ones(1,8), zeros(1,11)];  LC =llist*RR*llist';
  rlist= [zeros(1,8), ones(1,8) ,zeros(1,3)];   RC=rlist*RR*rlist';
  x=(RC-LC)/((RC+LC)/2);  % assymetry index 
  

  %_____ Your Task ....
  % 1) build the assymetry index for all the brainwaves' frequency bands and compare it
  % 2) try to build and inter-hemishesperic index  and compare it across frequency bands
  % 3) perform dependency analysis using the envelope of the signal ...
  
  %  Note: information about the brainwaves' frequency-bands definition can be found here  
  % https://choosemuse.com/blog/a-deep-dive-into-brainwaves-brainwave-frequencies-explained-2/
  
  
  %Task 1)
  % Finding delta band
  clear;
  load EEG_data
  [B,A]= butter(3,[1 4]/(Fs/2));  %%% Translating the band limits [1-4] to normalized frequencies
fdata=filtfilt(B,A,data')';

Sensor_i=10;
figure(1),clf,subplot(1,2,1),strips(data(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('raw')
              subplot(1,2,2),strips(fdata(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('band-limited')

         %%____ reminder about the Hilbert transform
 signal=fdata(3,1:1000); tt=1:1000;  
 Hsignal=hilbert(signal)  % notice the complex numbers

figure(2),clf,subplot(2,1,2),plot(tt,signal,tt,5*angle(Hsignal)+13,'r')
subplot(2,1,1),plot(tt,signal,tt,abs(Hsignal),'r')  % the signal envelope
  

%_____ Analysis of Signal-Dependency
  R=corrcoef(fdata(:,:)') ; 
figure(2),clf,subplot(1,2,1),imagesc(R),colormap(jet),colorbar,title('R(xi,xj)'),axis square   % negative correlations are missed
  subplot(1,2,2),imagesc(abs(R)),colormap(jet),colorbar,title('abs(R(xi,xj))'),axis square % both positive and negative correlations are considered
  
%____ ''order the channels''
figure(2), clf, subplot(1,2,1), plot(topoplot_data(:,1),topoplot_data(:,2),'ko','markerfacecolor','r'),axis square
                for i=1:19; text(topoplot_data(i,1),topoplot_data(i,2)+0.05,num2str(i)), end 
lsensors=[1 3 11 5 13 7 15 9]; rsensors=[2 4 12 6 14 8 16 10]; sensor_list=[lsensors,rsensors,17,18,19]
subplot(1,2,2),plot(topoplot_data(lsensors,1),topoplot_data(lsensors,2),'bo',topoplot_data(rsensors,1),topoplot_data(rsensors,2),'ro',topoplot_data(sensor_list(17:19),1),topoplot_data(sensor_list(17:19),2),'black*')
   plot(topoplot_data(lsensors,1),topoplot_data(lsensors,2),'bo-',topoplot_data(rsensors,1),topoplot_data(rsensors,2),'ro-',topoplot_data(sensor_list(17:19),1),topoplot_data(sensor_list(17:19),2),'black*-')

   %_____________  order the dependency-matrix
   Rordered=R(sensor_list',sensor_list');  
   RR=abs(Rordered); XY=topoplot_data(sensor_list',:); 
   figure(3), clf, subplot(1,2,1),imagesc(abs(R)),colorbar,title('abs(R(xi,xj))'),axis square
                    subplot(1,2,2),imagesc(RR), colormap('jet'),colorbar,title('abs( R(x[i],x[j]) )'),axis square
   
   %____ looking at the sensor-graph
   figure(4),clf, gplot(ones(19,19),XY,'-o')  % fully connected : no info about the dependencies
 
   %___ keeping the most important connections
   threshold =0.7;
   figure(4),clf, subplot(1,3,1),gplot(RR> threshold,XY,'-ok'),axis square, title('strongly connected') 
   subplot(1,3,2),plot(XY(:,1),XY(:,2),'ko'),hold, gplot(Rordered> threshold,XY,'-bo'),axis square,title('correlated') % correlated part
   subplot(1,3,3),plot(XY(:,1),XY(:,2),'ko'),hold, gplot(Rordered<- threshold,XY,'-ro'),axis square,title('anti-correlated') % anticorrelated part
  
  %____ within-hemisphere connectivity
  llist= [ones(1,8), zeros(1,11)];  LC =llist*RR*llist';
  rlist= [zeros(1,8), ones(1,8) ,zeros(1,3)];   RC=rlist*RR*rlist';
  x=(RC-LC)/((RC+LC)/2);  % assymetry index 
  
              
 %Finding theta band
 clear;
  load EEG_data
  [B,A]= butter(3,[4 8]/(Fs/2));  %%% Translating the band limits [4-8] to normalized frequencies
fdata=filtfilt(B,A,data')';

Sensor_i=10;
figure(1),clf,subplot(1,2,1),strips(data(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('raw')
              subplot(1,2,2),strips(fdata(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('band-limited')

              
              
   
 %Finding beta band             
 clear;
  load EEG_data
  [B,A]= butter(3,[13 30]/(Fs/2));  %%% Translating the band limits [13-30] to normalized frequencies
fdata=filtfilt(B,A,data')';

Sensor_i=10;
figure(1),clf,subplot(1,2,1),strips(data(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('raw')
              subplot(1,2,2),strips(fdata(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('band-limited')

              
              
              
 %Finding gamma band      
 clear;
  load EEG_data
  [B,A]= butter(3,[30 45]/(Fs/2));  %%% Translating the band limits [30-45] to normalized frequencies
fdata=filtfilt(B,A,data')';

Sensor_i=10;
figure(1),clf,subplot(1,2,1),strips(data(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('raw')
              subplot(1,2,2),strips(fdata(Sensor_i,1:100*Fs),10,Fs),xlabel('sec'),title('band-limited')
