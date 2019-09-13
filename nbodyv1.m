%n-body simulation
clear all
close all
format compact
format long

G=6.674e-11
global par
par=pparameters
com=0;
mr=0;
mtot=0
for k=2:size(par,1)
   mr=mr+par(k,7)*par(k,9);
   mtot=mtot+par(k,7);
end
com=mr/par(1,7)
mtot=mtot+par(1,7);

mm=par(:,7);
mm=sort(mm);
mm=mm(size(mm,1)-1)

[aa,bb]=find(par==mm)
ee=par(aa,11)

for k=1:size(par,1)
    if par(k,9)==0
        continue
    else
       par(k,15)=sqrt((4*pi^2*par(k,9)^3)/(G*(par(1,7)+par(k,7)))); 
    end
end
par(1,15)=par(aa,15)

for k=1:size(par,1)
    if par(k,9)==0
        par(k,1)=-com;
        par(k,14)=com;
        par(k,11)=ee;
        par(k,9)=com/(1+ee)
        par(k,13)=par(k,9)*(1-ee)
    else
        par(k,1)=par(k,14)*cos(par(k,12)*pi/180);
    end
    par(k,2)=0;
    par(k,3)=par(k,14)*sin(par(k,12)*pi/180);
    par(k,4)=0;
    if par(k,10)==0
        vv=(2*pi*par(k,9))/par(k,15)
        vvv=-sqrt(vv^2*(1-ee)/(1+ee));
        par(k,5)=-sqrt(vv^2*(1-ee)/(1+ee))
        par(k,11)=vv
    else
        par(k,5)=sqrt(G*(par(1,7)+par(k,7))*(2/par(k,14)-1/par(k,9)));
    end
    par(k,6)=0;
end

tfinal=60*60*24*365*1000
tint=60*60*24*200
tstart=0
plottype=1
ti=tstart
tf=tint
objects=size(par,1)
posvel0=[0,0,0,0,0,0] %x,y,z,vx,vy,vz
posvellast=zeros(objects,6)
%%posvelp=zeros(object,12)
%%posvelset=zeros(1,objects)
accuracy=zeros(1,objects)
tdy=sqrt((G*mtot)/(par(1,8)/2)) %sqrt(total mass * G)/1/2*system size)
options=odeset('RelTol',1e-8);
col=['k','r','g','b','m','c','y'];
count=1;
thres=ones(1,objects)
counts=ceil(tfinal/tint)
pp1=ones(counts,objects)
acc2=zeros(counts,objects);

cond=0;
while cond==0
    if objects>length(col)
        col=[col,col];
    else cond=1;
    end
end
col
for k=2:objects
    if k<=4
    thres(k)=(par(4,15)*0.99)/tint;
    else
    thres(k)=(par(k,15)*0.99)/tint;
    end
end
thres(1)=1e12;
global j;
j=1;
figure(1)
plot(par(1,8)*cos(pi*(0:0.01:2))-com,par(1,8)*sin(pi*(0:0.01:2)))
while ti<tfinal
    comx=0; comy=0; comz=0;
    xmr=0; ymr=0; zmr=0; 
    for k=1:size(par,1)
       xmr=xmr+par(k,1);
       ymr=ymr+par(k,2);
       zmr=zmr+par(k,3);
    end
    comx=xmr/mtot;
    comy=ymr/mtot;
    comz=zmr/mtot;
    figure(1)
    hold on
%     plot(comx,comy,'y')
%     hold off
%     figure(2)
%     plot(comx,comz,'y')
%     figure(3)
%     plot(comy,comz,'y')
%     figure(4)
%     plot3(comx,comy,comz,'y')
   while j<=objects
       j
       Ui=0.5*(par(j,4).^2+par(j,5).^2+par(j,6).^2)-(G*par(1,7))/((par(j,1).^2+par(j,2)^.2+par(j,3).^2))^(1/2)
       for p=1:6
           posvel0(p)=par(j,p);
       end
       posvel0
       [t,posvel]=ode45('planeteq',[ti,tf],posvel0,options);
       top=length(t);
       Uf=0.5*(posvel(top,4).^2+posvel(top,5).^2+posvel(top,6).^2)-(G*par(1,7))/((posvel(top,1).^2+posvel(top,2)^.2+posvel(top,3).^2))^(1/2)
       ac=abs(((Uf-Ui)/Ui));
       acc2(count,j)=ac;
       accuracy(1,j)=accuracy(1,j)+ac
       posvellast(j,1)=posvel(top,1)
       posvellast(j,2)=posvel(top,2);
       posvellast(j,3)=posvel(top,3);
       posvellast(j,4)=posvel(top,4);
       posvellast(j,5)=posvel(top,5);
       posvellast(j,6)=posvel(top,6);
       %plots
       j
       if thres(j)<(ti/tint)
          delete(pp1(count-ceil(thres(j))+1,j))
       end
       figure(1)
       hold on
       %x,y
       pp1(count,j)=plot([posvel(1,1),posvel(top,1)],[posvel(1,2),posvel(top,2)],col(j));
       grid on
       hold off
%        figure(2)
%        hold on
%        %x,z
%        plot(posvel(:,1),posvel(:,3),col(j))
%        grid on
%        hold off
%        figure(3)
%        hold on
%        %y,z
%        plot(posvel(:,2),posvel(:,3),col(j))
%        grid on
%        hold off
%        figure(4)
%        hold on
%        %y,z
%        plot3(posvel(:,1),posvel(:,2),posvel(:,3),col(j))
%        grid on
%        hold off
       j=j+1;
   end
   for k=1:objects
       par(k,1)=posvellast(k,1);
       par(k,2)=posvellast(k,2);
       par(k,3)=posvellast(k,3);
       par(k,4)=posvellast(k,4);
       par(k,5)=posvellast(k,5);
       par(k,6)=posvellast(k,6);
   end
   j=1;
   count=count+1
   counts
   ti=ti+tint
   tii=ti/(60*60*24*365.25)
   %pause(1e-9)
   tf=tf+tint;
end
figure(1)
grid on
xlabel('x')
ylabel('y')
% figure(2)
% grid on
% xlabel('x')
% ylabel('z')
% figure(3)
% grid on
% xlabel('y')
% ylabel('z')
% figure(4)
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('z')
tf