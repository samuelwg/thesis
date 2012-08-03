clear;
equator_separation=5; % in degrees
elevation_separation=10;
Nequator=360/equator_separation;

taula(18,Nequator+1)=0;

elev0(Nequator)=0;

n(19)=0;
delta(19)=0;

n(10)=1;
n(19)=1;
taula(10,1)=90;
taula(10,2)=0;
taula(19,1)=-90;
taula(19,2)=0;

for i=1:8
	n(i+1)=round(Nequator*cos( (10*i)/180*pi));
	n(i+10)=n(i+1);
	delta(i+1)=360/n(i+1);	
	delta(i+10)=delta(i+1);
endfor

% from i=1 to i=10, is from 0 to 90.
% from i=11 to i=19, is from -10 to -90.

for i=1:8
	taula(i+1,1)=10*i;
	taula(i+10,1)=-taula(i+1,1);
	for j=1:n(i+1)
		taula(i+1,j+1)=round((j-1)*delta(i+1));
		taula(i+10,j+1)=taula(i+1,j+1);
	endfor
endfor


n(1)=Nequator;
taula(1,1)=0;
for j=1:Nequator
	taula(1,j+1)=(j-1)*equator_separation;
endfor

nPairs=sum(n);
pairs(nPairs,2)=0;

nPairsTemp=0;
for i=1:19
	for j=1:n(i)
		nPairsTemp+=1;
		pairs(nPairsTemp,1)=taula(i,1);
		pairs(nPairsTemp,2)=taula(i,1+j);
	endfor
endfor

