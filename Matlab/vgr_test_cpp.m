clc;
clear variables;
close all;

cw = 1500;
cb = 2000;
rhob = 2;




%MP = [[0 cw cw 1 1]; [30 cw cb 1 rhob]; [150 cb 1600 rhob rhob]];
MP = [[0 cw cw 1 1]; [90 cw cb 1 rhob]];

z = 0:.25:600;

freqs = 20:5:100;
nmodmax = 5;
vgrs(1:length(freqs),1:nmodmax) = 0;
tic

performc = 0;

if performc
    for ii = 1:length(freqs)
        
        freq1 = freqs(ii);
        [wnum1, wmode1] = ac_modesr(z,MP,freq1);
        
        
        vgr1 = ModesGroupVelocities(z,freq1,wnum1,wmode1,MP);
        
        vgrs(ii, 1:min(nmodmax,length(wnum1)) ) = vgr1(1:min(nmodmax,length(wnum1)));
        
    end;
    toc
    
    dlmwrite('vgrs_matlab.txt',vgrs,'\t');
else
    
    vgrs = dlmread('vgrs_matlab.txt');
end;


figure; 
hold on;
for ii = 1:nmodmax
    
    plot( freqs, vgrs(:,ii),'r:','linewidth',2 );
end;

R = 3500;

dtimes = R./vgrs;
dtimes(dtimes>100) = -1;

dlmwrite('dtimes_synth.txt',[freqs.' dtimes],'\t');

vgrs1 = dlmread('mgv.txt');

for ii = 1:nmodmax
    
    plot( freqs, vgrs1(:,ii) );
end;

ylim([1250 1500]);

error = norm( vgrs1(:,1:nmodmax) - vgrs(:,1:nmodmax) );
disp(error);
