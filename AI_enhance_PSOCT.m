clc
clear
close all

addpath('PSOCT_simulation_toolbox_v1')
% ParentFolderPath = fileparts(pwd);
% 
% for index = 1:50
%     [currentPath, FolderName] = fileparts(ParentFolderPath);
%     if strcmp(FolderName,'Dropbox')     
%         toolboxPath = fullfile(currentPath,'Dropbox',...
%           'sspsoct_toolbox');
%         addpath(genpath(toolboxPath));
%         break;
%     end
%     ParentFolderPath = currentPath;
% end


kinterval = 0.0000001;
system = Device(kinterval,[1010 1110]);
see = Vision(system.k);
Jsource = system.source();
Jsource = system.polarizer(Jsource,pi/4);


Jref = Jsource/2;
Jsample = Jsource/2;

Jref = system.mirror(Jref);
Jref = system.space(Jref,100e3);
Jref = system.PC_and_polarizer(Jref,pi/4);


%Jsample1 = system.polarizer(Jsample,0);

Jsample = system.space(Jsample,500e3);


%%
phantom = Sample(system.k);
phantom.layered_sample(500e3,10);

calib = Sample(system.k);
calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],0e3,0);
%calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],100e3,0);
JsampleC = Jsample;
Jcal = calib.backscatter(JsampleC);
Fringe = 1000*system.dual_balanced_polarizaiton_diversity_detector(Jref,Jcal);
alinex = (10*log10(abs(fft(Fringe,4096*4))));
[~,ZeroP] = max(alinex(:,1));


calib = Sample(system.k);
%calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],0e3,0);
calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],100e3,0);
JsampleC = Jsample;
Jcal = calib.backscatter(JsampleC);
Fringe = 1000*system.dual_balanced_polarizaiton_diversity_detector(Jref,Jcal);
alinex = (10*log10(abs(fft(Fringe,4096*4))));
[~,M100P] = max(alinex(:,1));

ZeroP
M100P
zgrid = (M100P-ZeroP)*10*2;
%stop

%%
JsampleGT = Jsample;
Jsample = phantom.backscatter(Jsample);
Jtruth = phantom.backscatter_ground_truth(JsampleGT(2800,:));
Fringe = 1000*system.dual_balanced_polarizaiton_diversity_detector(Jref,Jsample);

% plot(Fringe)
% figure(2)
aline = fft(Fringe,4096*4);

aline = aline(ZeroP:(M100P-ZeroP)*5+ZeroP,:);
%plot(alinex);%(ZeroP:(M100P-ZeroP)*5+ZeroP,:)

% figure(3)
% plot(10*log10(abs(Jtruth)))

S = tom2Stokes(aline(:,1),aline(:,2));
Sn = normalizeStokes(S);
figure(1)
see.plotpoincare(Sn)

Stru = tom2Stokes(Jtruth(:,1),Jtruth(:,2));
Strun = normalizeStokes(Stru);
figure(2)
see.plotpoincare(Strun)





