clc
clear
close all

%addpath('PSOCT_simulation_toolbox_v1')

ParentFolderPath = fileparts(pwd);

for index = 1:50
    [currentPath, FolderName] = fileparts(ParentFolderPath);
    if strcmp(FolderName,'Dropbox')     
        toolboxPath = fullfile(currentPath,'Dropbox',...
          'sspsoct_toolbox');
        addpath(genpath(toolboxPath));
        break;
    end
    ParentFolderPath = currentPath;
end


ZeroP = 37
M100P = 60

fract = 5;
wnum = 2*fract - 1;
N = 4096;
W = round(N/fract);
scale = sqrt(sum(hanning(N).^2)/sum(hanning(W).^2));
for wind = 1:wnum
    window(:,wind) = circshift(cat(1,hanning(W),zeros(N-W,1)),round((wind-1)*N/(wnum+1)))*scale;
end
window = permute(shiftdim(window,-1),[2,1,3]);


for a_index = 1:100000


    a_index
    kinterval = 0.0000001;
    system = Device(kinterval,[1010 1110]);
    
    see = Vision(system.k);
    
    Jsource = system.source();
    
    Jsource = system.polarizer(Jsource,pi/4);
    fiber = Fiber(7,system.k);
    Jsource = fiber.through(Jsource);
    Jref = Jsource/2;
    Jsample = Jsource/2;
    
    Jref = system.mirror(Jref);
    Jref = system.space(Jref,100e3); 
    Jref = system.PC_and_polarizer(Jref,pi/4);
    
    
    %Jsample1 = system.polarizer(Jsample,0);
    
    Jsample = system.space(Jsample,500e3);
    
    
    %%
    
    
%     calib = Sample(system.k);
%     calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],0e3,0); % backscattering %0.04 %% 0.001
%     JsampleC = Jsample;
%     Jcal = calib.backscatter(JsampleC);
%     Fringe = 1000000*system.dual_balanced_polarizaiton_diversity_detector(Jref,Jcal);
%     alinex = (10*log10(abs(fft(Fringe,4096))));
% 
%    [~,ZeroP] = max(alinex(:,1));
%     
%     
%     calib = Sample(system.k);
%     %calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],0e3,0);
%     calib.add_layer([0.001 0.001],[0 0],[1.3 1.3],100e3,0);
%     JsampleC = Jsample;
%     Jcal = calib.backscatter(JsampleC);
%     Fringe = 1000*system.dual_balanced_polarizaiton_diversity_detector(Jref,Jcal);
%     alinex = (10*log10(abs(fft(Fringe,4096))));
% 
% 
% 
%     [~,M100P] = max(alinex(:,1));
%     ZeroP
%     M100P
    zgrid = (M100P-ZeroP)*10*2;
    
    %%
    
    
    phantom = Sample(system.k);
    phantom.layered_sample(500e3,10);
    
    
    JsampleGT = Jsample;
    
    Jsample = phantom.backscatter(Jsample);
    bw = round(W/2):round(W/2):round(W/2)*9;
    clear Jtruth
    for bin_index = 1:9
        Jtruth(:,:,bin_index) = phantom.backscatter_ground_truth(JsampleGT(bw(bin_index),:));
    end
    %Jtruth = phantom.backscatter_ground_truth(JsampleGT(2800,:));
    
    Fringe = 1000000*system.dual_balanced_polarizaiton_diversity_detector(Jref,Jsample);

    % plot(Fringe)
    %figure(1)
    aline = fft(Fringe.*window,4096);
    aline = aline(ZeroP:(M100P-ZeroP)*5+ZeroP,:,:);

    %plot(log(abs(aline)));%(ZeroP:(M100P-ZeroP)*5+ZeroP,:)
    
%     figure(3)
%     plot(10*log10(abs(Jtruth)))
    
    S = tom2Stokes(aline(:,1,:),aline(:,2,:));
    S_n = normalizeStokes(S);
    
    % 

%     
    Stru = tom2Stokes(Jtruth(:,1,:),Jtruth(:,2,:));
    S_tru_n = normalizeStokes(Stru);

%     figure(1)
%     see.plotpoincare(squeeze(S_tru_n(:,1,:)))
%     hold on
%     see.plotpoincare(squeeze(S_tru_n(:,5,:)))
%     see.plotpoincare(squeeze(S_tru_n(:,9,:)))



    
    x = linspace(1,length(S_tru_n),length(S_n));
    S_tru_n = interp1((1:length(S_tru_n))',S_tru_n,x);
    S_tru_n = normalizeStokes(S_tru_n);

    S_tru_n_array(:,:,:,a_index) = S_tru_n;
    S_n_array(:,:,:,a_index) = S_n;
    
%     figure(2)
%     see.plotpoincare(S_tru_n)
%     % 
%     figure(4)
%     plot(squeeze(S_n(:,1,:)))
%     hold on
%     plot(squeeze(S_tru_n(:,1,:)))
%     plot(squeeze(S_n(:,5,:)))
%     hold on
%     plot(squeeze(S_tru_n(:,5,:)))
% 
%     stop


end








