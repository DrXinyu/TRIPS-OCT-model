clear


kinterval = 0.0000001;
system = Device(kinterval,[1010 1110]);

see = Vision(system.k);


Light = system.source();
Light = system.polarizer(Light,0);



Light_ref = Light/2;
Light_sample = Light/2;

% reference path
Light_ref = system.space(Light_ref,400e3); %unit nanometer
Light_ref = system.mirror(Light_ref);
Light_ref = system.space(Light_ref,400e3); %unit nanometer

% particle 1
Light_sample1 = system.space(Light_sample,458e3); %unit nanometer
Light_sample1 = 0.001*system.mirror(Light_sample1);
Light_sample1 = system.space(Light_sample1,458e3); %unit nanometer

% particle 2
Light_sample2 = system.space(Light_sample,460e3); %unit nanometer
Light_sample2 = 0.001*0.2*system.mirror(Light_sample2);
Light_sample2 = system.space(Light_sample2,460e3); %unit nanometer

Light_sample = Light_sample1+Light_sample2;





% to detector
Light_combined = Light_sample+Light_ref;
Light_combined = system.space(Light_combined,800e3); %unit nanometer

%%

DetectorOutput = Light_combined.*conj(Light_combined);
DetectorOutput = DetectorOutput - Light_ref.*conj(Light_ref);



depth_profile = abs(fft(DetectorOutput,4096*8));

figure(1)
plot(depth_profile)

%see.waveform(Light);



