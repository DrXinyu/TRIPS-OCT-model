classdef Fiber < handle
    
    properties
        FiberSample
        backward_FiberSample
    end
    
    methods
        function self = Fiber(L,k)
            
            
            c = 3e8*1e9;
            b = 0.1*100; %ps/sqrt(km)
            PMD = b*sqrt(L/1000);
            
            NSeg = 50;
            meanh = (L/NSeg)*1e9;
            
            
            tPMDps = sqrt(3*pi/8)*((L/NSeg)/1000)*b;
            tPMD = tPMDps*1e-12;
            
            
            Fleng = L*1e9*rand(NSeg,1);
            Fleng(NSeg+1) = L*1e9;
            Fleng(1) = 0;
            Fleng = sort(Fleng);
            harray = diff(Fleng);
            
            aarray = rand(NSeg,1).*2*pi;
            detaNa = (tPMD*c)./harray;

            
            for Findex = 1:NSeg
                
                FiberStruct{Findex}.R = [0.00 0.00];
                FiberStruct{Findex}.n = [1.46-detaNa(Findex)/2 1.46+detaNa(Findex)/2];
                FiberStruct{Findex}.z = harray(Findex);
                FiberStruct{Findex}.a = aarray(Findex);  
                FiberStruct{Findex}.loss = [1 1];
            end           
            
            for Findex = NSeg:-1:1
                backward_FiberStruct{NSeg-Findex+1}.R = [0.00 0.00];
                backward_FiberStruct{NSeg-Findex+1}.n = [1.46-detaNa(Findex)/2 1.46+detaNa(Findex)/2];
                backward_FiberStruct{NSeg-Findex+1}.z = harray(Findex);
                backward_FiberStruct{NSeg-Findex+1}.a = pi - aarray(Findex);  
                backward_FiberStruct{NSeg-Findex+1}.loss = [1 1];  
            end
            
            self.FiberSample = Sample(k);
            self.FiberSample.Construct_from_Struct(FiberStruct)
            
            self.backward_FiberSample = Sample(k);
            self.backward_FiberSample.Construct_from_Struct(backward_FiberStruct);     
            

        end
        
        function Jout = through(self,Jin)
            Jout = self.FiberSample.through(Jin);
        end
        
        
        function Jout = backthrough(self,Jin)
            Jout = self.backward_FiberSample.through(Jin);
        end        

        
    end
end

