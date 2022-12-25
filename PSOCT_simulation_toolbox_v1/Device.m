classdef Device < handle

    properties  
        c 
        kinterval
        k
        meank
    end

    
    methods
 
        function self = Device(kinterval,waveband)
            self.c = 3e8*1e9;
            self.kinterval = kinterval;
            self.k = (2*pi/waveband(2):kinterval:2*pi/waveband(1))';
            self.meank = mean(self.k);
        end
        
        
        function Jout = source(self)
            
            spectra = gausswin(length(self.k),1);%%;hanning(length(self.k)); ones(length(self.k),1); 

            phase = 2*pi*rand(length(self.k),1);
            %Yphase = 2*pi*rand(length(self.k),1);
            
            E0x = 1;
            E0y = 1;

            Jout = [E0x.*exp(1i*(phase+self.k.*0)).*spectra, E0y.*exp(1i*(phase+self.k.*0)).*spectra];
            
            

        end
        
        function Jout = PC_and_polarizer(self,Jin,angleP)
            
            dfie = median(angle(Jin(:,1))-angle(Jin(:,2)));
            Jin(:,2) = Jin(:,2).*exp(1i*dfie);
            
            pa = atan(median(abs(Jin(:,2))./abs(Jin(:,1))));
            ra = pa - angleP;
            RotA = [cos(ra) sin(ra)
                    -sin(ra) cos(ra)];   
                
            Jin = (RotA*Jin')';
              
            RotA = [cos(angleP) sin(angleP)
                    -sin(angleP) cos(angleP)];
            iRotA = [cos(-angleP) sin(-angleP)
                -sin(-angleP) cos(-angleP)];
            Ji = (RotA*Jin')';
            Jtx = Ji(:,1);
            Jty = zeros(size(Ji(:,2)));
            Jt = [Jtx Jty];
            Jout = (iRotA*Jt')';
        end        
        
        
        
        function Jout = polarizer(self,Jin,angle)
            RotA = [cos(angle) sin(angle)
                    -sin(angle) cos(angle)];
            iRotA = [cos(-angle) sin(-angle)
                -sin(-angle) cos(-angle)];
            Ji = (RotA*Jin')';
            Jtx = Ji(:,1);
            Jty = zeros(size(Ji(:,2)));
            Jt = [Jtx Jty];
            Jout = (iRotA*Jt')';
        end
         
        function Jout = QWP(self,Jin,angle)
            RotA = [cos(angle) sin(angle)
                    -sin(angle) cos(angle)];
            iRotA = [cos(-angle) sin(-angle)
                -sin(-angle) cos(-angle)];
            Ji = (RotA*Jin')';
            Jtx = Ji(:,1)*exp(1i.*pi/2);
            Jty = Ji(:,2);
            Jt = [Jtx Jty];
            Jout = (iRotA*Jt')';
        end     
        
        function Jout = mirror(self,Jin)
            Jout = [-Jin(:,1) Jin(:,2)];
        end          
        
        
        
        
        function Jout = space(self,Jin,z)
            
            Jtx = Jin(:,1).*exp(1i.*(self.k.*z));
            Jty = Jin(:,2).*exp(1i.*(self.k.*z));
        
            Jout = [Jtx Jty];
            
        end          
        
        function Jout = depth_encoding(self,Jin,dz)
            
            Jtx = Jin(:,1).*exp(1i.*(self.k.*0));
            Jty = Jin(:,2).*exp(1i.*(self.k.*dz));
        
            Jout = [Jtx Jty];   
        end  
        
        function Jout = eo_modulator_forward(self,Jin,V)
            
            Jin = self.QWP(Jin,pi/8);
            angle = pi/4;
            ret = V*(pi/210);
            RotA = [cos(angle) sin(angle)
                    -sin(angle) cos(angle)];
            iRotA = [cos(-angle) sin(-angle)
                -sin(-angle) cos(-angle)];
            Ji = (RotA*Jin')';
            Jtx = Ji(:,1)*exp(1i.*ret);
            Jty = Ji(:,2);
            Jt = [Jtx Jty];
            Jout = (iRotA*Jt')';
            Jout = self.QWP(Jout,pi*(7/8));
            
        end          
        
        function Jout = eo_modulator_backward(self,Jin,V)
            Jin = self.QWP(Jin,pi/8);
            angle = pi-pi/4;
            ret = V*(pi/210);
            RotA = [cos(angle) sin(angle)
                    -sin(angle) cos(angle)];
            iRotA = [cos(-angle) sin(-angle)
                -sin(-angle) cos(-angle)];
            Ji = (RotA*Jin')';
            Jtx = Ji(:,1)*exp(1i.*ret);
            Jty = Ji(:,2);
            Jt = [Jtx Jty];
            Jout = (iRotA*Jt')';
            Jout = self.QWP(Jout,pi*(7/8));
        end   
        
        
        
        
        
        
        
        
        function Fout = dual_balanced_polarizaiton_diversity_detector(self,Jref,Jsample)
            
            realpn = 1400;
            resamplepn = linspace(1,length(self.k),4096);
            
            JcP = 0.717*Jref + 0.717*Jsample.*exp(1i.*pi);
            JcN = 0.717*Jref + 0.717*Jsample;
            Fout = JcP.*conj(JcP) - JcN.*conj(JcN);
            Fout(:,1) = smooth(Fout(:,1),round(length(self.k)/realpn*2));
            Fout(:,2) = smooth(Fout(:,2),round(length(self.k)/realpn*2));
            Fout = interp1((1:length(self.k))',Fout,resamplepn);
            
            Fout = Fout+0.0005*randn(size(Fout));

        end
        
        
        
    end
end

