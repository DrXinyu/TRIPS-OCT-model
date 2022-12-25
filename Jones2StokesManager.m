classdef Jones2StokesManager < handle

    properties
    end
    
    methods
        
        function self = Jones2StokesManager()
        end
        

        
        function [S1,S2] = exchange(self,J1,J2,JS1,JS2)
            [JS1,JS2] = self.JonesAmpAlign(J1,J2,JS1,JS2);
            S1 = tom2Stokes(J1,J2);
            S2 = tom2Stokes(JS1+J1,JS2+J2);
            
        end        

        
        function [JS1,JS2] = JonesAmpAlign(self,J1,J2,JS1,JS2)

            I = sqrt(abs(J1).^2 + abs(J2).^2);
            IS = sqrt(abs(JS1).^2 + abs(JS2).^2);  
            
            kernel = (gausswin(100,2)./sum(gausswin(100,2)))';
            
            I = imfilter(I,kernel,'replicate');
            IS = imfilter(IS,kernel,'replicate');
            
            mask = (self.SNRmask(I)>3)&(self.SNRmask(IS)>3);
            factor = median(I(mask)./IS(mask));
            JS1 = factor.*JS1;
            JS2 = factor.*JS2;
            
        end
        
        
        function Stokes = tom2Stokes(self,Hc,Vc)
            
            S1 = abs(Hc).^2 - abs(Vc).^2;
            S2 = 2*real(Hc.*conj(Vc));
            S3 = -2*imag(Hc.*conj(Vc));
            dimtom = size(Hc);
            Stokes = squeeze(reshape(cat(4,S1,S2,S3),[dimtom,3]));
            
        end
        
        function SNR = SNRmask(self,I)
            
            dim = size(I);    
            I = reshape(I,dim(1)*dim(2),[]);
            ground_noise = quantile(I,0.1);
            
            SNR = I./ground_noise;
            SNR = reshape(SNR,dim);

        end
    end
end

