classdef Sample < handle

    
    properties
        layer_num
        sstruct
        k
        zgrid
        
    end
    
    methods
        function self = Sample(k)
            self.layer_num = 0;      
            self.k = k;
            self.zgrid = 800;%800;%400;
            
        end
        
        function Construct_from_Struct(self,S)
            self.layer_num = length(S);
            self.sstruct =S;           
        end
        
        
        function  add_layer(self,reflectivity,loss,refraction,thickness,rotation)
            
            self.layer_num = 1 + self.layer_num ;
            self.sstruct{self.layer_num}.R = reflectivity;
            self.sstruct{self.layer_num}.loss = 1 - loss;
            self.sstruct{self.layer_num}.n = refraction;
            self.sstruct{self.layer_num}.z = thickness;
            self.sstruct{self.layer_num}.a = rotation;    

        end
        
        
        
        
        
        
        function Jout = through(self,Jin)
            Jt = Jin;
            for sindex = 1:self.layer_num   
%                 if(isfield(S{sindex},'cn'))
% 
%                     Cd = [0.5 -0.5*1i
%                         -0.5*1i 0.5];
%                     CJt = (Cd*Ji')';
% 
%                     CJtx = CJt(:,1).*exp(1i.*(k.*S{sindex}.z.*S{sindex}.cn(1)));
%                     CJty = CJt(:,2).*exp(1i.*(k.*S{sindex}.z.*S{sindex}.cn(2)));
% 
%                     CJt = [CJtx CJty];
%                     if(isfield(S{sindex},'cd'))
%                         CJt = [CJtx CJty].*S{sindex}.cd;
%                     end
%                     Jt = ([1 1i; 1i 1]*CJt')';
%                 end    
                RotA = [cos(self.sstruct{sindex}.a) sin(self.sstruct{sindex}.a)
                    -sin(self.sstruct{sindex}.a) cos(self.sstruct{sindex}.a)];

                iRotA = [cos(-self.sstruct{sindex}.a) sin(-self.sstruct{sindex}.a)
                    -sin(-self.sstruct{sindex}.a) cos(-self.sstruct{sindex}.a)];

                Jt = (RotA*Jt')';        


                Jtx = Jt(:,1).*exp(1i.*(self.k.*self.sstruct{sindex}.z.*self.sstruct{sindex}.n(1)));
                Jty = Jt(:,2).*exp(1i.*(self.k.*self.sstruct{sindex}.z.*self.sstruct{sindex}.n(2)));


                Jt = [(1-self.sstruct{sindex}.R(1)).*Jtx.*self.sstruct{sindex}.loss(1) (1-self.sstruct{sindex}.R(2)).*Jty.*self.sstruct{sindex}.loss(2)];        
%                 if(isfield(S{sindex},'ld'))
%                     Jt = Jt.*S{sindex}.ld;
%                 end
                Jout = (iRotA*Jt')';
            end        
        end

        function Jout = backscatter_ground_truth(self,Jin)

            Jt = Jin;
            k = self.k(round(length(self.k)/2));
            zgrid = self.zgrid;
            for Sindex = 1:round(length(self.sstruct))
                RotA = [cos(self.sstruct{Sindex}.a) sin(self.sstruct{Sindex}.a)
                        -sin(self.sstruct{Sindex}.a) cos(self.sstruct{Sindex}.a)];

                iRotA = [cos(-self.sstruct{Sindex}.a) sin(-self.sstruct{Sindex}.a)
                        -sin(-self.sstruct{Sindex}.a) cos(-self.sstruct{Sindex}.a)];

                Jt = (RotA*Jt')';

                Jtx = Jt(:,1).*exp(1i.*(k.*zgrid.*self.sstruct{Sindex}.n(1)));
                Jty = Jt(:,2).*exp(1i.*(k.*zgrid.*self.sstruct{Sindex}.n(2)));

                Jt =  [(1-self.sstruct{Sindex}.R(1)).*Jtx.*self.sstruct{Sindex}.loss(1) (1-self.sstruct{Sindex}.R(2)).*Jty.*self.sstruct{Sindex}.loss(2)];
                Jt = (iRotA*Jt')';

                Jrx = -self.sstruct{Sindex}.R(1).*Jtx.*exp(1i.*(k.*zgrid.*self.sstruct{Sindex}.n(1)));
                Jry = self.sstruct{Sindex}.R(2).*Jty.*exp(1i.*(k.*zgrid.*self.sstruct{Sindex}.n(2)));               

                Jr = [Jrx, Jry];  
                Jr = (RotA*Jr')'; 

                for iSindex = Sindex-1:-1:1
                   
                    RotA = [cos(self.sstruct{iSindex}.a) sin(self.sstruct{iSindex}.a)
                            -sin(self.sstruct{iSindex}.a) cos(self.sstruct{iSindex}.a)];
    
                    iRotA = [cos(-self.sstruct{iSindex}.a) sin(-self.sstruct{iSindex}.a)
                            -sin(-self.sstruct{iSindex}.a) cos(-self.sstruct{iSindex}.a)];                    
                    
                    Jr =  (iRotA*Jr')';      
                    Jrx = Jr(:,1).*exp(1i.*(k.*zgrid.*self.sstruct{iSindex}.n(1)));
                    Jry = Jr(:,2).*exp(1i.*(k.*zgrid.*self.sstruct{iSindex}.n(2)));   
    
                    Jr =  [(1-self.sstruct{iSindex}.R(1)).*Jrx.*self.sstruct{iSindex}.loss(1) (1-self.sstruct{iSindex}.R(2)).*Jry.*self.sstruct{iSindex}.loss(2)];
                    Jr = (RotA*Jr')';

                end
                Jout(Sindex,:)=Jr;
            end
            
        end
        function t = thickness(self)
            t = 0;
            for Sindex = 1:length(self.sstruct)
                t = t+self.sstruct{Sindex}.z;
            end
        end

        
        function Jout = backscatter(self,Jin)
            
            Jt = Jin;
            Rc = cell(length(self.sstruct),1);
            for Sindex = 1:length(self.sstruct)
                RotA = [cos(self.sstruct{Sindex}.a) sin(self.sstruct{Sindex}.a)
                        -sin(self.sstruct{Sindex}.a) cos(self.sstruct{Sindex}.a)];

                iRotA = [cos(-self.sstruct{Sindex}.a) sin(-self.sstruct{Sindex}.a)
                        -sin(-self.sstruct{Sindex}.a) cos(-self.sstruct{Sindex}.a)];

                Jt = (RotA*Jt')';

                Jtx = Jt(:,1).*exp(1i.*(self.k.*self.sstruct{Sindex}.z.*self.sstruct{Sindex}.n(1)));
                Jty = Jt(:,2).*exp(1i.*(self.k.*self.sstruct{Sindex}.z.*self.sstruct{Sindex}.n(2)));

                Jt = [(1-self.sstruct{Sindex}.R(1)).*Jtx.*self.sstruct{Sindex}.loss(1) (1-self.sstruct{Sindex}.R(2)).*Jty.*self.sstruct{Sindex}.loss(2)];
                Jt = (iRotA*Jt')';

                Jrx = -self.sstruct{Sindex}.R(1).*Jtx.*exp(1i.*(self.k.*self.sstruct{Sindex}.z.*self.sstruct{Sindex}.n(1)));
                Jry = self.sstruct{Sindex}.R(2).*Jty.*exp(1i.*(self.k.*self.sstruct{Sindex}.z.*self.sstruct{Sindex}.n(2)));               

                Jr = [Jrx, Jry];  
                Jr = (RotA*Jr')'; 
                Rc{Sindex}.Jr = Jr;
            end

            for Sindex = 1:length(self.sstruct)-1

                iSindex = length(self.sstruct)-Sindex;
                RotA = [cos(self.sstruct{iSindex}.a) sin(self.sstruct{iSindex}.a)
                        -sin(self.sstruct{iSindex}.a) cos(self.sstruct{iSindex}.a)];

                iRotA = [cos(-self.sstruct{iSindex}.a) sin(-self.sstruct{iSindex}.a)
                        -sin(-self.sstruct{iSindex}.a) cos(-self.sstruct{iSindex}.a)];

                Jr =  (iRotA*Jr')';      
                Jrx = Jr(:,1).*exp(1i.*(self.k.*self.sstruct{iSindex}.z.*self.sstruct{iSindex}.n(1)));
                Jry = Jr(:,2).*exp(1i.*(self.k.*self.sstruct{iSindex}.z.*self.sstruct{iSindex}.n(2)));   

                Jr = [(1-self.sstruct{iSindex}.R(1)).*Jrx.*self.sstruct{iSindex}.loss(1) (1-self.sstruct{iSindex}.R(2)).*Jry.*self.sstruct{iSindex}.loss(2)];
                %Jr = [Jrx Jry];
                Jr = (RotA*Jr')';

                Jr = Jr + Rc{iSindex}.Jr;
            end 
            
            Jout = Jr;
            
        end        
        
        
        
        
        function full_speckle_layer(self,scattering,thickness,deltan,rotation)
            
            N = round(thickness/self.zgrid);
            %deltan = deltan*(1e3/self.zgrid);
            Sleng = thickness*rand(N,1);
            Sleng(N+1) = thickness;
            Sleng(1) = 0;
            Sleng = sort(Sleng);
            harray = diff(Sleng);
            
            for Sindex = 1:N
                self.add_layer([scattering scattering],[0 0.000],[1.3 1.3+deltan],harray(Sindex),rotation)
            end
        end

        function layered_sample(self,thickness,Nlayer)

            Sleng = thickness*rand(Nlayer,1);
            Sleng(Nlayer+1) = thickness;
            Sleng(1) = 0;

            Sleng = sort(Sleng);
            harray = diff(Sleng);

            deltaN = rand(Nlayer,1)*0.001;
            oa = rand(Nlayer,1)*pi;
            inte = 10.^(rand(Nlayer,1)*3)./(10.^6);


            for ind = 1:Nlayer
                self.full_speckle_layer(inte(ind),harray(ind),deltaN(ind),oa(ind));
            end

        end

        
    end
           
end

