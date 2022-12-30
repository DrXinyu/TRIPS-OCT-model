classdef Vision < handle
    properties
        k
    end
    
    methods
        function self = Vision(k)
            self.k = k;
        end
        
        function waveform(self,Jin)
            
            zi = (0:10:1000*50)';
            for kindex = 1:length(Jin(:,1))
                Ex(:,kindex) = Jin(kindex,1).*exp(1i.*self.k(kindex).*zi);
                Ey(:,kindex) = Jin(kindex,2).*exp(1i.*self.k(kindex).*zi);
            end
            Ex = sum(Ex,2);
            Ey = sum(Ey,2);
            x = zeros(length(zi),1);
            y = zeros(length(zi),1);
            z = zi;
            u = real(Ex);
            v = real(Ey);
            w = zeros(length(zi),1);
            % clear 'Ex';
            % clear 'Ey';

            quiver3(x,y,z,u,v,w,0,'MaxHeadSize',0.001);
        end
        
        function P = powermeter(self,Jin)
            power = Jin.*conj(Jin);
            P = sum(power(:));
            
        end        
        function  plotpoincare(self,S)

            % Stokes vectors transform to cartesian x, y, z in a simple manner:
            warning off
            markersize = 2;
            [m n] = size(S);
            if m ==3 && n==1
                S = S.';
            end
            flg = 1;
            if (m==1)||(n==1)
                markersize = 15;
                flg = 0;
            end
            theta = 0:0.001:2*pi;
            m = length(theta);
            x = cos(theta);
            y = sin(theta);
            z = zeros(1,m);
            H = plot3(x,y,z,'-','color',[0.6 0.5 0.5],'linewidth',1);
            set(H,'linesmoothing','on');
            hold on;
            H = plot3(x,z,y,'color',[0.5 0.5 0.5],'linewidth',1);
            set(H,'linesmoothing','on');
            H = plot3(z,x,y,'color',[0.5 0.5 0.5],'linewidth',1);
            set(H,'linesmoothing','on');

            H = plot3(sqrt(1/3)*x,sqrt(2/3)*x,y,'color',[0.7 0.7 0.7],'linewidth',1);
            set(H,'linesmoothing','on');
            H = plot3(-sqrt(1/3)*x,sqrt(2/3)*x,y,'color',[0.7 0.7 0.7],'linewidth',1);
            set(H,'linesmoothing','on');
            z = repmat(sqrt(0.25),1,m);
            H = plot3(sqrt(0.75)*x,sqrt(0.75)*y,z,'color',[0.7 0.7 0.7],'linewidth',1);
            set(H,'linesmoothing','on');
            H = plot3(sqrt(0.75)*x,sqrt(0.75)*y,-z,'color',[0.7 0.7 0.7],'linewidth',1);
            set(H,'linesmoothing','on');
            [X,Y,Z] = sphere(200);
            re = [0.97 0.95 0.97];
            colormap(re);
            Hs = surf(X,Y,Z);
            % Hs = surf(X,Y,Z,'facecolor','w','edgecolor',[1 1 1]);  % set grid facecolor to white
            shading interp;
            caxis([1.0 1.01]);  % set grid to appear like all one color
            alpha(0.30); 
            axis equal;  % make the three axes equal so the ellipsoid looks like a sphere
            axis off;
            set(gcf,'Renderer','opengl');

            % Draw x- and y- and z-axes
            Hx = plot3([-1.1 1.1], [0 0], [0 0],'k-');
            set(Hx,'linewidth',2,'linestyle','-','color','k','linesmoothing','on');
            ht_x = text(1.3,0,0,'Q','fontweight','bold','fontsize',12,'fontname','arial');
            Hy = plot3([0 0], [-1.1 1.1], [0 0],'k-');
            set(Hy,'linewidth',2,'linestyle','-','color','k','linesmoothing','on');
            ht_y = text(0.1,1.3,0,'U','fontweight','bold','fontsize',12,'fontname','arial');
            Hz = plot3([0 0], [0 0], [-1.1 1.1],'k-');
            set(Hz,'linewidth',2,'linestyle','-','color','k','linesmoothing','on');
            ht_z = text(-0.05,0,1.1,'V','fontweight','bold','fontsize',12,'fontname','arial');
            % ht_lcp = text(-0.05,0.0,1.1,'RCP','fontweight','bold','fontsize',6,'fontname','arial','color','k');
            % Draw a bold circle about the equator (2*epsilon = 0)
            % x_e = (-1:.01:1);
            % for i = 1:length(x_e)
            % z_e(i) = 0;
            % y_e_p(i) = +sqrt(1 - x_e(i)^2);
            % y_e_n(i) = -sqrt(1 - x_e(i)^2);
            % end
            % He = plot3(x_e,y_e_p,z_e,'k-',x_e,y_e_n,z_e,'k-');
            % set(He,'linewidth',2,'color','k');
            % % Draw a bold circle about the prime meridian (2*theta = 0, 180)
            % y_pm = (-1:.01:1);
            % for i = 1:length(x_e)
            % x_pm(i) = 0;
            % z_pm_p(i) = +sqrt(1 - y_pm(i)^2);
            % z_pm_n(i) = -sqrt(1 - y_pm(i)^2);
            % end
            % Hpm = plot3(x_pm,y_pm,z_pm_p,'k-',x_pm,y_pm,z_pm_n,'k-');
            % set(Hpm,'linewidth',2,'color','k');

            % Now plot the polarimetry data
            x = S(:,1);
            y = S(:,2);
            z = S(:,3);
            H = plot3(x,y,z,'m-');
            H2 = plot3(x,y,z,'*');
            set(gca,'fontweight','bold','fontsize',26,'fontname','arial');
            set(H,'markersize',6,'markeredgecolor','m','markerfacecolor','m','color','m','linewidth',1,'linesmoothing','on');
            set(H2,'markersize',markersize,'linewidth',1);
            if flg
                x = S(1,1);
                y = S(1,2);
                z = S(1,3);
                plot3(x,y,z,'k*','markersize',12);
            end
            %H = plot3(x(100:120),y(100:120),z(100:120),'m-.');
            %set(gca,'fontweight','bold','fontsize',12,'fontname','arial');
            %set(H,'markersize',12,'markeredgecolor','m','markerfacecolor','m','color','m','linewidth',2);
            %view(135,20);  % change the view angle
            %alpha(0.70);  % set opacity of sphere to 70%
            %ylim(gca,[0 120]);
            %xlabel('S1');
            %ylabel('S2');
            %zlabel('S3');


            %legend('data');
            %line width of legend
            h = findobj('type', 'axes');  % Find all sets of axes
            set(h(1), 'linewidth',2)  % Making the vertical lines blend in with the background
            set(h(1), 'linewidth',2)  % Making the horizontal lines blend in with the background
            %set(H,'markersize',8,'color','m');

            % Get info for plotting text boxes on screen
            x_rng = xlim;   % range of horizontal axis
            y_rng = ylim;  % range of initial right vertical axis
            delta_x = x_rng(2) - x_rng(1);
            delta_y = y_rng(2) - y_rng(1);

        end


        
        
        
        
    end
end

