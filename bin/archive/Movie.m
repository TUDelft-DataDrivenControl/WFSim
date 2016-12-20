%% Reconstruction:
%clc;
load([Name,'\data'])
load([Name '\It',num2str(1),'.mat'])
cd(Name);
Vmax=max(u_Inf,max(max(u)))
scrsz = get(0,'ScreenSize');
figure('color',[0 166/255 214/255],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)], 'MenuBar','none','ToolBar','none','visible', 'off')
h=gcf; I=2;
aviobj = VideoWriter([Name,'.avi']);
aviobj.FrameRate=10;
open(aviobj)

subplot(4,4,4)
[C,hh]=contour(turbine.Pitch,turbine.TSR,max(turbine.Cp,0)); hold on; clabel(C,hh)
xlabel('Pitch angle')
ylabel('TSR [-]')
title('Cp')
axis([-2 10 4 10])
h33=plot(PitchA(:,:)',(1/2*Omega(:,:).*turbine.Drotor./Veffect(:,:))','k.');

subplot(4,4,8)
[C,hh]=contour(turbine.Pitch,turbine.TSR,max(turbine.Cdax,0)); hold on; clabel(C,hh)
xlabel('Pitch angle')
ylabel('TSR [-]')
title('Ct')
axis([-2 10 4 10])
h36=plot(PitchA(:,:)',(1/2*Omega(:,:).*turbine.Drotor./Veffect(:,:))','k.');

subplot(4,4,13)
h37=plot((2:1:I)*dt,Veffect(:,2:1:I)'); hold on
axis([0 NN*dt 0 u_Inf*1.1])
xlabel('time')
ylabel('Effective windspeed [m/s]')

subplot(4,4,14)
h38=plot((2:1:I)*dt,Power(:,2:1:I)'); hold on
axis([0 NN*dt 0 max(max(Power))*1.1])
xlabel('time')
ylabel('Power [W]')

subplot(4,4,15)
h39=plot((2:1:I)*dt,Omega(:,2:1:I)'); hold on
axis([0 NN*dt 0 max(max(Omega))*1.1])
xlabel('time')
ylabel('Rotor speed [rad/s]')

subplot(4,4,16)
h40=plot((2:1:I)*dt,PitchA(:,2:1:I)'); hold on
axis([0 NN*dt -3 max(max(20))*1.1])
xlabel('time')
ylabel('Pitch Angle [Deg]')

subplot(4,4,12)
h310=plot((2:1:I)*dt,sum(Power(:,2:1:I))'); hold on
axis([0 NN*dt 0 max(max(sum(Power(:,2:1:I))))*1.1])
xlabel('time')
ylabel('Total Power')

for I=1:NN
    load([Name '\It',num2str(I),'.mat'])
    
    
    subplot(4,4,[1 2 3 5 6 7  9 10 11 ])
    contourf(ldyy2(1,:),ldxx2(:,1)',min(sqrt(u.^2+v.^2),u_Inf*1.2),[0:0.05:u_Inf*1.2],'Linecolor','none');  colormap(jet); caxis([0 u_Inf*1.2]);  hold all; colorbar;
    axis equal; axis tight;
    
    for i=1:1:Number_Turbines
        Qy=(Cry(i)-0.5*turbine.Drotor):1:(Cry(i)+0.5*turbine.Drotor);
        plot(Qy, Crx(i)*ones(1,length(Qy)),'k','linewidth',1)
    end
    
    text(0,ldxx2(end,end)+80,['Time ', num2str(dt*I), 's']);
            if streamlines==1;
        for i=1:1:Number_Turbines
            XY = stream2(ldyy2(1,:),ldxx2(:,1)',v,u,[(Cry(i)-0.5*turbine.Drotor);(Cry(i)+0.5*turbine.Drotor)],[Crx(i);Crx(i)]);
            streamline(XY)
            XY1 = stream2(ldyy2(1,:),ldxx2(:,1)',-v,-u,[(Cry(i)-0.5*turbine.Drotor);(Cry(i)+0.5*turbine.Drotor)],[Crx(i);Crx(i)]);
            streamline(XY1)
        end
    end
    xlabel('y-direction [m]')
    ylabel('x-direction [m]')
%     if streamlines==1;
%         for i=1:1:Number_Turbines
%             XY = stream2(dy*clusters:dy*clusters:(Ny*clusters)*dy,dx:dx:(Nx)*dx,v,u,dy*[yline{i}(1);yline{i}(end)],dx*[xline(i,1);xline(i,1)]);
%             streamline(XY)
%             XY1 = stream2(dy*clusters:dy*clusters:(Ny*clusters)*dy,dx:dx:(Nx)*dx,(-v),(-u),dy*[yline{i}(1);yline{i}(end)],dx*[xline(i,1);xline(i,1)]);
%             streamline(XY1)
%         end
%     end
    hold off
    
    for i=1:1:Number_Turbines
        set(h37(i),'xdata',(2:1:I)*dt);
        set(h37(i),'ydata',Veffect(i,2:1:I)');
        set(h38(i),'xdata',(2:1:I)*dt);
        set(h38(i),'ydata',Power(i,2:1:I)');
        set(h39(i),'xdata',(2:1:I)*dt);
        set(h39(i),'ydata',Omega(i,2:1:I)');
        
        set(h40(i),'xdata',(2:1:I)*dt);
        set(h40(i),'ydata',PitchA(i,2:1:I)');
        
        set(h33(i),'xdata',PitchA(i,I));
        set(h33(i),'ydata',1/2*Omega(i,I).*turbine.Drotor./Veffect(i,I));
        
        
        set(h36(i),'xdata',PitchA(i,I));
        set(h36(i),'ydata',1/2*Omega(i,I).*turbine.Drotor./Veffect(i,I));
    end
    set(h310(1),'xdata',(2:1:I)*dt);
    set(h310(1),'ydata',sum(Power(:,2:1:I))');
    
    F=getframe(h);
    
    writeVideo(aviobj,F);
    
end

close(aviobj);

cd ..

