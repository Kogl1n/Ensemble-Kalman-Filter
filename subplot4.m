close all;
load('5w-h08-N272-001+param2-04812-0005-0005-Re011-sigma01','u_ens','w','w_gen','w_obs');

for i = 1:N_ens % kann parallelisiert werden
        u_last(1,:) = u_ens(:,1,8,i); % nehme den letzten Zustand des Ensembles wieder auf.
        u_last(2,:) = u_ens(:,2,8,i);
        u02 = @(locations) u_last;
        
        %delete(model.InitialConditions); % lösche dazu den Anfangswert
        setInitialConditions(model, u02); % und lege den zum letzten Zustand als neuen Anfangswert an
        tlist1 = 7:1:8;
        
%        global rho0_short;
%        global B0_short; 
%        global eta_short;
%        global A0_short;
%        global B_asterisk_short;
% 
%        rho0_short = 1.2*rho0;
%        B0_short = 1.2*B0; 
%        eta_short = 1.2*eta;
%        A0_short = 0.8*A0;
%        B_asterisk_short = 1.2*B_asterisk;
%        specifyCoefficients(model, 'm', 0, 'd', 1, 'c', @ccoeffunction, 'a', @acoeffunction, 'f', @fcoeffunction);
%         
        g=sprintf('-%d-', tlist1);
        fprintf('Berechne für Mitglied %d den Forecast im Zeitschritt %s\n',i,g);
        result2 = solvepde(model, tlist1);     
        
        u_ens(:,:,9,i) = result2.NodalSolution(:,:,2);
        % Modellrauschen + Clamp:
        u_ens(:,1,9,i) = max(min(u_ens(:,1,9,i) + mvnrnd(zeros(size_ugen(1),1), Q1)',1),0);
        u_ens(:,2,9,i) = max(min(u_ens(:,2,9,i) + mvnrnd(zeros(size_ugen(1),1), Q2)',1),0);
         
        clear result2;     
end

u_modell2(:,:,9) = mean(u_ens(:,:,9,:),4);

w_modell(:,:)=w(:,:);
w_modell(:,9)=u_modell2(:,1,9).*u_modell2(:,2,9);

for i=9
    figure

    %suptitle('Wahrscheinlichkeitsdichte zu T=8');

    ah(1)=subplot(2,2,1);
    pdeplot(model, 'xydata', w_gen(:,9), 'contour', 'off', 'colorbar', 'off');
    title(['wahr'])
    ylabel('y-Koord.')
    axis equal
    axis([0 10 0 10])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])

    
    ah(4)=subplot(2,2,4);
    pdeplot(model, 'xydata', w(:,9), 'contour', 'off', 'colorbar', 'off');
    title(['assimiliert'])
    xlabel('x-Koordinate')
    axis equal
    axis([0 10 0 10])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    ah(3)=subplot(2,2,3);
    pdeplot(model, 'xydata', w_modell(:,9), 'contour', 'off', 'colorbar', 'off');
    title(['Modell vor Update'])
    xlabel('x-Koordinate')
    ylabel('y-Koord.')
    axis equal
    axis([0 10 0 10])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[]) 
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    ah(2)=subplot(2,2,2);
    pdeplot(model, 'xydata', w_obs(:,9), 'contour', 'off', 'colorbar', 'off');
    title(['observiert'])
    axis equal
    axis([0 10 0 10])  
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

pos1 = get(ah(1),'Position');
pos2 = get(ah(2),'Position');
pos3 = get(ah(3),'Position');
pos4 = get(ah(4),'Position');


pos2(1) = pos1(1)+pos2(3)-0.077;
set(ah(2),'Position',pos2)

pos3(2) = pos1(2)-pos1(4)-0.04;
set(ah(3),'Position',pos3)


  
pos1 = get(ah(1),'Position');
pos2 = get(ah(2),'Position');
pos3 = get(ah(3),'Position');
pos4 = get(ah(4),'Position'); 
pos4(1) = pos2(1)+0.00020;
pos4(2) = pos3(2);
set(ah(4),'Position',pos4)

colorbar('Position', [pos1(1)+2*pos1(3)-0.1 pos3(2) 0.05 pos3(2)+2*pos1(3)-0.15]);
caxis([0, 1]);
colormap jet


axes( 'Position', [0, 0.98, 0.9, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Wahrscheinlichkeitsdichte zu T=8', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center') ;