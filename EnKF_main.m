%% MAIN SCRIPT
% Anwendung des Ensemble-Kalman-Filters auf das Modell von Short et al.

N_ens = 272; % Anzahl der Ensemble-Mitglieder

%% Anfangswerte und Parameter des Modells
B0 = 0.5; % Anfangswert, unverfälscht
rho0 = 0.1; % Anfangswert, unverfälscht

eta = 0.4;
A0 = 0.1;
B_asterisk = 0.5;

%% Zeitangaben
t_start = 0; % Startzeit
t_delta = 1; % Länge Zeitschritt
N_delta_obs = 4; % Zeitschritte zwischen einzelnen Observationen (std.4)
t_delta_obs = t_delta*N_delta_obs; % Zeit zwischen Observationen
N_obs = 3; % Gesamtzahl Observationen (std.2)

t_end = t_start+t_delta_obs*N_obs; % Endzeit
t_num = N_delta_obs*N_obs;

% Hier stand einmal sehr viel Nebenrechnung, da wir t_end zuerst vorgegeben
% hatten. Dies benötigt jedoch die (symbolische) Berechnung des kgV.

%% Wertegenerierung
% Generiere plausible Werte u_gen für "genuine/generated"

% Simulation durchführen mit Anfangswerten rho0 und B0
[model, u_gen] = Short(rho0, B0, eta, A0, B_asterisk, t_start, t_delta, t_end);
size_ugen = size(u_gen);
u_gen = max(min(u_gen, 1), 0);

%% Preprocessing EnKF
% Observationsmatrix
H = eye(272); % Diese würde bei Skalarwerten die Einheitsmatrix 
% (#Modellgl. x #Werte) sein, 
% da wir die Observationen von Hand erzeugt haben
% Bei einer 2D-Domain bei 2 Modellgleichungen
% H*x_obs


% Observations-Covarianz-Matrix
R = eye(272);
% Index 1 ist B, Index 2 ist Rho
H1=H;H2=H1;
R1 = 0.005*R; % 0.005
R2 = 0.005*R; % 0.005
% Modell-Covarianz-Matrix
Q1 = 0.01*R; % 0 oder 0.01
Q2 = Q1;

% Observationen sind die echten Werte, 
% versehen mit Rauschen der Form N(0,sigma)
for i=1:(t_num+1)
    u_obs(:,1,i) = u_gen(:,1,i)+mvnrnd(zeros(size_ugen(1),1), R1)'; %mvnrnd statt lognrnd
    u_obs(:,2,i) = u_gen(:,2,i)+mvnrnd(zeros(size_ugen(1),1), R2)'; %mvnrnd statt lognrnd
end
% sogenannte Clamp-Funktion, um professionellerweise seltene Dichtewerte
% unter 0 und über 1 auf 0 bzw. 1 zu setzen
% Wir brauchen ja auch negative Werte, sonst wäre lognrnd angebracht.
u_obs = max(min(u_obs, 1), 0);

% Lösung vordefinieren
u = zeros(size_ugen);
u(:,1,1) = u_gen(:,1,1); % Anf.werte
u(:,2,1) = u_gen(:,2,1);

%% Ensemble für EnKF erzeugen

sigma = 0.1; % Standardabweichung des Ensembles
u_ens = EnKF_ens_create(u_gen(:,1:2,1), N_ens, sigma); %
% u_ens ist ein 3D-Array mit in der Tiefe positionierten Ensemble-Matrizen 
% aus Anfangswertspalten des Modells
% bereits geclamped

%% Ensemble-Kalman-Filter: Hauptschleife
u(:,:,1) = mean(u_ens, 4); 

for t = 1:N_obs % Iteriere über die Observationen
    T=t_start+(t-1)*t_delta_obs;
    if t==1
        cpu_time = cputime; % Achtung: Bei paralleler Ausführung wird dies je Thread zusammengezählt.
    end
    % Simulation für die versch. Ensemblemitglieder
    
    for j = 1:N_delta_obs
    % falls Modellrauschen manuell eingefügt werden
    % soll, muss t_delta = 1 gelten und jeder Schritt einzeln berechnet
    % werden
    for i = 1:N_ens % kann parallelisiert werden
%         u_last(1,:) = u_ens(:,1,(t-1)*N_delta_obs+1,i); % nehme den letzten Zustand des Ensembles wieder auf.
%         u_last(2,:) = u_ens(:,2,(t-1)*N_delta_obs+1,i);

        u_last(1,:) = u_ens(:,1,(t-1)*N_delta_obs+j,i); % nehme den letzten Zustand des Ensembles wieder auf.
        u_last(2,:) = u_ens(:,2,(t-1)*N_delta_obs+j,i);
        u02 = @(locations) u_last;
        
        delete(model.InitialConditions); % lösche dazu den Anfangswert
        setInitialConditions(model, u02); % und lege den zum letzten Zustand als neuen Anfangswert an
        % Forecast: PDE-System-Lösung als Blackbox
%         tlist1 = T:t_delta:(T+t_delta_obs);
        tlist1 = T:t_delta:(T+t_delta);
        
       % für Verfälschung der Parameter
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
        
        g=sprintf('-%d-', tlist1);
        fprintf('Berechne für Mitglied %d den Forecast im Zeitschritt %s\n',i,g);
        result2 = solvepde(model, tlist1);     
        
%        u_ens(:,:,(t-1)*N_delta_obs+2:t*N_delta_obs+1,i) = result2.NodalSolution(:,:,2:end);
        u_ens(:,:,(t-1)*N_delta_obs+j+1,i) = result2.NodalSolution(:,:,2);
        % Modellrauschen + Clamp:
        u_ens(:,1,(t-1)*N_delta_obs+j+1,i) = max(min(u_ens(:,1,(t-1)*N_delta_obs+j+1,i) + mvnrnd(zeros(size_ugen(1),1), Q1)',1),0);
        u_ens(:,2,(t-1)*N_delta_obs+j+1,i) = max(min(u_ens(:,2,(t-1)*N_delta_obs+j+1,i) + mvnrnd(zeros(size_ugen(1),1), Q2)',1),0);
         
        clear result2;     
    end
    T=T+t_delta; % entfernen falls keine j-Schleife!
    end
    
    if t==1
        cpu_time2 = cputime-cpu_time;
    end
    % EnKF-Update:
    % Der derzeitige Zeitpunkt bzgl. t_delta ist ein Vielfaches vom
    % Observations-Zeitschritt bzgl. t_delta_obs, da wir 
    % über die Observationen iterieren. Sonst müssten wir dies über den
    % Rest einer Division prüfen
    
    u_modell(:,1,t*N_delta_obs+1,:) = u_ens(:,1,t*N_delta_obs+1,:);
    fprintf('Führe den EnKF-Update-Schritt aus...\n');
     [u_ens(:,1,t*N_delta_obs+1,:),K1,Pa_e1]=EnKF(u_ens(:,1,t*N_delta_obs+1,:),u_obs(:,1,t*N_delta_obs+1),H1,R1);
     [u_ens(:,2,t*N_delta_obs+1,:),K2,Pa_e2]=EnKF(u_ens(:,2,t*N_delta_obs+1,:),u_obs(:,2,t*N_delta_obs+1),H2,R2);
    if t==1
        cpu_time3 = cputime-cpu_time2-cpu_time;
    end   
%     for i=1:N_ens % Gauss-Filter zur Glättung
%         u_ens(:,1,t*N_delta_obs+1,i)=imgaussfilt(u_ens(:,1,t*N_delta_obs+1,i),2);
%         u_ens(:,2,t*N_delta_obs+1,i)=imgaussfilt(u_ens(:,2,t*N_delta_obs+1,i),2);
%     end
    % Bilde den Mittelwert der Ensemble-Werte
    u(:,:,(t-1)*N_delta_obs+2:t*N_delta_obs+1) = mean(u_ens(:,:,(t-1)*N_delta_obs+2:t*N_delta_obs+1,:), 4); 
    u_modell2(:,:,(t-1)*N_delta_obs+2:t*N_delta_obs+1) = mean(u_modell(:,:,(t-1)*N_delta_obs+2:t*N_delta_obs+1,:), 4); 
    % da Ensemble entlang der 4. Dimension liegt
end


%% Lösungen plotten

for i=0:t_num
    figure
    subplot(2,1,1) % Attraktivität
    pdeplot(model, 'xydata', u(:,1,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
    title(['assim. Attraktivität zu T=',num2str(t_end-i*t_delta)])
    xlabel('x-Koordinate')
    ylabel('y-Koordinate')
    axis equal
    caxis([0, 1]);
    
    subplot(2,1,2) % Einbrecher
    pdeplot(model, 'xydata', u(:,2,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
    title(['assim. Verbrecherdichte zu T=',num2str(t_end-i*t_delta)])
    xlabel('x-Koordinate')
    ylabel('y-Koordinate')
    axis equal
    caxis([0, 1]);
end
% Test mit 
% tlist1 = 0:1:5;
% result2 = solvepde(model, tlist1); 
% u2(:,:,1:6) = result2.NodalSolution;
% und ohne KF ergibt identisches Ergebnis

% assimilierte Wahrscheinlicheitsdichte
w(:,:)=u(:,1,:).*u(:,2,:);
% for i=0:t_num
%     figure
%     pdeplot(model, 'xydata', w(:,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
%     title(['assimilierte Wahrscheinlichkeitsdichte zu T=',num2str(t_end-i*t_delta)])
%     xlabel('x-Koordinate')
%     ylabel('y-Koordinate')
%     axis equal
%     caxis([0, 1]);
% end

%% Plotten und mittleren quadr. Fehler berechnen

% generierte Wahrscheinlicheitsdichte
w_gen(:,:) = u_gen(:,1,:).*u_gen(:,2,:);
% for i=0:t_num
%     figure
%     pdeplot(model, 'xydata', w_gen(:,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
%     title(['wahre Wahrscheinlicheitsdichte zu T=',num2str(t_end-i*t_delta)])
%     xlabel('x-Koordinate')
%     ylabel('y-Koordinate')
%     axis equal
%     caxis([0, 1]);
% end

% observierte Wahrscheinlicheitsdichte
w_obs(:,:) = u_obs(:,1,:).*u_obs(:,2,:);
% for i=0:t_num
%     figure
%     pdeplot(model, 'xydata', w_obs(:,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
%     title(['observierte Wahrscheinlicheitsdichte zu T=',num2str(t_end-i*t_delta)])
%     xlabel('x-Koordinate')
%     ylabel('y-Koordinate')
%     axis equal
%     caxis([0, 1]);
% end

for i=0:t_num
    figure
    subplot(3,1,1);
    pdeplot(model, 'xydata', w_gen(:,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
    title(['wahre Wahrscheinlichkeitsdichte zu T=',num2str(t_end-i*t_delta)])
    xlabel('x-Koordinate')
    ylabel('y-Koord.')
    axis equal
    axis([0 10 0 10])
    caxis([0, 1]);

    
    subplot(3,1,2);
    pdeplot(model, 'xydata', w(:,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
    title(['assimilierte Wahrscheinlichkeitsdichte zu T=',num2str(t_end-i*t_delta)])
    xlabel('x-Koordinate')
    ylabel('y-Koord.')
    axis equal
    axis([0 10 0 10])
    caxis([0, 1]);
    
    subplot(3,1,3);
    pdeplot(model, 'xydata', w_obs(:,t_num-i+1), 'contour', 'off', 'colormap', 'jet');
    title(['observierte Wahrscheinlichkeitsdichte zu T=',num2str(t_end-i*t_delta)])
    xlabel('x-Koordinate')
    ylabel('y-Koord.')
    axis equal
    axis([0 10 0 10])
    caxis([0, 1]);
end

w_ens=zeros(size_ugen(1),size_ugen(3),N_ens);
w_ens(:,:,:) = u_ens(:,1,:,:).*u_ens(:,2,:,:); % Dim Ort x Zeit x Ens

abwq = (w-w_gen).^2;

% MSE/Var ort-gemittelt
mseort = mean(abwq,1);

Pa_e = zeros(size(w_ens,1), size(w_ens,1), t_num+1);
hilfsvar = zeros(size(w_ens,1),size(w_ens,3));
for i=1:(t_num+1)
    hilfsvar(:,:) = w_ens(:,i,:);
    Pa_e(:,:,i) = cov(hilfsvar');
    mvar(i) = mean(diag(Pa_e(:,:,i)));
end
figure
plot([1:(t_num+1)],mseort(1,:),'b-o')
hold on
plot([1:(t_num+1)],mvar,'r-o')
title('zeitlicher Verlauf des MSE und Varianz über den Ort gemittelt');
ylabel('MSE bzw. Var')
xlabel('Zeitschrittindex')
legend('MSE','Var')

% MSE (elementweise) zeitlich gemittelt + 
% mittlere Update-Ensemble-Varianz

msezeit = mean(abwq,2);
Paezeit = mean(Pa_e,3);

figure
subplot(2,1,1)
pdeplot(model, 'xydata', msezeit, 'contour', 'off', 'colormap', 'jet');
title('(elementweiser) quadratischer Fehler zeitlich gemittelt')
xlabel('x-Koordinate')
ylabel('y-Koordinate')
axis equal
caxis([0, 1]);

subplot(2,1,2)
pdeplot(model, 'xydata', Paezeit, 'contour', 'off', 'colormap', 'jet');
title('(elementweise) Varianz zeitlich gemittelt')
xlabel('x-Koordinate')
ylabel('y-Koordinate')
axis equal
caxis([0, 1]);

% MSE/Var zeit- und ort-gemittelt
msezeitort = mean(msezeit,1);
varzeitort = mean(diag(Paezeit),1);

% Histogramm verschiedene absoluter Abweichungen, Zeit gemittelt
abw = (w-w_gen);
abszeit = mean(abw,2);
figure
histogram(abszeit,30);
title('Histogramm des absoluten Fehlers über die Zeit gemittelt');
ylabel('Anzahl')

% MSE/Var vs. Ensemble-Größe
% figure
% Zahl = [5,10];
% % vecmsezeitort = [0.0111,0.0097];
% % vecvarzeitort = [0.0226,0.0226];
% % plot(Zahl,vecmsezeitort,'b-o')
% hold on
% % plot(Zahl,vecvarzeitort,'r-o')
% xlabel('Mitglied')
% title('quadr. Fehler und Varianz (zeitlich und örtlich gemittelt) in Abh. von der Ens.-Größe');
% ylabel('MSE bzw. Var')
% legend('MSE','Var')


% Verlauf der quadratischen Abweichung und Varianz je Ort über Zeit
figure
subplot(2,1,1)
imagesc(abwq)
colormap('jet')
colorbar
caxis([0, 1]);
title('zeitlicher Verlauf der quadratischen Abweichung je Ort');
xlabel('Zeitschrittindex')
ylabel('Ort')

for i=1:(t_num+1)
    VarVerlauf(:,i) = diag(Pa_e(:,:,i));
end
subplot(2,1,2)
imagesc(VarVerlauf)
colormap('jet')
colorbar
caxis([0, 1]);
title('zeitlicher Verlauf der Varianz je Ort');
xlabel('Zeitschrittindex')
ylabel('Ort')

%-------------------------------------------------------------------------
% quadr. Abweichung per Ensemble-Mitglied
for i=1:N_ens
    abwqens(:,:,i) = (w_ens(:,:,i)-w_gen(:,:)).^2;
end
mseensort = mean(abwqens,1);
mseensortzeit = mean(mseensort,2);
plotmse(:) = mseensortzeit(1,1,:);
figure
plot([1:N_ens],plotmse,'b-o');
title('quadratischer Fehler über Zeit und Ort gemittelt für verschiedene Mitglieder');
ylabel('MSE')
xlabel('Mitglied')

%-------------------------------------------------------------------------
% benötigte Zeit für einen Zeitschritt je Ensemble-Größe
% cpu_time2: Zeit für den Forecast und cpu_time3: Zeit für das Update
% 01-38.03-377.24
% 03-110.53-442.42
% 05-191.92-447.09 % 623.28 bzw. 414.57 bei einem Schritt der Länge 4.
% 10-422.20-493.14
% 15-575.63-501.44
% 20-794.16-551.27
% Zahl = [1,3,5,10,15,20];
% Zeitf = [38.03,110.53,191.92,422.20,575.63,794.16];
% Zeita = [377.24,442.42,447.09,493.14,501.44,551.27];
% figure
% plot(Zahl,Zeitf,'b-o')
% hold on
% plot(Zahl,Zeita,'g-*')
% title('benötigte Zeit für verschiedene Ensemble-Größen')
% xlabel('Größe')
% ylabel('Zeit (s)')
% legend('Forecast','Update')
%-------------------------------------------------------------------------

%save('workspace.mat')
%load('workspace.mat')

%save_fig.m
