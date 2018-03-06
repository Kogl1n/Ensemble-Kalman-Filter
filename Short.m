function [model, u] = Short(rho0, B0, eta, A0, B_asterisk, tliststart, tlistdelta, tlistend)
%% Programm zur Modellierung des Systems von Short

% 2 Gleichungen
model=createpde(2);
% Gebiet erstellen 
% Als grundlegende Geometrie wurde immer ein Quadrat mit Ecken (0,0) und (1
% 0,10) gewählt, genannt rect1.

% für Quadrat: in erster Spalte 3 und 4.
% x-Koord.: Spalten 3-6
% y-Koord.: Spalten 7-10
rect1 = [3; 4; 0; 10; 10; 0; 0; 0; 10; 10];
%rect2 = [3; 4; 2.8; 3; 3; 2.8; 2.8; 2.8; 3; 3;];
%rect2 = [3; 4; 2; 3; 3; 2; 2; 2; 3; 3];
%rect2 = [3; 4; 5; 10; 10; 5; 10; 10; 15; 15;];
% decsg benötigt Stringnamen
% Vereinigung und Schnitt
%gd = [rect1, rect2];
%ns = char('rect1', 'rect2');
%ns = ns';
%sf = 'rect1+rect2';
%[dl, bt] = decsg(gd, sf, ns);
%[dl2, bt2] = csgdel(dl, bt);

% Geometrie für rect1 erzeugen:
dl = decsg(rect1);
% Geometrie ins Modell
geometryFromEdges(model, dl);
% Plot der Geometrie
% pdegplot(dl, 'EdgeLabels', 'on', 'subdomainLabels', 'on');
% Randbedingungen einsetzen, 2D: 'edge'.
% q, g Neumann RB, falls fehlend autom. noflux
Q0 = zeros(2);
G0 = zeros(2,1);
G1= [1; 1];
Q1 = [1 0; 0 1];
H0 = [0 0; 0 1];
R0 = [0, 1];
HM = [0 1];
% Randbedingungen:
applyBoundaryCondition(model, 'edge', [1, 2, 3, 4], 'q' ,Q0, 'g', G0);
%applyBoundaryCondition(model, 'edge', 3, 'q', Q0, 'g', -G0);
%applyBoundaryCondition(model, 'edge', 1, 'q', H0, 'g', G1);
applyBoundaryCondition(model, 'edge', 1, 'u', 1, 'EquationIndex', 1);
%applyBoundaryCondition(model, 'edge', 3, 'u', -1, 'EquationIndex', 2);
%applyBoundaryCondition(model, 'edge', [1, 2, 3, 4], 'q', Q0, 'g', G0);
%applyBoundaryCondition(model, 'edge', 1, 'h', H0, 'r', R0) ;
%applyBoundaryCondition(model, 'edge', [4], 'h', H0, 'r', [0 0]);
%applyBoundaryCondition(model, 'edge', [1, 2, 3, 4], 'q', Q0, 'g', G0);
%applyBoundaryCondition(model, 'edge', [1, 2, 3, 4, 5, 6, 7,], 'q',Q0, 'g', G0);

% Setze globale Variablen für den Gebrauch in PDE-Parameter-Funktionen
global rho0_short;
global B0_short; 
global eta_short;
global A0_short;
global B_asterisk_short;

rho0_short = rho0;
B0_short = B0; 
eta_short = eta;
A0_short = A0;
B_asterisk_short = B_asterisk;

% PDE-Koeffizienten festlegen; parabolisch =>Faktor bei 2.zeitl. Abl. 0;
% Skalierung 1. zeitl. Abl. 1

specifyCoefficients(model, 'm', 0, 'd', 1, 'c', @ccoeffunction, 'a',...
    @acoeffunction, 'f', @fcoeffunction);
% Anfangsbedingungen: function handle u0fun
u0 = @u0fun;
setInitialConditions(model, u0);
% Gitter für Numerik erzeugen
% Die maximalen Abstände zwischen Knotenpunkten können mit Hmax gesteuert werden.
generateMesh(model, 'Hmax', 0.8); %0.5

t_num = round((tlistend-tliststart)/tlistdelta); % Anzahl der Zeitschritte
tlistend = tliststart+t_num*tlistdelta;
% Zeitschritte angegeben
tlist = tliststart:tlistdelta:tlistend;
% Lösen: autom. nichtlin. Löser, falls nichtlin. System
result = solvepde(model, tlist);
% Abspeichern von der Lösung u sowie den Gradienten.
u = result.NodalSolution;
ux = result.XGradients;
uy = result.YGradients;

% Lösungen für alle Zeitschritte plotten.
% for i=1:length(tlist)
%     figure
%     %set(gca, 'CLim', [0 1]); % Colorbars harmonisieren
%     subplot(2,1,1) % Attraktivität
%     pdeplot(model, 'xydata', u(:,1,length(tlist)-i+1), 'contour', 'off', 'colormap', 'jet');
%     %pdeplot(model, 'xydata', u(:,1,length(tlist)-i+1), 'flowdata', ...
%     %[ux(:,1,length(tlist)-i+1), uy(:,1,length(tlist)-i+1)], 'contour', 'off', 'colormap', 'cool');
%     title(['Attraktivität zu T=',num2str(tlistend-(i-1)*tlistdelta)])
%     xlabel('x-Koordinate')
%     ylabel('y-Koordinate')
%     axis equal
%     caxis([0, 1]);
%     
%     subplot(2,1,2) % Einbrecher
%     pdeplot(model, 'xydata', u(:,2,length(tlist)-i+1), 'contour', 'off', 'colormap', 'jet');
%     %pdeplot(model, 'xydata', u(:,2,length(tlist)-i+1), 'flowdata', ...
%     %[ux(:,2,length(tlist)-i+1), uy(:,2,length(tlist)-i+1)], 'contour', 'off', 'colormap', 'cool');
%     title(['Verbrecherdichte zu T=',num2str(tlistend-(i-1)*tlistdelta)])
%     xlabel('x-Koordinate')
%     ylabel('y-Koordinate')
%     axis equal
%     caxis([0, 1]);
% end

end

