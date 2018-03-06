function [psia_e,K_e,Pa_e] = EnKF(psif_e1, d1, H, R)
% EnKF: Ensemble Kalman-Filter
%
%   psif_e:     Ensemble der Prediction (Forecast); 
%               Tiefendim: Ensemble; 2 Werte-Spalten; 
%               #Werte x #Modellgleichungen x #Ensemble-Mitglieder
%   d:          Observationen #Werte x #Modellgleichungen
%   H:          Observationsmatrix #Modellgleichungen x #Werte
%   R:          Observations-Covarianzmatrix
%   psia_e:     Ensemble des Updates; überlagerte Ensemblematrizen 
%               aus Werte-Spalten
%               #Werte x #Modellgleichungen x #Ensemble-Mitglieder
[sizeh, ~] = size(H);
psif_e(:,:)=psif_e1(:,1,1,:);
% Dimension der Observation (zu einer bestimmten Zeit)
N_obs = size(d1,1);

d=zeros(N_obs,1);
d(:)=d1(:,1,1); % hier ggf. syst. Fehler addieren +0.1

% Anzahl der Ensemble-Mitglieder, d.h. die Anzahl in der Tiefendimension von xf_e
N_ens = size(psif_e,2);
M = size(psif_e, 1);
% Covarianz-Matrix der Prediction (Forecast); Benennung nach eve03a
Pf_e = cov(psif_e');

% Observations-Ensemble
d_e = zeros(N_obs, N_ens);
for i = 1:N_ens
    d_e(:,i) = d+mvnrnd(zeros(1,N_obs), 1.0*R)';
    d_e(:,i) = max(min(d_e(:,i), 1), 0); % Clamp-Funktion
end

% Observations-Covarianzmatrix
R_e = 1.00*cov(d_e'); % hier ggf. verändern, z.B. 1000* oder +100*eye(709)
% Die Observationen sind die Zeilen, die Var. die Observationen

% Kalman gain
K_e = Pf_e*H'*pinv(H*Pf_e*H'+R_e); % (23) in eve03a
% pinv bringt beste Lösung im Vergleich zu "\", aber jenes bringt geringste
% nicht-null Komponenten

% Update-Ensemble erzeugen
psia_e = zeros(M,1,1,N_ens);

% Update-Schritt ausführen
for i = 1:N_ens
    psia_e(:,1,1,i) = psif_e(:,i)+K_e*(d_e(:,i)-H*psif_e(:,i));
end
% (20) mit (23) eingesetzt in eve03a

% Update-Covarianz
Pa_e = (eye(sizeh)-K_e*H)*Pf_e;

end