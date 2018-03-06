function [u_ens] = EnKF_ens_create(u_gen0, n_ens, sigma)
%% EnKF_ens_create: Ensemble erzeugen
%   u_gen0: Anfangswert
%   sigma:  Standardabweichung der Normalverteilung N(0,sigma^2)
%   n_ens:  Anzahl der Ensemblemitglieder
%   u_ens:  Ensemble; 3D-Array mit in der Tiefe positionierten
%           Ensemble-Matrizen aus Anfangswertspalten des Modells

[dim_val, dim_mod] = size(u_gen0);
u_ens = zeros(dim_val,dim_mod,1,n_ens);
if n_ens==1
    u_ens(:,:,1,1) = u_gen0;
    u_ens(:,:,1,1) = max(min(u_ens(:,:,1), 1), 0); % Clamp-Funktion
    return 
end
for i = 1:n_ens
    u_ens(:,:,1,i) = u_gen0+normrnd(0, sigma, [dim_val,dim_mod]); 
    % normrnd unterstützt leider kein Festlegen der kryptograph. Seed
    % normrnd, nicht mvnrnd
    
    % Alternative: multiplikative statt additiver Abweichung
    u_ens(:,:,1,i) = max(min(u_ens(:,:,i), 1), 0); % Clamp-Funktion
end

end