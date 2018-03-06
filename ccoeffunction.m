function cmatrix = ccoeffunction(region, state)
%% Koeffizienten von C festlegen:
global A0_short;
global eta_short;


%eta = (3*B/(A0+B)+1-sqrt(12*B/(A0+B)))/(A0+B)-0.01;
%eta = (-0.002+3*B/(A0+B)+1-sqrt(12*B/(A0+B)))/(A0+B);

% Da Matlab mit einem Tensor bei der Multiplikation mit C arbeitet, wird
% eine 4x4 Matrix benötigt, bzw. ein 16x1 Vektor. Zusätzlich wird
% eine weitere Dimension für jeden Auswertungspunkt benötigt. Diese wird
% durch den Befehl numel(region.x) erzeugt und angehängt.
n1 = 16;
M = numel(region.x); % Anzahl der Gesamtmeshpunkte
% d.h. size(model.Mesh.Nodes,2)
cmatrix = zeros(n1, M);
cmatrix(1,:) = eta_short*ones(1, M);
%cmatrix(1,:) = eta*round(abs((region.subdomain-2)))+0.001;
cmatrix(2,:) = 0*ones(1,M);
cmatrix(3,:) = cmatrix(2,:);
cmatrix(4,:)= eta_short*ones(1,M);
%cmatrix(4,:) = eta*round(abs((region.subdomain-2)))+0.001;
cmatrix(5,:) = -2*state.u(2,:)./(state.u(1,:)+A0_short);
%cmatrix(5,:) = -2*state.u(2,:)./(state.u(1,:)+0.018*region.y+0.01);
cmatrix(6,:) = cmatrix(2,:);
cmatrix(7,:) = cmatrix(2,:);
%cmatrix(8,:) = -2*state.u(2,:)./(state.u(1,:)+0.018*region.y+0.01);
cmatrix(8,:) = -2*state.u(2,:)./(state.u(1,:)+A0_short);
cmatrix(9,:) = cmatrix(2,:);
cmatrix(10,:) = cmatrix(2,:);
cmatrix(11,:) = cmatrix(2,:);
cmatrix(12,:) = cmatrix(2,:);
cmatrix(13,:) = ones(1,M);
cmatrix(14,:) = cmatrix(2,:);
cmatrix(15,:) = cmatrix(2,:);
cmatrix(16,:) = ones(1,M);

end