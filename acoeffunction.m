function amatrix = acoeffunction(region, state)
%% Koeffizienten von A festlegen:
% freie Parameter festlegen
s = 2;
global A0_short;

% Größe definieren.
% state.u(1,:) ist die Auswertung von u1.
n1 = 4;
M = numel(region.x); % Anzahl der Gesamtmeshpunkte
% d.h. size(model.Mesh.Nodes,2)
amatrix = zeros(n1,M);
amatrix(1,:) = +s*ones(1,M);
amatrix(2,:) = zeros(1,M);
amatrix(3,:) = s*(-state.u(1,:)-A0_short);
amatrix(4,:) = s*(state.u(1,:)+A0_short);
%amatrix(3,:) = s*(-state.u(1,:)-0.018*region.y+0.01);
%amatrix(4,:) = s*(state.u(1,:)+0.018*region.y+0.01);

end
