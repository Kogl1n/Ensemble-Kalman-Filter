function f = fcoeffunction(region, state)
%% Koeffizienten von f festlegen:
% freie Parameter festlegen

global B_asterisk_short;
N = 2; % Anzahl der Gleichungen
M = length(region.x); % Anzahl der Gesamtmeshpunkte
% d.h. size(model.Mesh.Nodes,2)
f = zeros(N,M); % Erstellung von f
s = 2;

% Hier wird f genau angegeben.
f(1,:) = 0;
f(2,:) = s*B_asterisk_short;

end
