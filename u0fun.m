function uinit = u0fun(locations)
%% Anfangswerte festlegen
% freie Parameter festlegen
global rho0_short;
global B0_short;

% uinit
M = length(locations.x); % Dies ist identisch mit der Anzahl der Gesamtmeshpunkte,
% d.h. size(model.Mesh.Nodes,2)
uinit = zeros(2,M);
% genaue Werte für uinit:

uinit(1,:) = B0_short;%+normrnd(0, 0.1, [1,M]);%+normrnd(0, 0.1, [1,M]) %+0.01*rand(1,M);
uinit(2,:) = rho0_short;%+normrnd(0, 0.1, [1,M]);%normrnd(0, 0.1, [1,M]) %+0.01*rand(1,M);
%uinit(1,:) = locations.y*0.05+0.25;
%uinit(2,:) = 0.018*locations.y+0.01;

end