
solveropts.returnDiv = true;      %return an E field divergence operator for bound charge calculation below
inspect_only = false;
addpath(genpath('~/Documents/MATLAB/maxwellfdfd-master'));

%resolution
dl = 0.01;
c0= 1;

%geometric parameters
L = 1;              %box width
src = 0.25*L;       %source distance to box edge

%frequency range
Nf = 100;                        %number of freqs
freqs = linspace(0.1,.3,Nf);    %frequency range (units of c/L)
freqs = 1;
Nf = length(freqs);

%simulation domain
X = L/2+src+L+dl*10;    %X half width
Y = X;                  %Y half width
    
%box material
eps = 16%11.56;            %box permittivity


%arrays to store sensitivities
Is = zeros(Nf,1);

for fi = (1:Nf)
    B = Box([-L/2 L/2; -L/2 L/2; 0 1]);
    
    perc_done = fi/Nf*100

    freq = freqs(fi);
    lam = c0/freq;
    omega = 2*pi*freq;
    k0 = omega/c0;
    % original calculation
    [E, H, obj_array, src_array, extra, epsTyler] = maxwell_run(...
        'OSC', 1, lam, ...
        'DOM', {'vacuum', 'none', 1.0}, [-X X; -Y Y; 0 1], [dl,dl,1] , BC.p, [dl*10 dl*10 0], ...
        'OBJ', ...
            {'box','k',eps}, ...
                B,...                            
        'SRCM', PointSrc(Axis.z,[-src-L/2,0,0],1), ...
        solveropts, inspect_only);

    % calculate H_z at observation point
    Ex = E{Axis.x}.data_original;
    Ey = E{Axis.y}.data_original;
    Hz = H{Axis.z}.data_original;

    [NHx, NHy] = size(Hz);    
    H_obs = Hz(floor(NHx/2+(L/2+src)/dl)+1,floor(NHy/2)+1);
    I = abs(H_obs)^2;
    Is(fi) = I;
    
end
figure(2); clf;
imagesc(real(Hz));

%% Plot results
figure(1)
clf; hold all;

plot(freqs,Is);
xlabel('frequency (c/L)')
ylabel('normalized sensitivity (a.u.)')
legend('Direct','AVM')
title('Hz : 2d : varying device length')
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','normal')
set(gca,'FontSize',25,'fontWeight','normal')
grid