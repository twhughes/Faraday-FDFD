clear all;
solveropts.returnDiv = true;      %return an E field divergence operator for bound charge calculation below
inspect_only = false;

num_trials = 10;
ns = (100:20:800);
N = length(ns);
times = zeros(N,num_trials);
lam = 50;

for j = (1:num_trials)
for i = (1:N)
    tic();
    n = floor(ns(i));
    L1 = floor(n/3); 
    L2 = floor(2*n/3);
    B = Box([L1 L2; L1 L2; 0 1]);
    
    X = floor(n);    
    [E, H, obj_array, src_array, extra, epsTyler] = maxwell_run(...
        'OSC', 1, lam, ...
        'DOM', {'vacuum', 'none', 1.0}, [0 X; 0 X; 0 1], [1,1,1] , BC.p, [20 20 0], ...
        'SOBJ', ...
            {'box','k',3}, ...
                B,...                            
        'SRCM', TFSFPlaneSrc([40 n-40; 40 n-40; 0 1], Axis.x, Axis.y), ...
        solveropts, inspect_only);
    
    % calculate H_z at observation point
    Hz = H{Axis.z}.data_original;
    Ex = E{Axis.x}.data_original;
    Ey = E{Axis.y}.data_original;
    time = toc();
    times(i,j) = time;
end
end
%%
figure(1); clf; hold all;
plot(ns.*ns./1000,times);
plot(ns_f.*ns_f./1000,times_f);
set(findall(gcf,'type','text'),'FontSize',25,'fontWeight','normal');
set(gca,'FontSize',25,'fontWeight','normal');
xlabel('number of grid points (thousands)');
ylabel('simulation time (seconds)');
legend('Maxwell-FDFD','Faraday-FDFD');
