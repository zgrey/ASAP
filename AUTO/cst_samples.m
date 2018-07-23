%% A wrapper script for generating CST airfoil samples and SU2 meshes
close all; clearvars;

%% discretize airfoil domain
nl = 250; l = [linspace(0,0.005,50),linspace(0.005+.005/50,1,nl)];

%% Generate Example Nominal (fit to NACA 0012)
% number of parameters (even number)
m = 10;
% dv(1:m/2) correspond to lower surface coefficients, dv(m/2+1:m) are upper surface coefficients
% optimize to find nominal values
dv = fminunc(@(dv) coord_obj(l,dv), [-0.2*ones(1,m/2), 0.2*ones(1,m/2)]);
% compute coordinates for optimal coefficients
[coordU0, coordL0] = cst_airfoil(l',dv(1:m/2)',dv(m/2+1:end)',0);
% plot Nominal Airfoil
% plot(l,coordU0(2,:),'b','LineWidth',2,'color',0.75*ones(3,1)); hold on;
% plot(l,coordL0(2,:),'b','LineWidth',2,'color',0.75*ones(3,1)); axis equal;
% title 'Example Nominal Airfoil Shape'

%% Random Perturbation from Hypercube
% Sample uniform hypercube:
N = 1000; rng(47);

% Define upper and lower bounds
pct = 0.2;
% lower surface
lb0(1:m/2) = (1+pct)*dv(1:m/2)'; ub0(1:m/2) = (1-pct)*dv(1:m/2);
% upper surface
ub0(m/2+1:m) = (1+pct)*dv(m/2+1:end)'; lb0(m/2+1:m) = (1-pct)*dv(m/2+1:end);

% print bounds
disp('----- Parameter Bounds -----')
disp(['Param.', '   ','lb','      ', 'ub']);
disp([['a1 ';'a2 ';'a3 ';'a4 ';'a5 ';'a6 ';'a7 ';'a8 ';'a9 ';'a10',], ...
    [' ';' ';' ';' ';' ';' ';' ';' ';' ';' '],num2str(lb0'), ...
    [' ';' ';' ';' ';' ';' ';' ';' ';' ';' '],num2str(ub0')]);

X = 2*rand(N,length(dv))-1;
X0 = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(X+1)));

% initialize admissible set count and indices
NF = 0; IP = zeros(N,1);
% initialize thickness metric
max_thk = zeros(N,1); I_maxthk = max_thk; L1 = zeros(N,1);
for i=1:N
    % Evaluate CST
    [coordU, coordL] = cst_airfoil(l',X0(i,(1:length(dv)/2)),X0(i,(length(dv)/2+1:end)),0);
    % common engineering airfoil metrics
    [max_thk(i), I_maxthk(i)] = max(coordU(2,:)-coordL(2,:));
    L1(i) = mean(abs(coordU(2,:) + coordL(2,:)));
    % check feasibility
    if min(coordU(2,:)-coordL(2,:)) >= -1e-3
        IP(i) = 1;
        %  verify perturbations with plots
%         h1 = plot(l,coordU(2,:),'k--','LineWidth',2); h2 = plot(l,coordL(2,:),'k--','LineWidth',2);
%         title(['Example Nominal Airfoil Shape',' & Pert i = ',num2str(i)])
%         pause(0.25); delete([h1,h2]);
    else
        NF = NF + 1;
        scatter([l,l],[coordU(2,:), coordL(2,:)],'r.'); axis equal;
        pause(0.25)
    end
end
IP = logical(IP);
disp(['Number of feasible designs: ',num2str(sum(IP))])

% Reassign indices according to feasible airfoils
X  =  X(IP,:); X0 =  X0(IP,:);

%% Generate Meshes
for i=1:N-NF
    
    % Evaluate CST
    [coordU, coordL] = cst_airfoil(l',X0(i,(1:length(dv)/2)),X0(i,(length(dv)/2+1:end)),0);    
    mesh_coords(coordU,coordL,i,N,NF);
    
end
%% Save design workspace
save(['./designs/','designs_m',num2str(m),'_N',num2str(N),'_pm',num2str(pct*100),'pct.mat'],'X','X0','NF','IP','lb0','ub0','l','max_thk','I_maxthk','L1','pct','N','m');
