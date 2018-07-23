%% A wrapper script for generating CST airfoil samples and SU2 meshes
clc; close all; clearvars;

%% Discretize airfoil domain
np = 250; l = [linspace(0,0.005,50),linspace(0.005+.005/50,1,np)];

%% Generate Example Nominal (fit to NACA 0012)
% number of parameters (even number)
m = 10;
% dv(1:m/2) correspond to lower surface coefficients, dv(m/2+1:m) are upper surface coefficients
% optimize to find nominal values
dv = fminunc(@(dv) coord_obj(l,dv), [-0.2*ones(1,m/2), 0.2*ones(1,m/2)]);
% compute coordinates for optimal coefficients
[coordU0, coordL0] = cst_airfoil(l',dv(1:m/2)',dv(m/2+1:end)',0);
% plot Nominal Airfoil
% plot(l,coordU0(2,:),'LineWidth',2,'color',0.75*ones(3,1)); hold on;
% plot(l,coordL0(2,:),'LineWidth',2,'color',0.75*ones(3,1)); axis equal;
% title 'Example Nominal Airfoil Shape'

%% Random sweeps
% number of sweeps
N = 10;
% number of points over sweep
T = 50; t = linspace(0,1,T);

% Define upper and lower bounds
pct = 0.6;
% lower surface
lb0(1:m/2) = (1+pct)*dv(1:m/2); ub0(1:m/2) = (1-pct)*dv(1:m/2);
% upper surface
ub0(m/2+1:m) = (1+pct)*dv(m/2+1:end); lb0(m/2+1:m) = (1-pct)*dv(m/2+1:end);

% random samples from [-1,1] hypercube
X = 2*rand(2*N,m) - 1;
% normalize by sup-norm to scale to boundary
X = X./repmat(max(abs(X),[],2),1,m);
% build convex combinations
X  = kron(X(1:N,:)',(1-t))' + kron(X(N+1:end,:)',t)';
X0 = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(X+1)));

% initialize admissible set count and indices
NF = 0; IP = zeros(N*T,1);
% initialize thickness metric
max_thk = zeros(N*T,1); I_maxthk = max_thk; L1 = zeros(N*T,1);
for i=1:N*T
    % Evaluate CST
    [coordU, coordL] = cst_airfoil(l',X0(i,(1:m/2)),X0(i,(m/2+1:end)),0);
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
for i=1:( (N*T)-NF )
    
    % Evaluate CST
    [coordU, coordL] = cst_airfoil(l',X0(i,(1:length(dv)/2)),X0(i,(length(dv)/2+1:end)),0);
    % run meshing routines
    mesh_coords(coordU,coordL,i,N*T,NF);
    
end
%% Save sweeps workspace
save(['./sweeps/','sweeps_m',num2str(m),'_N',num2str(N*T),'_pm',num2str(pct*100),'pct.mat'],'X','X0','NF','IP','lb0','ub0','l','max_thk','I_maxthk','L1','N','T','t','pct');
