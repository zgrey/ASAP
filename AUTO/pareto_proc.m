% Pareto processing routines
%% Load AS Matlab Environment
addpath '~/RESEARCH/active_subspaces/matlab/'
addpath '~/RESEARCH/active_subspaces/matlab/Subspaces'
addpath '~/RESEARCH/active_subspaces/matlab/ResponseSurfaces'
addpath '~/RESEARCH/active_subspaces/matlab/Plotters'
addpath '~/RESEARCH/active_subspaces/matlab/Domains'
addpath './paretopoints'
warning('off','MATLAB:dispatcher:pathWarning')

%% Regenerate inputs according to cst_samples.m
rng(45); % Used for 5000 Sample Designs
m = 10;
X = 2*rand(5000,m)-1;

%% Read Data
dir_str = dir('./forces');
[CL.F] = readsu2_forces(dir_str,2);
[CD.F] = readsu2_forces(dir_str,3);

% AoA = 2.5 data error fix 
% X    = [X(1:499,:);X(501:998,:);X(1000:4999,:)];
% CL.F = [CL.F(1:499);CL.F(501:998);CL.F(1000:4999)];
% CD.F = [CD.F(1:499);CD.F(501:998);CD.F(1000:4999)];

%% Compute drag subspace
AS.sub = compute(X,CD.F,[],[],6,0,0);
AS.X   = X;
if size(AS.sub.W1,2) < 2
    AS.sub.W1= [AS.sub.W1, AS.sub.W2(:,1)];
elseif size(AS.sub.W1,2) > 2
    AS.sub.W1 = AS.sub.W1(:,1:2);
end
Y      = AS.X*AS.sub.W1;

%% Compute pareto
[Pk,rowi] = prtp([CD.F,-CL.F]);
Yk = AS.X(rowi,:)*AS.sub.W1;

%% Find extremals
[~,CD.I] = min(CD.F);
[~,CL.I] = max(CL.F);
CD.y     = AS.X(CD.I,:)*AS.sub.W1;
CL.y     = AS.X(CL.I,:)*AS.sub.W1;

%% Lift summary
AS.F = CL.F;
% Eigenvectors
opts.markersize = 10;
eigenvectors(AS.sub.W1,opts);

% Subspace Errors
% subspace_errors(AS.sub.sub_br,opts);

% Plot eigenvalue decay
eigenvalues(AS.sub.eigenvalues,AS.sub.e_br,opts);
grid on

% Summary plot
opts.ylabel = 'QoI';
opts.markersize = 30;
sufficient_summary(AS.X*AS.sub.W1,AS.F,opts);
hold on; scatter(Yk(:,1),Yk(:,2),200,'k.');
scatter(CL.y(1), CL.y(2),100,'r','filled');
scatter(CD.y(1), CD.y(2),100,'r','linewidth',2.5);
grid off; axis equal;

%% Drag summary
AS.F = CD.F;
% Eigenvectors
opts.markersize = 10;
eigenvectors(AS.sub.W1,opts);

% Subspace Errors
% subspace_errors(AS.sub.sub_br,opts);

% Plot eigenvalue decay
eigenvalues(AS.sub.eigenvalues,AS.sub.e_br,opts);
grid on

% Summary plot
opts.ylabel = 'QoI';
opts.markersize = 30;
sufficient_summary(AS.X*AS.sub.W1,AS.F,opts);
hold on; scatter(Yk(:,1),Yk(:,2),200,'k.');
scatter(CD.y(1), CD.y(2),100,'r','filled');
scatter(CL.y(1), CL.y(2),100,'r','linewidth',2.5);
grid off; axis equal;

%% zonotope corners (pick one)

[~,zi] = min(Y(:,2));
yp     = Y(zi,:);

% [~,zi] = max(Y(:,1));
% yp     = Y(zi,:);

% [~,zi] = max(Y(:,2));
% yp     = Y(zi,:);

% [~,zi] = min(Y(:,1));
% yp     = Y(zi,:);

%% Parameterize path
t = linspace(0,1,50000);
p = [t'*yp + (1-t)'*CD.y; t'*CL.y + (1-t)'*yp];
plot(p(:,1),p(:,2),'r--','linewidth',3);

%% Approximate Pareto
CD.q = 3; CL.q = 2;
CD.coef = poly_train(AS.X*AS.sub.W1,CD.F,CD.q);
CL.coef = poly_train(AS.X*AS.sub.W1,CL.F,CL.q);
CD.P = poly_predict(p,CD.coef,CD.q);
CL.P = poly_predict(p,CL.coef,CL.q);

plot(CD.P,CL.P,'Linewidth',2,'color',[0,0.447,0.741]);
scatter(Pk(:,1),-Pk(:,2),'k.');
%% 
i = 1;
Y = X*CL.W(:,:,i); [~,zi] = min(Y(:,2)); yp     = Y(zi,:);
p = [t'*yp + (1-t)'*CD.y{i}; t'*CL.y{i} + (1-t)'*yp];
CD.cp{i} = poly_predict(p,poly_train(X*CD.W(:,:,i),CD.F{i},2),CD.q);
CL.cp{i} = poly_predict(p,poly_train(X*CD.W(:,:,i),CL.F{i},1),CL.q);
figure(1); hold on
% plot(CD.cp{i},CL.cp{i},'Linewidth',2,'color',[0,0.447,0.741]); 
Pk = P{i}; scatter(Pk(:,1),Pk(:,2),'k.');