% Pareto processing routines
clc; clearvars;
%% Load AS Matlab Environment
% Assign your custom active_subspaces home directory
AS_HOME = '~/active_subspaces';
addpath([AS_HOME,'/matlab/']);
addpath([AS_HOME,'/matlab/Subspaces']);
addpath([AS_HOME,'/matlab/ResponseSurfaces']);
addpath([AS_HOME,'/active_subspaces/matlab/Plotters']);
warning('off','MATLAB:dispatcher:pathWarning');
addpath './paretopoints';

%% Load designs workspace
disp('LOADING designs.mat...')
load('designs.mat');
% Check for consitency from designs.mat
if sum(IP) ~= N
    disp('WARNING: Total evaluations does not match initial number of feasible designs')
end

%% Read Data
dir_str = dir('./forces');
[CL.F] = readsu2_forces(dir_str,2);
[CD.F] = readsu2_forces(dir_str,3);

%% Compute drag subspace
AS.sub = compute(X,CD.F,[],[],6,0,0);
AS.X   = X;
% Assign two dimensional subspace
if size(AS.sub.W1,2) < 2
    AS.sub.W1 = [AS.sub.W1, AS.sub.W2(:,1)];
elseif size(AS.sub.W1,2) > 2
    AS.sub.W1 = AS.sub.W1(:,1:2);
end
% Project to active subspace
Y      = AS.X*AS.sub.W1;
% Assign subspace
CD.W   = AS.sub.W1; CL.W = AS.sub.W1;

%% Compute pareto
[Pk,rowi] = prtp([CD.F,-CL.F]);
Yk = AS.X(rowi,:)*AS.sub.W1;

%% Find extremals
[~,CD.I] = min(CD.F); [~,CL.I] = max(CL.F);
CD.y = AS.X(CD.I,:)*AS.sub.W1; CL.y = AS.X(CL.I,:)*AS.sub.W1;

%% Lift summary
AS.F = CL.F;
% Summary plot
opts.title = 'Lift Summary'; opts.ylabel = 'QoI'; opts.markersize = 30;
sufficient_summary(AS.X*AS.sub.W1,AS.F,opts);
hold on; scatter(Yk(:,1),Yk(:,2),200,'k.');
scatter(CL.y(1), CL.y(2),100,'r','filled');
scatter(CD.y(1), CD.y(2),100,'r','linewidth',2.5);
grid off; axis equal;

%% Drag summary
AS.F = CD.F;
% Summary plot
opts.title = 'Drag Summary'; opts.ylabel = 'QoI'; opts.markersize = 30;
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
figure;
plot(CD.P,CL.P,'Linewidth',2,'color',[0,0.447,0.741]); hold on;
scatter(Pk(:,1),-Pk(:,2),'k.');