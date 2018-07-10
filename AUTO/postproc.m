% Read Monte-Carlo sampling results
clc; clearvars;
%% Load AS Matlab Environment
% Assign your custom active_subspaces home directory
AS_HOME = '~/active_subspaces';
addpath([AS_HOME,'/matlab/']);
addpath([AS_HOME,'/matlab/Subspaces']);
addpath([AS_HOME,'/matlab/ResponseSurfaces']);
addpath([AS_HOME,'/matlab/Plotters']);
warning('off','MATLAB:dispatcher:pathWarning')

%% Read SU2 data files
% pick QoI
QOI = 3; % CL/CD QOI = 1, CL QOI = 2, CD QOI = 3 
% read force data
[F, N, I, Fmax, Fmin, Imax, Imin] = readsu2_forces(dir('./forces'),QOI);

%% Construct Active Subspace
disp('LOADING designs.mat...')
load('designs.mat');

% Check for consitency from designs.mat
if sum(IP) ~= N
    disp('WARNING: Total evaluations does not match initial number of feasible designs')
end

% Compute qphd active subspace
disp('Approximating active subspace...')
nboot = 500; sub = compute(X,F,[],[],6,0,nboot);
W1 = sub.W1;
if size(W1,2) > 2
    disp('NOTICE: Gap detected beyond two-dimensions. Reassigning eigenvectors for visualization...')
    % reassign one-dimensional subspace
%     W1 = sub.W1(:,1);
    % reassign two-dimensional subspace
    W1 = sub.W1(:,1:2);
end

%% Visualize
close all;
% Plot eigenvector(s) entries
if QOI == 1
    opts.title = 'C_l / C_d Subspace';
elseif QOI == 2
    opts.title = 'C_l Subspace';
elseif QOI == 3
    opts.title = 'C_d Subspace';
end
opts.markersize = 10;
eigenvectors(W1,opts);

% Plot subspace distance
subspace_errors(sub.sub_br);

% Plot eigenvalue decay
eigenvalues(sub.eigenvalues,sub.e_br); grid on;

% Summary plot
opts.markersize = 30;
sufficient_summary(X*W1,F,opts);
% Polynomial approximation and plot
% Order of polynomial
if QOI == 2
    p = 1;
elseif QOI == 3
    p = 2;
end
% Train polynomial and predict
[Coef, B, ~, ~, f_hat, r] = poly_train(X*W1, F, p);
NY = 50;

if size(W1,2) == 1
    % visualize one-dimensional approximation
    lb = min(X*W1); ub = max(X*W1);
    Y1 = linspace(lb,ub,NY)';
    [Fh,df] = poly_predict(Y1,Coef,p);
    figure(4); hold on; plot(Y1,Fh,'k--','linewidth',2);
elseif size(W1,2) == 2
    % visualize two-dimensional approximation
    % note, the following bounds extrapolate beyond the zonotope
    lb = [min(min(X*W1(:,1))), min(min(X*W1(:,2)))]; ub = [max(max(X*W1(:,1))), max(max(X*W1(:,2)))];
    [Y1,Y2] = meshgrid(linspace(lb(1),ub(1),NY),linspace(lb(2),ub(2),NY));
    [Fh,df] = poly_predict([reshape(Y1,NY*NY,1), reshape(Y2,NY*NY,1)],Coef,p);
    figure(5); hold on; contour(Y1,Y2,reshape(Fh,NY,NY),35,'--');
end


%% Save computations
save(['./subspaces/AS_m',num2str(m),'_N',num2str(N),'_pm',num2str(pct*100),'pct_','QOI',num2str(QOI),'.mat']);
