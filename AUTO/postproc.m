% Read SHDP Results
clc; clearvars; close all;
addpath '~/RESEARCH/active_subspaces/matlab/'
addpath '~/RESEARCH/active_subspaces/matlab/Subspaces'
addpath '~/RESEARCH/active_subspaces/matlab/ResponseSurfaces'
addpath '~/RESEARCH/active_subspaces/matlab/Plotters'
warning('off','MATLAB:dispatcher:pathWarning')

%% Read SU2 data files
% pick QoI
QOI = 2; % QOI = 1: CL/CD, QOI = 2: CL, QOI = 3: CD 
% read force data
[F, N, I, Fmax, Fmin, Imax, Imin] = readsu2_forces(dir('./forces'),QOI);

%% Construct Active Subspace
disp('LOADING designs.mat...')
if exist('designs.mat','file')
    load('designs.mat');
else
    disp('ERROR: designs.mat was not found')
    return
end

% filter designs if necessary by FU
% FU = 0.1; ind = (F < FU);
% X = X(ind,:); F = F(ind);

% Check for consitency from designs.mat
if sum(IP) ~= N
    disp('WARNING: Total evaluations does not match initial number of feasible designs')
end

% Compute qphd active subspace
disp('Approximating active subspace...')
nboot = 500; sub = compute(X,F,[],[],6,0,nboot);
W1 = sub.W1;
if size(W1,2) > 2
    disp('WARNING: Gap not detected, reassigning eigenvectors')
    W1 = sub.W1(:,1:2);
end
e  = sub.eigenvalues;

%% Plot
% Eigenvector(s)
if QOI == 1
    opts.title = 'C_l / C_d Subspace';
elseif QOI == 2
    opts.title = 'C_l Subspace';
elseif QOI == 3
    opts.title = 'C_d Subspace';
end
opts.markersize = 10;
eigenvectors(W1,opts);
% print('-f1','./figs/eig_vec','-djpeg')

% Subspace Errors
subspace_errors(sub.sub_br);
% print('-f2','./figs/sub_err','-djpeg')

% Plot eigenvalue decay
eigenvalues(e,sub.e_br);
grid on
% print('-f3','./figs/eig_val','-djpeg')

% Summary plot
opts.markersize = 30;
sufficient_summary(X*W1,F,opts);
% print('-f4','./figs/summary','-djpeg')

if size(W1,2) == 2
%     print('-f5','./figs/summary2','-djpeg')
end

% Approximate g() for 2D Case
% if size(W1,2) == 2
%     p = 2;
% 
%     [Coef, B, ~, ~, f_hat, r] = poly_train(X*W1, F, p);
%     NY = 50;
%     lb = [min(min(X*W1(:,1))), min(min(X*W1(:,2)))];
%     ub = [max(max(X*W1(:,1))), max(max(X*W1(:,2)))];
%     [Y1,Y2] = meshgrid(linspace(lb(1),ub(1),NY),linspace(lb(2),ub(2),NY));
%     [Fh,df] = poly_predict([reshape(Y1,NY*NY,1), reshape(Y2,NY*NY,1)],Coef,p);
%     contour(Y1,Y2,reshape(Fh,NY,NY),35);
%     hold on;
%     print('-f4','./figs/summary2','-djpeg')
% end

