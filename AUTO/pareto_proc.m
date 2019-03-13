% Pareto processing routines
clc; clearvars;
%% Load AS Matlab Environment
% Assign your custom active_subspaces home directory
AS_HOME = '~/active_subspaces';
addpath([AS_HOME,'/matlab/']);
addpath([AS_HOME,'/matlab/Subspaces']);
addpath([AS_HOME,'/matlab/ResponseSurfaces']);
addpath([AS_HOME,'/matlab/Plotters']);
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

%% Compute lift & drag subspace
CD.sub = compute(X,CD.F,[],[],6,0,0);
CL.sub = compute(X,CL.F,[],[],6,0,0);
CD.X   = X;
% Assign two dimensional subspaces
if size(CD.sub.W1,2) < 2
    CD.sub.W1 = [CD.sub.W1, CD.sub.W2(:,1)];
elseif size(CD.sub.W1,2) > 2
    CD.sub.W1 = CD.sub.W1(:,1:2);
end
% Assign two dimensional subspace
if size(CL.sub.W1,2) < 2
    CL.sub.W1 = [CL.sub.W1, CL.sub.W2(:,1)];
elseif size(CL.sub.W1,2) > 2
    CL.sub.W1 = CL.sub.W1(:,1:2);
end
% Project to drag active subspace
Y      = CD.X*CD.sub.W1;
% Assign subspace
CD.W   = CD.sub.W1; CL.W = CD.sub.W1;

%% Compute pareto
[Pk,rowi] = prtp([CD.F,-CL.F]);
Yk = CD.X(rowi,:)*CD.sub.W1;

%% Find extremals
[~,CD.I] = min(CD.F); [~,CL.I] = max(CL.F);
CD.y = CD.X(CD.I,:)*CD.sub.W1; CL.y = CD.X(CL.I,:)*CD.sub.W1;

% %% Lift summary
% CD.F = CL.F;
% % Summary plot
% opts.title = 'Lift Summary'; opts.ylabel = 'QoI'; opts.markersize = 30;
% sufficient_summary(CD.X*CD.sub.W1,CD.F,opts);
% hold on; scatter(Yk(:,1),Yk(:,2),200,'k.');
% scatter(CL.y(1), CL.y(2),100,'r','filled');
% scatter(CD.y(1), CD.y(2),100,'r','linewidth',2.5);
% grid off; axis equal;
% 
% %% Drag summary
% CD.F = CD.F;
% % Summary plot
% opts.title = 'Drag Summary'; opts.ylabel = 'QoI'; opts.markersize = 30;
% sufficient_summary(CD.X*CD.sub.W1,CD.F,opts);
% hold on; scatter(Yk(:,1),Yk(:,2),200,'k.');
% scatter(CD.y(1), CD.y(2),100,'r','filled');
% scatter(CL.y(1), CL.y(2),100,'r','linewidth',2.5);
% grid off; axis equal;

%% Mix two-dimensional subspaces over a Grassmann geodesic
% "Mix" the subspaces over the Grassmann Manifold from lift to drag AS
% compute projection onto the orthogonal complement of lift AS
CL.Proj = eye(m) - CL.sub.W1*CL.sub.W1';
% compute the reduced (thin) SVD
[U,S,V] = svd(CL.Proj*CD.sub.W1*inv(CL.sub.W1'*CD.sub.W1),0);
% parameterize the exponential map between subspaces
Exp = @(t) CL.sub.W1*V*diag(cos(atan(diag(S))*t)) + U*diag(sin(atan(diag(S))*t));

% make a gif of the changing subspaces
gifname = ['./subspaces/grassmann_m',num2str(m),'_N',num2str(N),'_pm',num2str(pct*100),'pct','.gif']; fig = figure;
T = 100; t = linspace(0,1,T);
for ii=1:T
    W = Exp(t(ii));
    if size(W,2) == 2
        % plot lift over first two Grassmann eigenvectors
        subplot(1,2,1), hold on; h1 = scatter(X*W(:,1),X*W(:,2),100,'filled','cdata',CL.F);
        title 'C_l'; axis equal; colorbar; fig.CurrentAxes.DataAspectRatio = [1 1 3];
        % plot drag over first two Grassmann eigenvectors
        subplot(1,2,2), hold on; h2 = scatter(X*W(:,1),X*W(:,2),100,'filled','cdata',CD.F);
        title 'C_d'; axis equal; colorbar; fig.CurrentAxes.DataAspectRatio = [1 1 3];
        % project the pareto points onto the Grassmann geodesic subspace
        subplot(1,2,1), hold on, h3 = scatter(X(rowi,:)*W(:,1),X(rowi,:)*W(:,2),200,'k.');
        subplot(1,2,2), hold on, h4 = scatter(X(rowi,:)*W(:,1),X(rowi,:)*W(:,2),200,'k.');
    elseif size(W,2) == 1
        subplot(1,2,1), h1 = scatter(X*W,CL.F,50,'filled','MarkerEdgeColor','k');
        ylabel 'C_l'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
        subplot(1,2,2), h2 = scatter(X*W,CD.F,50,'filled','MarkerEdgeColor','k');
        ylabel 'C_d'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    end
    % build gif
    figure(fig); frame = getframe(fig); 
    [A,map] = rgb2ind(frame2im(frame),256);
    if ii == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.1);
    end
    % delete plots
    delete([h1,h2,h3,h4]);
end
% replot an intermediate subspace
W = Exp(0.5);
if size(W,2) == 2
    % plot lift over first two Grassmann eigenvectors
    subplot(1,2,1), hold on; h1 = scatter(X*W(:,1),X*W(:,2),100,'filled','cdata',CL.F);
    title 'C_l'; axis equal; colorbar; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    % plot drag over first two Grassmann eigenvectors
    subplot(1,2,2), hold on; h2 = scatter(X*W(:,1),X*W(:,2),100,'filled','cdata',CD.F);
    title 'C_d'; axis equal; colorbar; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    % project the pareto points onto the Grassmann geodesic subspace
    subplot(1,2,1), hold on, h3 = scatter(X(rowi,:)*W(:,1),X(rowi,:)*W(:,2),200,'k.');
    subplot(1,2,2), hold on, h4 = scatter(X(rowi,:)*W(:,1),X(rowi,:)*W(:,2),200,'k.');
elseif size(W,2) == 1
    subplot(1,2,1), h1 = scatter(X*W,CL.F,50,'filled','MarkerEdgeColor','k');
    ylabel 'C_l'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    subplot(1,2,2), h2 = scatter(X*W,CD.F,50,'filled','MarkerEdgeColor','k');
    ylabel 'C_d'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
end

% %% Approximate Pareto
% CD.q = 3; CL.q = 2;
% CD.coef = poly_train(CD.X*CD.sub.W1,CD.F,CD.q);
% CL.coef = poly_train(CD.X*CD.sub.W1,CL.F,CL.q);
% CD.P = poly_predict(p,CD.coef,CD.q);
% CL.P = poly_predict(p,CL.coef,CL.q);
% figure;
% plot(CD.P,CL.P,'Linewidth',2,'color',[0,0.447,0.741]); hold on;
% scatter(Pk(:,1),-Pk(:,2),'k.');