%% post-process sweeps
clc; clearvars; close all;

%% Read assocaited workspace generated with cst_sweeps.m
disp('LOADING sweeps.mat...')
if exist('sweeps.mat','file')
    load('sweeps.mat');
else
    disp('ERROR: sweeps.mat was not found')
    return
end
m = size(X,2);

%% Read SU2 data files
% pick QoI
QOI = 2; % QOI = 1: CL/CD, QOI = 2: CL, QOI = 3: CD 
% read force data
F = readsu2_forces(dir('./forces'),QOI);

%% Plot and visualize sweeps
% more points over sweep for geometry visualization
TT = 25; tt = linspace(0,1,TT);

% Index set (columns are sweep indices, rows are point indices in sweep)
I = reshape(linspace(1,N*T,N*T)',T,N);
fig1 = figure; hold on; fig2 = figure; hold on; fig3 = figure; hold on; fig4 = figure; hold on;  
filename = ['./sweeps/sweeps_m',num2str(m),'_N',num2str(N*T),'_pm',num2str(pct*100),'pct_','QOI',num2str(QOI),'.gif'];
for i=1:N
    % Euclidean distance
    d = cumsum([0; sqrt(sum((X(I(2:end,i),:)-X(I(1:end-1,i),:)).^2,2))]);
    % Euclidean distance plot
    figure(fig1);
    plot(d,F(I(:,i)),'o-');
%     text(d,F(I(:,i)),num2str(I(:,i)));
    % parametric distance plot
    figure(fig2);
    plot(t,F(I(:,i)),'o-');
%     text(t,F(I(:,i)),num2str(I(:,i)));
    % projection over longest length
    w = 2*ones(m,1); w = w/norm(w);
    % plot projection
    figure(fig3);
    plot(X(I(:,i),:)*w,F(I(:,i)),'o-');
%     text(X0(I(:,i),:)*w,F(I(:,i)),num2str(I(:,i)));
    %% plot geometry over sweeps
    % build finer convex combinations    
    XX = kron(X(I(1,i),:)',(1-tt))' + kron(X(I(end,i),:)',tt)';
    XX0 = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(XX+1)));
    % build spline of sweep
    spl = pchip(X(I(:,i),:)*w,F(I(:,i)));
    % build figure
    figure(fig4);
    subplot(2,1,1), plot(X(I(:,i),:)*w,F(I(:,i)),'k.-','MarkerSize',8); hold on;
%     text(X(I(:,i),:)*w,F(I(:,i)),num2str(I(:,i)));
    for ii=1:TT
        % generate refined airfoil over sweep
        [coordU,coordL] = cst_airfoil(l',XX0(ii,1:m/2),XX0(ii,m/2+1:m),0);
        % sweep plot
        figure(fig4); subplot(2,1,1), h1 = scatter(XX(ii,:)*w,ppval(spl,XX(ii,:)*w),50,'filled','r');
        % airfoil plots
        figure(fig4); subplot(2,1,2), h2 = plot(l,coordU(2,:),'b','LineWidth',2); hold on; 
        figure(fig4); subplot(2,1,2), h3 = plot(l,coordL(2,:),'b','LineWidth',2); axis equal;        
        fig4.CurrentAxes.Visible = 'off';
        % build gif
        figure(fig4); frame = getframe(fig4);
        [A,map] = rgb2ind(frame2im(frame),256);
        if i*ii == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
        % delete plots
        delete([h1,h2,h3]);
    end
end
% string for legend
sweepstr = [repmat('sweep-',N,1),num2str((1:N)')];
% label plots
figure(fig1); xlabel('x(t)'); ylabel('f(x(t))'); legend(sweepstr);
figure(fig2); xlabel('t'); ylabel('f(t)'); legend(sweepstr);
figure(fig3); xlabel('x(t)^Tw'); ylabel('f(x(t)^Tw)'); legend(sweepstr);