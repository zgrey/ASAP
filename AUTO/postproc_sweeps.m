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
QOI = 3; % QOI = 1: CL/CD, QOI = 2: CL, QOI = 3: CD 
% read force data
F = readsu2_forces(dir('./forces'),QOI);

%% Plot and visualize sweeps
% more points over sweep for geometry visualization
TT = 25; tt = linspace(0,1,TT);

% Index set (columns are sweep indices, rows are point indices in sweep)
I = reshape(linspace(1,N*T,N*T)',T,N);
fig1 = figure; hold on; fig2 = figure; hold on; fig3 = figure; hold on; fig4 = figure; hold on;
filename = 'sweeps.gif';
for i=1:N
    % Euclidean distance
    d = cumsum([0; sqrt(sum((X(I(2:end,i),:)-X(I(1:end-1,i),:)).^2,2))]);
    % Euclidean distance plot
    figure(fig1);
    plot(d,F(I(:,i)),'o-');
    text(d,F(I(:,i)),num2str(I(:,i)));
    % parametric distance plot
    figure(fig2);
    plot(t,F(I(:,i)),'o-');
    text(t,F(I(:,i)),num2str(I(:,i)));
    % projection over longest length
    w = 2*ones(m,1); w = w/norm(w);
    % plot projection
    figure(fig3);
    plot(X(I(:,i),:)*w,F(I(:,i)),'o-');
%     text(X0(I(:,i),:)*w,F(I(:,i)),num2str(I(:,i)));
    %% plot geometry over sweeps
%     % build finer convex combinations    
%     XX0 = kron(X0(I(1,i),:)',(1-tt))' + kron(X0(I(end,i),:)',tt)';
%     % build spline of sweep
%     spl = pchip(X0(I(:,i),:)*w,F(I(:,i)));
%     % build figure
%     figure(fig4);
%     subplot(2,1,1), plot(X0(I(:,i),:)*w,F(I(:,i)),'k.-','MarkerSize',8); hold on;
%     text(X0(I(:,i),:)*w,F(I(:,i)),num2str(I(:,i)));
%     for ii=1:TT
%         [s,~,M] = nasa_nozzle(linspace(0,10,100),[XX0(ii,1:3)'.^2*pi; XX0(ii,4)],0);
%         % compute aspect ratio of convex cells
%         r = sqrt(s(1:2:end)/pi); AR = repmat(100./(99*r'),1,5);
%         dx = max(diff(M(:,1))); dMy = diff(reshape(M(:,2),50,5)); dy = repmat(dMy(1:5:end,:),5,1);
%         % sweep plot
%         subplot(2,1,1), h = scatter(XX0(ii,:)*w,ppval(spl,XX0(ii,:)*w),50,'filled','r');
%         % nozzle plots
%         subplot(2,1,2), h6 = surf(M(1:5:end,1),M(:,2),dy/dx); view([0,0,1]); shading interp; colorbar; caxis([1,10]);
% %         subplot(2,1,2), h6 = surf(reshape(M(:,1),50,5),reshape(M(:,2),50,5),dy/dx); view([0,0,1]); shading interp; colorbar;
% %         subplot(2,1,2), h6 = scatter(reshape(M(:,1),50,5),reshape(M(:,2),50,5),100,'cdata',dy/dx);
%         subplot(2,1,2), hh = plot(linspace(0,10,100),sqrt(s/pi),'k-','linewidth',2); hold on; axis equal; axis([0,10,-xu(1),xu(1)]);
%         subplot(2,1,2), hhh = plot(linspace(0,10,100),-sqrt(s/pi),'k-','linewidth',2);
% %         subplot(2,1,2), hhhh = plot(reshape(M(:,1),50,5),reshape(M(:,2),50,5),'k.');
%         subplot(2,1,2), hhhhh = plot(reshape(M(:,1),50,5),-reshape(M(:,2),50,5),'k.'); fig4.CurrentAxes.Visible = 'off';        
%         % hollywood shit
%         frame = getframe(fig4);
%         [A,map] = rgb2ind(frame2im(frame),256);
%         if i*ii == 1
%             imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
%         else
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
%         end
%         % delete plots
%         delete(h); delete(hh); delete(hhh); delete(hhhh); delete(hhhhh);
%         delete(h6);
%     end
end
% string for legend
sweepstr = [repmat('sweep-',N,1),num2str((1:N)')];
% label plots
figure(fig1); xlabel('x(t)'); ylabel('f(x(t))'); legend(sweepstr);
figure(fig2); xlabel('t'); ylabel('f(t)'); legend(sweepstr);
figure(fig3); xlabel('x(t)^Tw'); ylabel('f(x(t)^Tw)'); legend(sweepstr);