% Create gif over active subspace
clearvars; close all;

% Assign your custom active_subspaces home directory
AS_HOME = '~/active_subspaces';
addpath([AS_HOME,'/matlab/']);
addpath([AS_HOME,'/matlab/ResponseSurfaces']);
% select subspace file to load
load ./subspaces/FINE/AS_m10_N1000_pm60pct_QOI2.mat
% pick parametrization
PARSEC = 0;

%% discretize and compute coordinates
% only sweep over the first active variable for drag
W = [sub.W1, sub.W2];
W1 = W(:,1); W2 = W(:,2:end);
% discretize over min and max observed active variable
T = 50; Y = linspace(min(X*W1),max(X*W1),T)'; Y = [Y; Y(end:-1:1)];
% compute active coordinates (uncomment inactive variables of interest)
Z = zeros(2*T,m-1); Zstr = '_Z0';
% rng(47); Z = repmat(2*rand(1,m-1) - 1,2*T,1); Zstr = '_Zrnd1'; %WARNING: This will extrapolate beyond bounds
% rng(48); Z = repmat(2*rand(1,m-1) - 1,2*T,1); Zstr = '_Zrnd2'; %WARNING: This will extrapolate beyond bounds
XX = Y*W1' + Z*W2'; XX0 = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(XX+1)));
% determine valid designs
indu = max(( XX(1:T,:) - ones(T,m) > 0 ),[],2); indl = max(( -ones(T,m) - XX(1:T,:) > 0),[],2);
% shadow plot
fig = figure;
subplot(2,1,1), scatter(X*W1,F,30,'filled'); hold on;
% function approximation plot
% visualize one-dimensional approximation
[Coef] = poly_train(X*W1, F, p);
lb = min(X*W1); ub = max(X*W1);
Y1 = linspace(lb,ub,NY)';
p = 1; [Fh,df] = poly_predict(Y1,Coef,p);
subplot(2,1,1), plot(Y1,Fh,'k--','linewidth',2);
% plot extrapolated points
subplot(2,1,1), plot([min(XX(indu,:)*W1),max(XX(indu,:)*W1)],...
                     poly_predict([min(XX(indu,:)*W1); max(XX(indu,:)*W1)],Coef,p),'rx-','LineWidth',2,'MarkerSize',10); 
subplot(2,1,1), plot([min(XX(indl,:)*W1),max(XX(indl,:)*W1)],...
                     poly_predict([min(XX(indl,:)*W1); max(XX(indl,:)*W1)],Coef,p),'rx-','LineWidth',2,'MarkerSize',10);
% compute nominal shape
M  = diag(2./(ub0-lb0)); b = -(M*lb0' + ones(m,1));
x0 = -inv(M)*b; 
if PARSEC == 1
        bf = zeros(length(l),6);
        for ii = 1:6
            bf(:,ii) = (l').^(ii-1/2);
        end
        a = parsec([x0(1:7)',0,0,x0(8:end)']');
        fhU = bf*a(1:6); coordU0 = [l; fhU']; 
        fhL = -bf*a(7:end); coordL0 = [l; fhL'];
else
    [coordU0,coordL0] = cst_airfoil(l',x0(1:m/2),x0(m/2+1:m),0);
end


%% plot geometry over active subspace
filename = ['./subspaces/AS_sweep_m',num2str(m),'_N',num2str(N),'_pm',num2str(pct*100),'pct_','QOI',num2str(QOI),Zstr,'.gif'];
% plot nominal shape
subplot(2,1,2), plot(l,coordU0(2,:),'--','LineWidth',2,'color',0.5*ones(3,1)); hold on; axis equal;
ax = subplot(2,1,2); plot(l,coordL0(2,:),'--','LineWidth',2,'color',0.5*ones(3,1));
for i=1:2*T
    % generate airfoil over sweep
    if PARSEC == 1
        bf = zeros(length(l),6);
        for ii = 1:6
            bf(:,ii) = (l').^(ii-1/2);
        end
        a = parsec([XX0(i,1:7),0,0,XX0(i,8:end)]');
        fhU = bf*a(1:6); coordU = [l; fhU']; 
        fhL = -bf*a(7:end); coordL = [l; fhL'];
    else
        [coordU,coordL] = cst_airfoil(l',XX0(i,1:m/2),XX0(i,m/2+1:m),0);
    end
    % sweep plot
    set(0,'defaulttextInterpreter','latex')
    subplot(2,1,1), h1 = scatter(Y(i),poly_predict(Y(i),Coef,p),50,'filled','r'); xlabel '$$w_1^Tx$$'; ylabel '$$C_{\ell}$$';
    % airfoil plots
    if  max(( XX(i,:) - ones(1,m) > 0 ),[],2) == 0 && max(( -ones(1,m) - XX(i,:) > 0),[],2) == 0
        subplot(2,1,2), h2 = plot(l,coordU(2,:),'b','LineWidth',2); hold on; 
        subplot(2,1,2), h3 = plot(l,coordL(2,:),'b','LineWidth',2); axis([ax.XLim, ax.YLim]);
        fig.CurrentAxes.Visible = 'off'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    else
        subplot(2,1,2), h2 = plot(l,coordU(2,:),'r','LineWidth',2); hold on; 
        subplot(2,1,2), h3 = plot(l,coordL(2,:),'r','LineWidth',2); axis([ax.XLim, ax.YLim]);
        fig.CurrentAxes.Visible = 'off'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    end
     % build gif
    figure(fig); frame = getframe(fig); 
    [A,map] = rgb2ind(frame2im(frame),256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    % delete plots
    delete([h1,h2,h3]);    
end
