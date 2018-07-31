% Create gif over active subspace
clearvars; close all;
% select subspace file to load
load ./subspaces/FINE/AS_m10_N1000_pm60pct_QOI2.mat

%% discretize and compute coordinates
% only sweep over the first active variable for drag
W = [sub.W1, sub.W2];
W1 = W(:,1); W2 = W(:,2:end);
% discretize over min and max observed active variable
T = 50; Y = linspace(min(X*W1),max(X*W1),T)'; Y = [Y; Y(end:-1:1)];
% compute active coordinates (uncomment inactive variables of interest)
Z = zeros(2*T,9); Zstr = '_Z0';
% rng(47); Z = repmat(2*rand(1,9) - 1,2*T,1); Zstr = '_Zrnd1'; %WARNING: This will extrapolate beyond previous bounds
% rng(48); Z = repmat(2*rand(1,9) - 1,2*T,1); Zstr = '_Zrnd2'; %WARNING: This will extrapolate beyond previous bounds
XX = Y*W1' + Z*W2'; XX0 = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(XX+1)));
% shadow plot
fig = figure;
subplot(2,1,1), scatter(X*W1,F,30,'filled'); hold on;
% function approximation plot
subplot(2,1,1), plot(Y1,Fh,'k--','linewidth',2);

%% plot geometry over active subspace
filename = ['./subspaces/AS_sweep_m',num2str(m),'_N',num2str(N),'_pm',num2str(pct*100),'pct_','QOI',num2str(QOI),Zstr,'.gif'];
for ii=1:2*T
    % generate airfoil over sweep
    [coordU,coordL] = cst_airfoil(l',XX0(ii,1:m/2),XX0(ii,m/2+1:m),0);
    % sweep plot
    subplot(2,1,1), h1 = scatter(Y(ii),poly_predict(Y(ii),Coef,p),50,'filled','r');
    % airfoil plots
    subplot(2,1,2), h2 = plot(l,coordU(2,:),'b','LineWidth',2); hold on; 
    subplot(2,1,2), h3 = plot(l,coordL(2,:),'b','LineWidth',2);        
    fig.CurrentAxes.Visible = 'off'; fig.CurrentAxes.DataAspectRatio = [1 1 3];
    % build gif
    figure(fig); frame = getframe(fig); 
    [A,map] = rgb2ind(frame2im(frame),256);
    if ii == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    % delete plots
    delete([h1,h2,h3]);    
end
