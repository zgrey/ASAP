% Mix all samples
% Read SHDP Results
clc;
clear all;
close all;

%% Load AS Matlab Environment
addpath '~/active_subspaces/matlab/'
addpath '~/active_subspaces/matlab/Subspaces'
addpath '~/active_subspaces/matlab/ResponseSurfaces'
addpath '~/active_subspaces/matlab/Plotters'
addpath '~/active_subspaces/matlab/Domains'
warning('off','MATLAB:dispatcher:pathWarning')

%% Read ALL SU2 data files
% Pick QoI
QOI = 3; % QOI = 1: CL/CD, QOI = 2: CL, QOI = 3: CD 
sstype = 6;
Nbs = 500;

%% CST Parameterization
disp('CST Parameterization Post Processing...')
F1 = []; F2 = []; F3 = [];
X1 = []; X2 = []; X3 = [];
for i=1:3
    if i == 1
        !rm -rf forces
        !unzip -qq ./500_Samples/CST_forces.zip 
        dir_str = dir('./forces');
        [F1] = readsu2_forces(dir_str,QOI);
        !cp ./500_Samples/CST_designs*.mat ./designs.mat
        load('designs.mat');
        X1 = X;
    elseif i == 2
        !rm -rf forces
        !unzip -qq ./1000_Samples/CST_forces.zip
        dir_str = dir('./forces');
        [F2] = readsu2_forces(dir_str,QOI);
        !cp ./1000_Samples/CST_designs*.mat ./designs.mat
        load('designs.mat');
        X2 = X;
    elseif i == 3
        !rm -rf forces
        !unzip -qq ./5000_Samples/CST_forces.zip
        dir_str = dir('./forces');
        [F3] = readsu2_forces(dir_str,QOI);
        !cp ./5000_Samples/CST_designs*.mat ./designs.mat
        load('designs.mat');
        X3 = X;
    end
end

CST.F = [F1;F2;F3];
CST.X = [X1; X2; X3];
disp('Approximating CST active subspace...')
CST.sub = compute(CST.X,CST.F,[],[],sstype,0,Nbs); 

%% CNACA Parameterization
disp('CNACA Parameterization Post Processing...')
F1 = []; F2 = []; F3 = [];
X1 = []; X2 = []; X3 = [];
for i=1:3
    if i == 1
        !rm -rf forces
        !unzip -qq ./500_Samples/CNACA_forces.zip
        dir_str = dir('./forces');
        [F1] = readsu2_forces(dir_str,QOI);
        !cp ./500_Samples/CNACA_designs*.mat ./designs.mat
        load('designs.mat');
        X1 = X;
    elseif i == 2
        !rm -rf forces
        !unzip -qq ./1000_Samples/CNACA_forces.zip
        dir_str = dir('./forces');
        [F2] = readsu2_forces(dir_str,QOI);
        !cp ./1000_Samples/CNACA_designs*.mat ./designs.mat
        load('designs.mat');
        X2 = X;
    elseif i == 3
        !rm -rf forces
        !unzip -qq ./5000_Samples/CNACA_forces.zip
        dir_str = dir('./forces');
        [F3] = readsu2_forces(dir_str,QOI);
        !cp ./5000_Samples/CNACA_designs*.mat ./designs.mat
        load('designs.mat');
        X3 = X;
    end
end

CNACA.F = [F1;F2;F3];
CNACA.X = [X1; X2; X3];
disp('Approximating CNACA active subspace...')
CNACA.sub = compute(CNACA.X,CNACA.F,[],[],sstype,0,Nbs); 

%% PARSEC Parameterization
disp('PARSEC Parameterization Post Processing...')
F1 = []; F2 = []; F3 = [];
X1 = []; X2 = []; X3 = [];
for i=1:3
    if i == 1
        !rm -rf forces
        !unzip -qq ./500_Samples/PARSEC_forces.zip
        dir_str = dir('./forces');
        [F1] = readsu2_forces(dir_str,QOI);
        !cp ./500_Samples/PARSEC_designs*.mat ./designs.mat
        load('designs.mat');
        X1 = X;
    elseif i == 2
        !rm -rf forces
        !unzip -qq ./1000_Samples/PARSEC_forces.zip
        dir_str = dir('./forces');
        [F2] = readsu2_forces(dir_str,QOI);
        !cp ./1000_Samples/PARSEC_designs*.mat ./designs.mat
        load('designs.mat');
        X2 = X;
    elseif i == 3
        !rm -rf forces
        !unzip -qq ./5000_Samples/PARSEC_forces.zip
        dir_str = dir('./forces');
        [F3] = readsu2_forces(dir_str,QOI);
        !cp ./5000_Samples/PARSEC_designs*.mat ./designs.mat
        load('designs.mat');
        X3 = X;
    end
end

PARS.F = [F1;F2;F3];
PARS.X = [X1; X2; X3];
% Change Indices to Match Paper
Xorig = PARS.X;
PARS.X(:,1) = Xorig(:,2);
PARS.X(:,2) = Xorig(:,5);
PARS.X(:,4) = Xorig(:,6);
PARS.X(:,5) = Xorig(:,8);
PARS.X(:,6) = Xorig(:,9);
PARS.X(:,7) = Xorig(:,10);
PARS.X(:,8) = Xorig(:,11);
PARS.X(:,9) = Xorig(:,4);
PARS.X(:,10) = Xorig(:,7);
PARS.X(:,11) = Xorig(:,1);
disp('Approximating PARSEC active subspace...')
PARS.sub = compute(PARS.X,PARS.F,[],[],sstype,0,Nbs); 

%% Plot Active Subspace
if exist('AS','var')
    clearvars AS;
    close all;
end

% Pick a parameterization
AS = CST;

if size(AS.sub.W1,2) > 2
    disp('WARNING: Gap not detected, reassigning eigenvectors')
    AS.sub.W1 = AS.sub.W1(:,1:2);
end

% Options
if exist('opts','var')
    clearvars opts;
end
opts.fontsize = 16;
if QOI == 1
    opts.title = 'C_l / C_d Subspace';
elseif QOI == 2
    opts.title = [];
elseif QOI == 3
    opts.title = [];
end

% Eigenvectors
opts.markersize = 10;
eigenvectors(AS.sub.W1,opts);
print('-f1','./figs/eig_vec','-djpeg')

% Subspace Errors
subspace_errors(AS.sub.sub_br,opts);
print('-f2','./figs/sub_err','-djpeg')

% Plot eigenvalue decay
eigenvalues(AS.sub.eigenvalues,AS.sub.e_br,opts);
grid on
print('-f3','./figs/eig_val','-djpeg')

% Summary plot
opts.ylabel = 'QoI';
opts.markersize = 30;
sufficient_summary(AS.X*AS.sub.W1,AS.F,opts);
grid off
print('-f4','./figs/summary','-djpeg')
print('-f5','./figs/summary2','-djpeg')
figure(5);
% map = colormap;
% [cmin,cmax] = caxis;
    
%% Approximate h(W1'*x)
if size(AS.sub.W1,2) == 2
    % Approximate h() for 2D Case
    figure()
    p = 2;

    [Coef, B, ~, ~, f_hat, r] = poly_train(AS.X*AS.sub.W1, AS.F, p);
    NY = 50;
    lb = [min(min(AS.X*AS.sub.W1(:,1))), min(min(AS.X*AS.sub.W1(:,2)))];
    ub = [max(max(AS.X*AS.sub.W1(:,1))), max(max(AS.X*AS.sub.W1(:,2)))];
    [Y1,Y2] = meshgrid(linspace(lb(1),ub(1),NY),linspace(lb(2),ub(2),NY));
    [Fh,df] = poly_predict([reshape(Y1,NY*NY,1), reshape(Y2,NY*NY,1)],Coef,p);
    
    % Make Contour 2D sufficient summary plot.
    scatter(AS.X*AS.sub.W1(:,1), AS.X*AS.sub.W1(:,2),100, 'filled', 'cdata', AS.F, 'marker', 'o')
    hold on
    contour(Y1,Y2,reshape(Fh,NY,NY),35,'LineWidth',2);
%     colormap(map)
%     caxis([0,cmax])
    colorbar;
%     title([opts.title,' \approx h(W_{1}^{T}x)'],'fontsize',16)
    ylabel('Active Variable 2','fontsize',16);
    xlabel('Active Variable 1','fontsize',16);
    set(gca,'fontsize',16);
    print('-f6','./figs/summary_fit','-djpeg')
    
end

%% Subspace Error Convergence
k = 2.^[0:6]*100;
% PARSEC
disp('Computing PARSEC Subspace Errors...')
[PARS.serr,figerr,h3] = sub_conv(PARS,k,sstype,Nbs,[]);
hold on;

% % CNACA
% disp('Computing CNACA Subspace Errors...')
% [CNACA.serr,~,h2] = sub_conv(CNACA,k,sstype,Nbs,figerr);

% CST
disp('Computing CST Subspace Errors...')
[CST.serr,~, h1] = sub_conv(CST,k,sstype,Nbs,figerr);

% Reference Convergence
% loglog(k,1./sqrt(k),'--')

% Plot Format
figure(figerr.Number)
% h = legend([h1,h2,h3],{'CST','CNACA','PARSEC'});
h = legend([h1,h3],{'CST','PARSEC'});
if QOI == 1
    title(['C_l/C_d Subspace Convergence'])
elseif QOI == 2
    title(['C_l Subspace Convergence'])
elseif QOI == 3
    title(['C_d Subspace Convergence'])
end
print(['-f',num2str(figerr.Number)],'./figs/sub_conv','-djpeg')

%% Z - samples
clc;
Y = CST.X*CST.sub.W1;
% Pareto perturbations
NP = 100; % Try 50 Y_p conditionals

y_opt1 = [0, min(Y(:,2))];
y_opt2 = [min(Y(:,1)),0];
t = linspace(0,1,NP);
Y_p = bsxfun(@times,t',y_opt1)+bsxfun(@times,(1-t)',y_opt2);

% X | Z = 0
% Xpert = Y_p*CST.sub.W1'; % Must be drag subspace!!!!
NZ = 10;
Xpert = zeros(NP*NZ,size(CST.X,2));
k = 1;
j = zeros(NP,1);
for i = 1:NP
    Z = hit_and_run_z(NZ,Y_p(i,:)',CST.sub.W1,CST.sub.W2);
    if sum(sum(Z)) ~= 0
        j(i) = 1;
        for ii = 1:NZ
            Xpert(k,:) = Y_p(i,:)*CST.sub.W1' + Z(ii,:)*CST.sub.W2';
            k = k + 1;
        end
    end
end
j = logical(j);
Y_p = Y_p(j,:);
Xpert = Xpert(1:(k-1),:);