close all; clearvars;
plt_flg = 0;

% Mesh
np = 250;
l = [linspace(0,0.005,50),linspace(0.005+.005/50,1,np)];

%% Generate Example Nominal
% Input is a vector of PARSEC parameters p=[p1, p2, ...pn] where dv=p s.t.
% p1=rle         
% p2=Xup
% p3=Yup
% p4=YXXup
% p5=Xlow
% p6=Ylow
% p7=YXXlow
% p8=yte
% p9=delta yte (t.e. thickness)
% p10=alpha te
% p11=beta te
% NACA 0012  from Sean Wu (WRONG)
dv  = [0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,0,0,-2.7791,9.2496];
m =length(dv) - 2;
% Sharp Trailing Edge
% dv  = [0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,0,0,-2.7791,9.2496];

a  = parsec(dv);
bf = zeros(length(l),6);
for i = 1:6
    bf(:,i) = (l').^(i-1/2);
end

%% Plot Nominal Airfoil
if plt_flg == 1
    plot(l',bf*a(1:6),'LineWidth',4); hold on
    plot(l',-bf*a(7:end),'LineWidth',4)
    axis equal;
end

%% Random Perturbation from Hypercube
% Sample uniformly from box constraints in 2*N+2 dimensions
N = 1000;
% rng(43); % Used for 500 Sample Designs
rng(44); % Used for 1000 Sample Designs
% rng(45); % Used for 5000 Sample Designs

% Define Reasonable upper and lower bounds
pct = 0.2;
lb0    = (1-pct)*[dv(1:7),dv(10:end)];
% lb0(2) = 0.2;
ub0    = (1+pct)*[dv(1:7),dv(10:end)];
% ub0(2) = 0.6;

%disp('----- Parameter Bounds -----')
%disp(['Param.', '      ','lb','      ', 'ub']);
%disp([['rle  ';'X1   ';'Y1   ';'Yxx1 ';'X2   ';'Y2   ';'Yxx2 ';'Yte  ';'thk  ';'alpha';'beta '], ...
%    [' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' '],num2str(lb0'), ...
%    [' ';' ';' ';' ';' ';' ';' ';' ';' ';' ';' '],num2str(ub0')]);

X = 2*rand(N,length(dv)-2)-1;
% X = [XC1';XC2';XC3';XC4'];
X0 = bsxfun(@plus,lb0,bsxfun(@times,ub0-lb0,0.5*(X+1)));

% Sample from constraint polytope
NF = 0;
IP = zeros(N,1);
for i=1:N
    % Evaluate PARSEC
    a = parsec([X0(i,1:7),0,0,X0(i,8:end)]');
    fhU = bf*a(1:6); coordU = [l; fhU']; 
    fhL = -bf*a(7:end); coordL = [l; fhL'];
    if min(fhU-fhL) >= -1e-3
        IP(i) = 1;
        %  verify perturbations with plots
        if plt_flg == 1
        h1 = plot(l,coordU(2,:),'k--','LineWidth',2); hold on; h2 = plot(l,coordL(2,:),'k--','LineWidth',2);
        title(['Example Nominal Airfoil Shape',' & Pert i = ',num2str(i)])
        pause(0.25); delete([h1,h2]);
        end
    else
        NF = NF + 1;
    end
end
IP = logical(IP);
disp(['Number of feasible designs: ',num2str(sum(IP))])

% Reassign sample to feasible polytope
X  =  X(IP,:);
X0 =  X0(IP,:);

%% Generate Meshes
for i=1:N-NF
    
    % Evaluate PARSEC
    a = parsec([X0(i,1:7),0,0,X0(i,8:end)]');
    fhU = bf*a(1:6);
    fhL = bf*a(7:end);
    coordU = [l;fhU'];
    coordL = [l;-fhL'];
    
    mesh_coords(coordU,coordL,i,N,NF);
    
end

%% Save design workspace
save(['./designs/','PARSEC_designs_m',num2str(m),'_N',num2str(N),'_pm',num2str(pct*100),'pct.mat'],'X','X0','NF','IP','lb0','ub0','l','pct','N','m');
