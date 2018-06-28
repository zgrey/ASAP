% Mesh Generic Coordinates

function [] = mesh_coords(coordU,coordL,i,NN,NF)

%% Write Coords. to *.dat file
M = [coordU(:,end:-1:1)'; [0 0]; coordL'];
M = [M, zeros(size(M,1),1)];
dlmwrite('airfoil.dat',M,'delimiter','\t');

%% Mesh airfoil using gmsh 2.10.1
fprintf('calling dat2gmsh\n');
!python dat2gmsh.py airfoil.dat
fprintf('calling gmsh\n');
% may require sudo for encrypted hard drives
!gmsh airfoil.dat.geo -2 -o airfoil.mesh

%% Convert Mesh
tic;
fprintf('mesh %i of %i pre-processing/conversion to .su2 format...\n',i,NN-NF);
% Run DARPA EQUiPS SEQUOIA team codes:
meshGMF = ReadGMF('airfoil.mesh');
meshGMF = meshPrepro(meshGMF);
meshSU2 = convertGMFtoSU2(meshGMF);
WriteSU2_airfoil(meshSU2, ['./meshes/airfoil_',num2str(i),'.su2']);

tocp = toc;
% print time stats
disp(['Finished in... ',num2str(toc),' sec | Remaining time ~',num2str((NN-NF-i)*(toc+tocp)/120),' min']);