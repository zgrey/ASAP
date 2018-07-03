% Function for reading SU2 ''forces*.dat'' in specified file directories
% Read su2 forces files format
% dir_str is the structure for directory of interest
% Quantity of Interest (QOI) QOI = 1: CL/CD, QOI = 2: CL, QOI = 2: CD

function [F, N, I, Fmax, Fmin, Imax, Imin] = readsu2_forces(dir_str,QOI)

% determine number of runs (should be consistent with designs_*.mat params
N = size(dir_str,1)-2;
% precondition vectors
F = zeros(N,1); I = zeros(N,1);
for i=1:N
    fname = ['./forces/',dir_str(i+2).name];
    copyfile(fname,'./dummy.dat')
    
    if QOI == 1
        % Lift over Drag
        if i==1, disp('READING LIFT/DRAG DATA...'); end
        !grep 'Total CL/CD:'  dummy.dat > grepdummy.dat
    elseif QOI == 2
        % Lift
        if i==1, disp('READING LIFT DATA...'); end
        !grep 'Total CL:'  dummy.dat > grepdummy.dat
    elseif QOI == 3
        % Drag
        if i==1, disp('READING DRAG DATA...'); end
        !grep 'Total CD:'  dummy.dat > grepdummy.dat
    end
    
    delete dummy.dat
    
    % Import data from text file
    try 
        fid   = fopen('grepdummy.dat');
        fline = fgetl(fid);
        p1 = strfind(fline,':');
        p1 = p1(1);
        p2 = strfind(fline,'|');
        p2 = p2(1);
        p3 = strfind(dir_str(i+2).name,'.');
        p3 = p3(1);
        
        % Results files must have 'airfoil' prefix
        I(i) = str2double(dir_str(i+2).name(9:p3-1));
        F(I(i)) = str2double(fline(p1+1:p2-1));
    catch
        disp(['ERROR: File Index ', dir_str(i+2).name(9:p3-1)])
        I(i) = str2double(dir_str(i+2).name(9:p3-1));
        F(I(i)) = -10;
    end
    
    % close file
    fclose(fid);
end
% print summary statistics
% Check for min and max design
[Fmax,Imax] = max(F);
[Fmin,Imin] = min(F);
disp(['Max Design: # ',num2str(Imax),', Max QoI = ',num2str(Fmax)]);
disp(['Min Design: # ',num2str(Imin),', Min QoI = ',num2str(Fmin)]);
fclose all; delete grepdummy.dat;