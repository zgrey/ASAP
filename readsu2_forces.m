% Function for reading SU2 ''forces*.dat'' in specified file directories

function [F, N, I, Fmax, Fmin, Imax, Imin] = readsu2_forces(dir_str,QOI)

% dir structure for directory of interest
% Quantity of Interest (QOI) QOI = 1: CL/CD, QOI = 2: CL, QOI = 2: CD

N = size(dir_str,1)-2;
disp('READING FORCES DATA...');
F = zeros(N,1);
I = zeros(N,1);
Fmax = -1e6;
Fmin = 1e6;
Imax = 0;
Imin = 0;
for i=1:N
    fname = ['./forces/',dir_str(i+2).name];
    copyfile(fname,'./dummy.dat')
    
    if QOI == 1
        % Lift over Drag
        !grep 'Total CL/CD:'  dummy.dat > grepdummy.dat
    elseif QOI == 2
        % Lift
        !grep 'Total CL:'  dummy.dat > grepdummy.dat
    elseif QOI == 3
        % Drag
        !grep 'Total CD:'  dummy.dat > grepdummy.dat
    end
    
    delete dummy.dat
    
    % Import data from text file
    fid   = fopen('grepdummy.dat');
    fline = fgetl(fid);
    p1 = strfind(fline,':');
    p1 = p1(1);
    p2 = strfind(fline,'|');
    p2 = p2(1);
    p3 = strfind(dir_str(i+2).name,'.');
    p3 = p3(1);
    
    % Results files must have 'airfoil' prefix
    I(i) = str2double(dir_str(i+2).name(8:p3-1));
    F(I(i)) = str2double(fline(p1+1:p2-1));
    
    fclose(fid);
    
    % Check for min and max design
    if max(F) > Fmax
        Fmax = max(F);
        Imax = I(i);
    end
    if min(F) < Fmin
        Fmin = min(F);
        Imin = I(i);
    end
end
disp(['Max Design: # ',num2str(Imax)]);
disp(['Min Design: # ',num2str(Imin)]);
fclose all;
delete grepdummy.dat;