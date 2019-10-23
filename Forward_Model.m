% FORWARD_MODEL generates a true model formed by 2 dense bodies and 
% calculates the observed gravity perturbed with noise
% Authors: Zulima Fernández-Muñiz and Juan Luis Fernández-Martínez
clear all
G = 6.67*10^-06; % Gravity constant for obtaining g in mGal
%
% Stations (measured points)
nstation = 20; % number of measure points
rgrid.xki = 200; % coordinate x of the first point of measurement
rgrid.xkf = 800; % coordinate x of the last point of measurement
rgrid.zki = 0; % coordinate z of the first point of measurement
rgrid.zkf = 0; % coordinate z of the last point of measurement
rgrid.xk = linspace(rgrid.xki,rgrid.xkf,nstation); % coordinates x of the measuremet points
rgrid.zk = linspace(rgrid.zki,rgrid.zkf,nstation); % coordinates z of the measuremet points
%
% Grid
rgrid.nx = 1000; % number of cells of the grid in the x direction 
rgrid.nz = 200; % number of cells of the grid in the z direction 
rgrid.xmin = 0; % first coordinate x of the grid
rgrid.xmax = 1000; % last coordinate x of the grid
rgrid.zmin = 1; % first coordinate z of the grid
rgrid.zmax = 200; % last coordinate z of the grid
rgrid.xm = linspace(rgrid.xmin,rgrid.xmax,rgrid.nx+1); % x coordinates of the grid
rgrid.zm = linspace(rgrid.zmin,rgrid.zmax,rgrid.nz+1); % z coordinates of the grid
dx = (rgrid.xmax-rgrid.xmin)/rgrid.nx ; % length of each cell in the x direction
dz = (rgrid.zmax-rgrid.zmin)/rgrid.nz ; % length of each cell in the z direction
%
% Dense Bodies
nvertex = 4; % Number of vertex of the anomalies shape
nanomalies = 2; % Number of anomalies
% 1
rgrid.rhok1 = 3500; % body 1 density
rgrid.xa1 = [160 201 299 251]; % x coordinates of the anomaly 1
rgrid.za1 = [90 70 110 140]; % z coordinates of the anomaly 1
rgrid.xmin1 = min(rgrid.xa1);
rgrid.xmax1 = max(rgrid.xa1);
rgrid.zmin1 = min(rgrid.za1);
rgrid.zmax1 = max(rgrid.za1);
colum1 = ceil([rgrid.xmin1/dx,rgrid.xmax1/dx]); 
row1 = ceil([rgrid.zmin1/dz,rgrid.zmax1/dz]);
% 2
rgrid.rhok2 = 1000; % body 2 density
rgrid.xa2 = [570 750 800 600]; % x coordinates of the anomaly 2
rgrid.za2 = [130 100 170 150]; % z coordinates of the anomaly 2
rgrid.xmin2 = min(rgrid.xa2);
rgrid.xmax2 = max(rgrid.xa2);
rgrid.zmin2 = min(rgrid.za2);
rgrid.zmax2 = max(rgrid.za2);
colum2 = ceil([rgrid.xmin2/dx,rgrid.xmax2/dx]);
row2 = ceil([rgrid.zmin2/dz,rgrid.zmax2/dz]);
%
% Background
rgrid.rhob = 2720; % mean of background density
rhomodel = zeros(rgrid.nz,rgrid.nx);
rhomodel(1:rgrid.nz,1:rgrid.nx) = rgrid.rhob;
%
rhomodel(row1(1):row1(2),colum1(1):colum1(2)) = rgrid.rhok1;
rhomodel(row2(1):row2(2),colum2(1):colum2(2)) = rgrid.rhok2;
%
cxmodel1 = rgrid.xmin1:rgrid.xmax1; % x coordinates of the rectangular model 1
czmodel1 = rgrid.zmin1:rgrid.zmax1;
cxmodel2 = rgrid.xmin2:rgrid.xmax2; % x coordinates of the rectangular model 2
czmodel2 = rgrid.zmin2:rgrid.zmax2;
%
for i = 1:length(cxmodel1)
    x = cxmodel1(i);
    for j=1:length(czmodel1)
        z = czmodel1(j);
        vecx = rgrid.xa1-x;
        vecz = rgrid.za1-z;
        v1 = [vecx(1) vecz(1) 0];
        v2 = [vecx(2) vecz(2) 0];
        v3 = [vecx(3) vecz(3) 0];
        v4 = [vecx(4) vecz(4) 0];
        col = ceil(x/dx);
        row = ceil(z/dz);
        signo = [sign(cross(v1,v2)),sign(cross(v2,v3)),sign(cross(v3,v4)),sign(cross(v4,v1))];
        signo = [signo(3) signo(6) signo(9) signo(12)];
        if isequal(signo,ones(1,4))==1 | isequal(signo,-ones(1,4))==1
            rhomodel(row,col) = rgrid.rhok1;
        else
            rhomodel(row,col) = rgrid.rhob;
        end
    end
end
%
for i = 1:length(cxmodel2)
    x = cxmodel2(i);
    for j=1:length(czmodel2)
        z = czmodel2(j);
        vecx = rgrid.xa2-x;
        vecz = rgrid.za2-z;
        v1 = [vecx(1) vecz(1) 0];
        v2 = [vecx(2) vecz(2) 0];
        v3 = [vecx(3) vecz(3) 0];
        v4 = [vecx(4) vecz(4) 0];
        col = ceil(x/dx);
        row = ceil(z/dz);
        signo = [sign(cross(v1,v2)),sign(cross(v2,v3)),sign(cross(v3,v4)),sign(cross(v4,v1))];
        signo = [signo(3) signo(6) signo(9) signo(12)];
        if isequal(signo,ones(1,4))==1 | isequal(signo,-ones(1,4))==1
            rhomodel(row,col) = rgrid.rhok2;
        else
            rhomodel(row,col) = rgrid.rhob;
        end
    end
end
rho = rhomodel';
rgrid.rho = rho(:); % density values for all cells (dense bodies and background)
rgrid.rhomodel = rhomodel;
%
% Calculating noisy gravimetry anomaly
F = matrixF(nstation,rgrid);
gcal = G*F*(rgrid.rho-rgrid.rhob); % theoretical gravity anomaly
percentage = 0.05; % percentage of noise
gobs = gcal(:)+randn(length(gcal),1)*mean(gcal)*percentage; % noisy gravity anomaly
%
% Storing results
opfun.G = G;  
opfun.F = F;
opfun.nanomalies = nanomalies; 
opfun.nstation = nstation;
opfun.nvertex = nvertex;
opfun.nx = rgrid.nx;
opfun.nz = rgrid.nz;
data.gobs = gobs;
opfun.dx = dx;
opfun.dz = dz;
%
