function [gcal,mref] = solve_model(model,opfun)
% SOLVE_MODEL calculates the density in each cell of the grid and the
% gravity anomaly provides for each dense body
% Authors: Zulima Fernández-Muñiz and Juan Luis Fernández-Martínez
nanomalies = opfun.nanomalies;
nstation = opfun.nstation;
G = opfun.G;
nx = opfun.nx;
nz = opfun.nz;
dx = opfun.dx;
dz = opfun.dz;
F = opfun.F;
density = zeros(nz,nx);
rhob = model(1);
density(:,:) = rhob-rhob;
for k=1:nanomalies
    rhoi = model(5*(k-1)+2)-model(1);
    nx1 = ceil(model(5*(k-1)+3)/dx);
    nx2 = floor(model(5*(k-1)+4)/dx);
    nz1 = ceil(model(5*(k-1)+5)/dz);
    nz2 = floor(model(5*(k-1)+6)/dz);
    if nx1 > nx2
        temp = nx2;
        nx2 = nx1;
        nx1 = temp;
    end
    if nz1 > nz2
        temp = nz2;
        nz2 = nz1;
        nz1 = temp;
    end
    density(nz1:nz2,nx1:nx2) = rhoi;
end
dens = density';
gcal = G*F*dens(:);
mref = density+rhob;