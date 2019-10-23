function [misfit] = fcost(swarm,data,opciones,opfun)
% FCOST  function calculates the density value in each cell of the domain
% and the error misfit data 
% Authors: Zulima Fernández-Muñiz and Juan Luis Fernández-Martínez
rhob = swarm(:,1);
nanomalies = opfun.nanomalies;
nstation = opfun.nstation;
G = opfun.G;
nx = opfun.nx;
nz = opfun.nz;
dx = opfun.dx;
dz = opfun.dz;
F = opfun.F;
for i=1:size(swarm,1)
    gcal = zeros(1,nstation);
    density = zeros(nz,nx);
    rhob = swarm(i,1);
    density(:,:) = rhob-rhob;
    for k=1:nanomalies
        rhoi= swarm(i,5*(k-1)+2)-swarm(i,1);
        nx1 = ceil(swarm(i,5*(k-1)+3)/dx);
        nx2 = floor(swarm(i,5*(k-1)+4)/dx);
        nz1 = ceil(swarm(i,5*(k-1)+5)/dz);
        nz2 = floor(swarm(i,5*(k-1)+6)/dz);
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
    misfit(i) = norm(gcal(:)-data.gobs(:),opfun.norm)/norm(data.gobs(:),opfun.norm);
end