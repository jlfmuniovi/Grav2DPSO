function [F] = matrixF(nstation,grid)
% MATRIXF function calculates the continuous forward operator at each 
% station point 
% Authors: Zulima Fernández-Muñiz and Juan Luis Fernández-Martínez
ncels = grid.nz*grid.nx;
zm = grid.zm;
xm = grid.xm;
xk = grid.xk;
zk = grid.zk;
F = zeros(nstation,ncels);

for k=1:nstation
    f = zeros(grid.nz,grid.nx);
    for j=1:grid.nz
        for i=1:grid.nx
            A = xk(k)-xm(i);
            B = xk(k)-xm(i+1);
            C = zk(k)-zm(j);
            D = zk(k)-zm(j+1);
            T1 = A*(log(A^2+D^2)-log(A^2+C^2));
            T2 = B*(log(B^2+D^2)-log(B^2+C^2));
            T3 = 2*D*(atan(A/D)-atan(B/D));
            T4 = 2*C*(atan(A/C)-atan(B/C));
            f(j,i) = T1-T2+T3-T4;
        end
    end
    f = f'; 
    F(k,:) = f(:);
end
