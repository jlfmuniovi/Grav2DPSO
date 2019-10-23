function [parent]=initialpop(talla,lowlimit,upperlimit);
%..........................................................................
%AIM: generation of a random initial population  in [lowlimit,upperlimit]
%              
%VARIABLES DE ENTRADA:
%[nparam]  :    Number of model parameters.
% talla    :    Number of models
%[lowlimit]:    Low limit     
%[upperlimit]:  Upper limit
%[bits]:        Bits array
%VARIABLES DE SALIDA:
%[parent]:      population
%
%SINTAXIS:
%..........................................................................
%               [parent]=initialpob(nparam,talla,lowlimit,upperlimit);
%..........................................................................
% Random initial population in [lowlimit, upperlimit]
nparam=length(lowlimit);
rango=upperlimit-lowlimit;
parent=(ones(talla,1)*rango).*(rand(talla,nparam))+(ones(talla,1)*lowlimit);