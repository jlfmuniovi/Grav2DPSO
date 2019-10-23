function [] = graphics(rgrid,data,mod,IQRm,mr,gp)
% Authors: Zulima Fernández-Muñiz and Juan Luis Fernández-Martínez
x = rgrid.xmin:rgrid.xmax;
z = rgrid.zmin:rgrid.zmax;
xa1d = [rgrid.xa1 rgrid.xa1(1)];
za1d = [rgrid.za1 rgrid.za1(1)];
xa2d = [rgrid.xa2 rgrid.xa2(1)];
za2d = [rgrid.za2 rgrid.za2(1)];
%
figure
subplot(2,1,1)
imagesc(x,z,mod)
hold on
plot(xa1d,za1d,'r')
plot(xa2d,za2d,'r')
grid on
colorbar
title ('Median model')
subplot(2,1,2)
imagesc(x,z,IQRm)
hold on
plot(xa1d,za1d,'k')
plot(xa2d,za2d,'k')
grid on
colorbar
title ('IQR')
%
figure
subplot(2,1,1)
imagesc(x,z,rgrid.rhomodel)
hold on
plot(xa1d,za1d,'r')
plot(xa2d,za2d)
title('True model')
grid on
colorbar
subplot(2,1,2)
imagesc(x,z,mr)
hold on
plot(xa1d,za1d,'r')
plot(xa2d,za2d,'r')
title ('Best model')
colorbar
grid on
%
figure
plot(rgrid.xk,data.gobs,'k-.')
hold on
plot(rgrid.xk,gp)
xlabel('Stations');
ylabel('g^{obs} , g^{pred} (mGal)');
grid on
legend('g observed','g predicted')
