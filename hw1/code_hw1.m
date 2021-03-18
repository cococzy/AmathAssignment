% Clean workspace
clear; close all; clc

load subdata.mat

L = 10; % spatial domain
n = 64; % Fourier modes

x2 = linspace(-L,L,n+1); 
x = x2(1:n); y = x; 
z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; 
ks = fftshift(k);

[x1,y1,z1] = meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks,ks,ks);

Uave = zeros(n,n,n);
for j = 1:49
    Un(:,:,:) = fftn(reshape(subdata(:,j),n,n,n));
    Uave = Uave + Un;   
end

Uave = fftshift(Uave)./49;
[m, Index] = max(Uave(:));
[Ix, Iy, Iz] = ind2sub(size(Uave), Index);
x0 = Kx(Ix, Iy, Iz);  %  5.34
y0 = Ky(Ix, Iy, Iz);  % -6.91
z0 = Kz(Ix, Iy, Iz);  %  2.19

%%
x1 = []; y1 = []; z1 = [];
tau = -0.2;
filter=exp(tau*(((Kx-x0).^2)+((Ky-y0).^2)+((Kz-z0).^2)));

for j = 1:49
    Un = fftn(reshape(subdata(:,j),n,n,n));
    Unt = fftshift(Un);
    Unft = filter.*Unt;
    Unf = ifftn(Unft);
    [Max,idx] = max(Unf(:));
    [b,a,c] = ind2sub(size(Unf),idx);
    x1(j) = a;
    y1(j) = b;
    C(j) = c;
end

plot3(x(x1),y(y1),z(C))
title('Submarine Movement')
xlabel('x')
ylabel('y')
zlabel('z')
result= [x(x1); y(y1);z(C)]';