clear all; 
close all; 
clc

% Load Data Import
v1 = VideoReader('monte_carlo_low.mp4');
v2 = VideoReader('v2_drop_low.mp4')

dt = 1/v1.Framerate;
t = 0:dt:v1.Duration;
vid1 = read(v1);
num = get(v1, 'numberOfFrames');

dt = 1/v2.Framerate;
t = 0:dt:v2.Duration;
vid2 = read(v2);
num = get(v2, 'numberOfFrames');

for i = 1:num
    mov(i).cdata = vid2(:,:,:,i);
    mov(i).colormap = [];
end

images = [];

images1 = images(:,1:end-1);
images2 = images(:,2:end);

[U,Sigma,V] = svd(images1,'econ');

energy = 0;
modes = 0;
while energy < 0.9
    modes = modes + 1
    energy = energy + (sig(modes).^2/sum(sig.^2))  
end

r = 2;
Ur = U(:,1:r);
Sigmar = Sigma(1:r,1:r);
Vr = V(:,1:r);

S = Ur'*images2*Vr*diag(1./diag(Sigmar));


[eV,D] = eig(S);
mu = diag(D); 
omega = log(mu)/dt;
Phi = Ur*eV;


y0 = Phi\images1(:,1); 

imagesmodes = zeros(length(y0),length(t)-1);
for iter = 1:(length(t)-1)
    imagesmodes(:,iter) = y0.*eimagesp(omega*t(iter));
end
imadmd = Phi*imagesmodes;


sparse = images1-abs(imadmd);

R = sparse.*(sparse<0);

images_low = R + abs(imadmd);
imagesp = sparse-R;
images_r = images_low + imagesp;

% image_con = reshape(images_r, [135,240,378]);
% image_back = reshape(images_dmd, [135,240,378]);
% image_fore = reshape(images_sparse, [135,240,378]);
% image_ori = reshape(images1, [135,240,378]);

image_con=reshape(images_r, [135,240,453]);
image=uint8(real(image_con));
imshow(image);drawnow


image_back=reshape(imadmd, [135,240,453]);
image=uint8(real(image_back));
imshow(image);drawnow


image_fore=reshape(imagesp, [135,240,453]);
image=uint8(real(image_fore));
imshow(image);drawnow


image_ori=reshape(images1, [135,240,453]);
image=uint8(real(image_ori));
imshow(image);drawnow

