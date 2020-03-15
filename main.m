close all; clear all; clc % Clean all variables.

Undata=load('/Users/yuhongliu/Downloads/Testdata.mat');
Undata=Undata.Undata; % Obatin the data from the struct.

L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); 
% We only need the first n, since it's periodic.
y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k); 
% Rescaling by 2pi/L since fft assumes 2pi period.
[X,Y,Z]=meshgrid(x,y,z); % Create spatial domain.
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); % Create spectral domain.

%%

% In spatial domain.
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n); % The data was reshaped to a row vector.
close all, isosurface(X,Y,Z,abs(Un),0.4) % Visualize the data.
axis([-20 20 -20 20 -20 20]), grid on, drawnow, pause(1)
xlabel('x'), ylabel('y'), zlabel('z')
title('Visualisation for spatial domain data')
end

%%

% In frequency domain.
for j = 1:20
    Un(:,:,:) = reshape(Undata(j,:),n,n,n);
    Unt = fftn(Un); % The data is multivariable.
    Untshift = fftshift(Unt)/max(max(max(abs(Unt)))); 
    % the result of FFTn is shifted, to view the mathematically correct
    % data in the following line, we need to fftshift it first and also
    % rescale it by its maximum value.
    close all, isosurface(Kx,Ky,Kz,abs(Untshift),0.4)
    % The value input should be absolute value since there are imaginary parts.
    axis([-8 8 -8 8 -8 8]), grid on, drawnow, pause(1)
    xlabel('kx'), ylabel('ky'), zlabel('kz')
    title('Visualisation for frequency domain data')
end


%%

% Averaging the signals.
Unave=zeros(n,n,n);
Untave=zeros(n,n,n);
for i = 1:20
    Un(:,:,:) = reshape(Undata(i,:),n,n,n);
    Unave=Unave+Un;
    Unt = fftn(Un);
    Untave=Untave+Unt;
end
Unave=Unave/10.0; % The value was too small to be seen, so 2.0 is been multiplied.
Untave=fftshift(Untave)/20.; % Like before, FFTn shifted the result.

figure(1)
% The averaged spatial data.
isosurface(X,Y,Z,abs(Unave),0.4) % To avoid the imaginary part.
axis([-20 20 -20 20 -20 20]), grid on, drawnow
xlabel('x'), ylabel('y'), zlabel('z')
title('Visualisation for averaged spatial domain data')

figure(2)
% The averaged frequency data.
isosurface(Kx,Ky,Kz,abs(Untave)/max(max(max(abs(Untave)))),0.4)
axis([-8 8 -8 8 -8 8]), grid on, drawnow
xlabel('kx'), ylabel('ky'), zlabel('kz')
title('Visualisation for averaged frequency domain data')

figure(3)
% The averaged frequency data transformed back to spatial domain.
isosurface(X,Y,Z,abs(ifftn(fftshift(Untave)))*2.0,0.4) 
% The value was too small to show on the graph.
% We need to re-fftshift it back.
axis([-20 20 -20 20 -20 20]), grid on, drawnow
xlabel('x'), ylabel('y'), zlabel('z')
title('Visualisation for averaged frequency domain data in spatial domain')

abs_Untave = abs(Untave)/max(max(max(abs(Untave))));
[omega_x_ind, omega_y_ind, omega_z_ind] = ind2sub(size(abs_Untave), find(abs_Untave == max(abs_Untave(:))));
% The results are all index.
omega_x = ks(omega_y_ind); 
% We look up the coordinate in ks, the mathematically correct.
omega_y = ks(omega_x_ind);
omega_z = ks(omega_z_ind);

% The data still seems pretty noisy.
% But according to the frequency domain, we should filter around omega_x,
% omega_y, omega_z = 1.89, -1.05, 0.0.

%%

% We try to filter it using Gaussian filter.

filter=exp(-.2*((fftshift(Kx)-omega_x).^2+(fftshift(Ky)-omega_y).^2+(fftshift(Kz)-omega_z).^2));
% 3D filter.

marble_pos=zeros(20, 3);
for i=1:20
    Un(:,:,:) = reshape(Undata(i,:),n,n,n);
    Unt = fftn(Un);
    Unft = filter.*Unt; % Piecewise multiplication.
    Unf=ifftn(Unft); % No need to shift since ifftn takes that into consideration.
    
    % close all, isosurface(X,Y,Z,abs(Unf),0.4) 
    % axis([-20 20 -20 20 -20 20]), grid on, drawnow, pause(1)
    isosurface(X,Y,Z,abs(Unf),0.4),
    xlabel('x'), ylabel('y'), zlabel('z')
    title('Marble path in spatial domain')
    axis([-20 20 -20 20 -20 20]), grid on, drawnow
    
    %close all, isosurface(Kx,Ky,Kz,abs(fftshift(Unft))/max(max(max(abs(fftshift(Unft))))),0.4) 
    % xlabel('kx'), ylabel('ky'), zlabel('kz')
    % title('Marble frequency')
    %axis([-8 8 -8 8 -8 8]), grid on, drawnow, pause(1)
    
    % Now we try to extract the path.
    % We use the peak value's location as the marble's location.
    abs_Unf = abs(Unf);
    [pos_x, pos_y, pos_z] = ind2sub(size(abs_Unf), find(abs_Unf == max(abs_Unf(:))));
    marble_pos(i,:)=[x(pos_y) x(pos_x) x(pos_z)]; % The x and y are flipped somehow.
end

% We found the marble!

%%

% Plot the path by itself.
figure(3)
plot3(marble_pos(:,1), marble_pos(:,2), marble_pos(:,3), 'r-', 'Linewidth',[5])
xlabel('x'), ylabel('y'), zlabel('z')
title('Marble path in spatial domain')
grid on

%%

% The final position of the marble at the 20th data measurement is:
% x = -5.6250, y = 4.2188, z=-6.0938.
disp(marble_pos(20,:));

