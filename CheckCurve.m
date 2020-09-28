close all;clear;clc;
% warning off;

%% Calculation
global selectPlane NselectPlane err;
selectPlane = [1,2]; % select plane 1-z 2-rho 3-pz 4-prho
NselectPlane = setdiff(1:4,selectPlane);
err=1e-4; % calculation error
z = 0;r = 1.652;pz = 0.0225;pr = 0; % periodic
% z = 0;r = 1.616;pz = 0.0352;pr = 0; % periodic
% z = 0;r = 1.63;pz = 0.01;pr = 0;
te = 1e4;
%% Calculation
y0=[z,r,pz,pr]';
ts = 0:0.1:te;
odeFunc = @(t,x)chao_rhs_ext(t,x);
opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-10);
[ts,y]=ode45(odeFunc,ts,y0,opts);
fig = figure;
set(fig, 'position', get(0,'ScreenSize')); % Fullscreen
plot(y(:,selectPlane(1)),y(:,selectPlane(2)));
hold on;
grid on;
title(['\{',num2str(z),',',num2str(r),',',num2str(pz),',',num2str(pr),'\}'],'FontSize',24);
set(gca,'FontSize',20);
names = {'z','\rho','p_z','p_\rho'};
xlabel(names{selectPlane(1)});
ylabel(names{selectPlane(2)});

%% subfunction
function f=chao_rhs_ext(t,x)
z=x(1); r=x(2); pz=x(3); pr=x(4);
f=[ pz, pr, (3*r*z*(r/(r^2 + z^2)^(3/2) - 1/r))/(r^2 + z^2)^(5/2), -(r/(r^2 + z^2)^(3/2) - 1/r)*(1/(r^2 + z^2)^(3/2) - (3*r^2)/(r^2 + z^2)^(5/2) + 1/r^2)]';
end