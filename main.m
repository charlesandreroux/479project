

function [dat,j,time] = main(N, Tl, Tr, lambda, k, tstep, deltaM, maxStep, beta)

t = tic();

if ~exist('N', 'var')
    N = 1024;
end
if ~exist('Tl', 'var')
    Tl = 1.2;
end
if ~exist('Tr', 'var')
    Tr = 0.8;
end
if ~exist('lambda', 'var')
    lambda = 0.8;
end
if ~exist('k', 'var')
    k = 1;
end
if ~exist('tstep', 'var')
    tstep = 0.005;
end
if ~exist('deltaM', 'var')
    deltaM = 0;
end
if ~exist('maxStep', 'var')
    maxStep = 1e7;
end
if ~exist('beta', 'var')
    beta = 1;
end

kB = 1.380649e-23;

% arrays of positions, velocities, and accelerations
x = zeros(N + 2, 1);
v = zeros(N + 2, 1);
a = zeros(N + 2, 1);
v(2) = Tl;
v(N+1) = Tr;
xil0 = 2.*lambda.*kB.*Tl;
xir0 = 2.*lambda.*kB.*Tr;

%array of masses
m = 1 - deltaM + (2*deltaM)*rand(N + 2, 1);
m(1) = inf;
m(N+2) = inf;

newa = zeros(N+2,1);
xl = zeros(N+2,1);
xr = zeros(N+2,1);
fl = zeros(N+2,1);
fstep = zeros(N+2,1); 
cubel = zeros(N+2,1);
cuber = zeros(N+2,1);

fcum = zeros(N+2, 1);

for step = 1:1000

    x = x + v.*tstep + 0.5.*a.*tstep.^2;

    xl = circshift(x, 1);
    xl(N+2) = 0;

    cubel = (xl - x);
    cubel = cubel.*cubel.*cubel;
    cuber = -circshift(cubel,-1);
    cuber(N+1) = (x(N+1) - x(N))^3;
    cuber(1) = 0;

    xr = circshift(x,-1);
    xr(1) = 0;

%     x = gpuArray(x);
%     xl = gpuArray(xl);
%     xr = gpuArray(xr);
%     fl = gpuArray(fl);
%     cubel = gpuArray(cubel);
%     cuber = gpuArray(cuber);
%     v = gpuArray(v);
 
    fl = k.*xl + beta.*cubel;
    fstep = fl + k.*(- 2.*x + xr) + beta.*cuber;

    xiL = xil0*randn;
    xiR = xir0*randn;

    fstep(2) = fstep(2) + xiL - lambda.*v(2);
    fstep(N+1) = fstep(N+1) + xiR - lambda.*v(N+1);

    newa = fstep./m;


    % calculating v(t + delta_t)
    v = v + 0.5.*(a + newa).*tstep;

    a = newa;
    fcum = fcum + fl.*v;

    
    
end

%disp(fstep);
%disp(fcum);
%disp(j);
%disp(sum(fl)/N);

j = sum(fcum)/(maxStep.*(N-1));
time = toc(t);
dat = [x,v,a];

end



