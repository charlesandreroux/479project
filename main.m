

function [dat,j,m,time,T, Tl1, Tr1] = main(deltaM, beta, plotl, save, kB)

t = tic();

if ~exist('N', 'var')
    N = 1024;
end
if ~exist('Tl', 'var')
    Tl = 0.8;
end
if ~exist('Tr', 'var')
    Tr = 1.2;
end
if ~exist('lambda', 'var')
    lambda = 0.8;
end
if ~exist('k', 'var')
    k = 1;
end
if ~exist('tstep', 'var')
    tstep = 0.05;
end
if ~exist('deltaM', 'var')
    deltaM = 0.2;
end
if ~exist('maxStep', 'var')
    maxStep = 1e6;
end
if ~exist('beta', 'var')
    beta = 1;
end
if ~exist('nDataPoints', 'var')
    nDataPoints = 100;
end
if ~exist('save', 'var')
    save = true;
end
if ~exist('plotl', 'var')
    plotl = true;
end
if ~exist('kB', 'var')
    kB = 1.380649e-23;
end

nPoints = maxStep/nDataPoints;

% arrays of positions, velocities, and accelerations
x = zeros(N + 2, 1);
v = zeros(N + 2, 1);
a = zeros(N + 2, 1);

% chi terms
xil0 = sqrt(2.*lambda.*kB.*Tl./tstep);
xir0 = sqrt(2.*lambda.*kB.*Tr./tstep);

% array of masses
m = 1 - deltaM + (2*deltaM)*rand(N + 2, 1);
m(1) = inf;
m(N+2) = inf;

% initial conditions
%v(2) = sqrt(kB*Tl/m(2));
%v(N+1) = sqrt(kB*Tr/m(N+2));
%v(2) = Tl;
%v(N+1) = Tr;
% pre-allocation
newa = zeros(N+2,1);
xl = zeros(N+2,1);
xr = zeros(N+2,1);
fl = zeros(N+2,1);
fstep = zeros(N+2,1); 
cubel = zeros(N+2,1);
cuber = zeros(N+2,1);

fcum = zeros(N+2, 1);

j = zeros(nDataPoints,1);
T = zeros(nDataPoints,1);
Tl1 = zeros(nDataPoints,1);
Tr1 = zeros(nDataPoints,1);

tcuml = 0;
tcumr = 0;
tcumall = zeros(N,1);

for i = 1:nDataPoints

    for step = 1:nPoints
    
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
     
        fl = k.*(xl-x) + beta.*cubel;
        fstep = fl + k.*(xr - x) + beta.*cuber;
    
        xiL = xil0*randn;
        xiR = xir0*randn;
    
        fstep(2) = fstep(2) + xiL - lambda.*v(2);
        fstep(N+1) = fstep(N+1) + xiR - lambda.*v(N+1);
    
        newa = fstep./m;
    
    
        % calculating v(t + delta_t)
        v = v + 0.5.*(a + newa).*tstep;
    
        a = newa;
        fcum = fcum + fl.*v;
    
        tcuml = tcuml + v(2)*v(2);
        tcumr = tcumr + v(N+1)*v(N+1);

        if i == 20
            tcumall = tcumall + v(2:N+1).*v(2:N+1);
        end
        
    end

    m(1) = 0;
    m(N+2) = 0;
    j(i) = sum(fcum)./(nPoints.*(N-1));
    T(i) = sum(m.*v.*v)./(kB*N);
    m(1) = Inf;
    m(N+2) = Inf;
    fcum = zeros(N+2, 1);

    Tl1(i) = m(2).*tcuml./(kB.*nPoints);
    Tr1(i) = m(N+1).*tcumr./(kB.*nPoints);
    tcuml = 0;
    tcumr = 0;

end

tcumall = m(2:N+1).*tcumall./(kB.*nPoints);

%disp(fstep);
%disp(fcum);
%disp(j);
%disp(sum(fl)/N);

time = toc(t);
dat = [x,v,a];

nIterations = linspace(0, maxStep, nDataPoints);

j = j/kB;

%{
disp(T)
disp(m)
disp(v)
disp(kB)
disp(N)
%}
if plotl
    plotTitle = "{\Delta}m = " + deltaM + ", {\beta} = " + beta;
    
    plotT = figure();
    plot(nIterations, T);
    hold on
    plot(nIterations, Tl1);
    plot(nIterations, Tr1);
    hold off
    xlabel('Number of Iterations');
    ylabel('Temperature');
    title(plotTitle);
    legend('Average Temperature', 'Left endpoint', 'Right endpoint');
    
    
    plotJ = figure();
    plot(nIterations, j);
    xlabel('Number of Iterations');
    ylabel('Heat Current {j}');
    title(plotTitle);

    plotTall = figure();
    plot(tcumall);
    xlabel('Index of Mass');
    ylabel('Temperature');
    title(plotTitle);
end

if save
    saveas(plotT, "figures/temperature" + "dm" + deltaM + "beta" + beta + ".png");
    saveas(plotJ, "figures/heatcurrent" + "dm" + deltaM + "beta" + beta +".png");
    saveas(plotTall, "figures/alltemp20" + "dm" + deltaM + "beta" + beta +".png");
end


end





