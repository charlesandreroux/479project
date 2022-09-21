

function xva = main(N, Tl, Tr, lambda, k, tstep, deltaM, maxStep, beta)
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
    maxStep = 10e5;
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

%array of masses
m = zeros(N+2, 1);
m = 1 - deltaM + (2*deltaM)*rand(N + 2, 1);
m(1) = inf;
m(N+2) = inf;

xiL = sqrt(2*lambda*kB*Tl);
xiR = sqrt(2*lambda*kB*Tr);


for step = 1:maxStep
    
    newx = zeros(N+2, 1);
    newx = x + v.*tstep + 0.5.*a.*tstep.^2;
    newa = zeros(N+2, 1);

    for i = 2:N+1
        
        % calculating a(t + delta_t)
        newa(i) = k*(newx(i-1) - 2*newx(i) + newx(i+1))...
            + beta.*((newx(i-1) - newx(i))^3 + (newx(i+1) - newx(i)).^3);
        
    end

    newa(2) = newa(2) + xiL - lambda.*v(2);
    newa(N+1) = newa(N+1) + xiR - lambda.*v(N+1);
    newa = newa./m;


    % calculating v(t + delta_t)
    newv = zeros(N+2, 1);
    newv = v + 0.5.*(a + newa).*tstep;
    
    x = newx;
    v = newv;
    a = newa;
    
end

xva = [x,v,a];
toc(t);

end



