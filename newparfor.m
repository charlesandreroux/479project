
% param m: the mass distribution
% param desc: optional string argument to describe the mass distribution
%             and save the plots if present

function newparfor(m, desc)

% make sure m is a column vector
m = m(:);

% append an infinite mass to both ends
m(length(m) + 1) = inf;
m = flipud(m);
m(length(m) + 1) = inf;
m = flipud(m);

t = tic();

% parameters for the calculations
N = length(m)-2;
Tl = 0.8;
Tr = 1.2;
lambda = 0.8;
k = 1;
beta = 1;
kB = 1;
tStep = 0.05;
maxStep = 1e6;
nGraphDataPoints = 50;
nPoints = maxStep/nGraphDataPoints;
nIterations = linspace(1, maxStep, nGraphDataPoints);
nPars = 500;

% results with original temperatures (j+)
jps = zeros(nPars,nGraphDataPoints);
Tlps = zeros(nPars,nGraphDataPoints);
Trps = zeros(nPars,nGraphDataPoints);
Tps = zeros(nPars,nGraphDataPoints);
tcumallps = zeros(nPars,N);

% results with swapped temperatures (j-)
jms = zeros(nPars,nGraphDataPoints);
Tlms = zeros(nPars,nGraphDataPoints);
Trms = zeros(nPars,nGraphDataPoints);
Tms = zeros(nPars,nGraphDataPoints);
tcumallms = zeros(nPars,N);

% parallel processing
parfor i = 1:nPars

    [jpi, Tlpi, Trpi, Tpi, tcumallpi] = compute(m,N,Tl,Tr,lambda,beta,k,kB, ...
    tStep,nGraphDataPoints,nPoints);
    jps(i,:) = jpi;
    Tlps(i,:) = Tlpi;
    Trps(i,:) = Trpi;
    Tps(i,:) = Tpi;
    tcumallps(i,:) = tcumallpi;

end

% same, but with temperatures swapped
parfor i = 1:nPars

    [jmi, Tlmi, Trmi, Tmi, tcumallmi] = compute(m,N,Tr,Tl,lambda,beta,k,kB, ...
    tStep,nGraphDataPoints,nPoints);
    jms(i,:) = jmi;
    Tlms(i,:) = Tlmi;
    Trms(i,:) = Trmi;
    Tms(i,:) = Tmi;
    tcumallms(i,:) = tcumallmi;

end

% averaging over parallel processing results
jp = mean(jps(:,1:nGraphDataPoints));
Tlpg = mean(Tlps(:,1:nGraphDataPoints));
Trpg = mean(Trps(:,1:nGraphDataPoints));
Tp = mean(Tps(:,1:nGraphDataPoints));
tcumallp = mean(tcumallps(:,1:N));

jm = mean(jms(:,1:nGraphDataPoints));
Tlmg = mean(Tlms(:,1:nGraphDataPoints));
Trmg = mean(Trms(:,1:nGraphDataPoints));
Tm = mean(Tms(:,1:nGraphDataPoints));
tcumallm = mean(tcumallms(:,1:N));

% take the absolute values to compute f_r accurately
jp = abs(jp);
jm = abs(jm);

% designate j+ as the higher value of the two
swap = false;
if mean(jp) < mean(jm)
    temp = jp;
    jp = jm;
    jm = temp;
    swap = true;
end

%%%%%%%%%% not sure which one is best between
% graph of f_r values
fr = (jp - jm)./jm.*100;
% single f_r value from the mean j's
frm = (mean(jp)-mean(jm))/mean(jm)*100;

% plotting
plotTitle1 = "{T_l} = " + Tl + ", {T_r} = " + Tr;
plotTitle2 = "{T_l} = " + Tr + ", {T_r} = " + Tl;

pltOT = figure();
plot(nIterations, Tp);
hold on
plot(nIterations, Tlpg);
plot(nIterations, Trpg);
hold off
xlabel('Number of Iterations');
ylabel('Temperature');
title(plotTitle1);
legend('Average Temperature', 'Left endpoint', 'Right endpoint');

pltST = figure();
plot(nIterations, Tm);
hold on
plot(nIterations, Tlmg);
plot(nIterations, Trmg);
hold off
xlabel('Number of Iterations');
ylabel('Temperature');
title(plotTitle2);
legend('Average Temperature', 'Left endpoint', 'Right endpoint');

pltJ = figure();
plot(nIterations, jp);
hold on
plot(nIterations, jm)
hold off
xlabel('Number of Iterations');
ylabel('Heat Current ({j})');
lbls = ["Left to Right", "Right to Left"];

% swap the labels if the higher j was in reverse direction
if swap
    lbls = fliplr(lbls);
end

lbls(1) = lbls(1) + " {j_+}";
lbls(2) = lbls(2) + " {j_-}";
legend(lbls);
titj = "$\bar{j_+}$ = " + mean(jp) + ", $\bar{j_-}$ = " + mean(jm)...
    + ", $f_r$ = " + frm;
title(titj, 'interpreter', 'latex');

pltFR = figure();
plot(nIterations, fr);
xlabel('Number of Iterations');
ylabel('Rectification Factor {f_r}');
titfr = "$\bar{f_r}$ = " + mean(fr);
title(titfr, 'interpreter', 'latex');

pltTALL = figure();
plot(tcumallp);
hold on
plot(tcumallm);
hold off
xlabel('Index of Mass');
ylabel('Temperature');
legend(plotTitle1, plotTitle2);

if exist('desc', 'var')
    saveas(pltOT, "rect/originalT_averageTemps_" + desc + ".png");
    saveas(pltST, "rect/swappedT_averageTemps_" + desc + ".png");
    saveas(pltJ, "rect/J_" + desc + ".png");
    saveas(pltFR, "rect/Fr_" + desc + ".png");
    saveas(pltTALL, "rect/initialTemps_" + desc + ".png");
end

toc(t);

end

function [j, Tl1, Tr1, T, tcumall] = compute(m,N,Tl,Tr,lambda, beta, k, kB, ...
    tstep, nDataPoints, nPoints)

% arrays of positions, velocities, and accelerations
x = zeros(N + 2, 1);
v = zeros(N + 2, 1);
a = zeros(N + 2, 1);

% chi terms
xil0 = sqrt(2.*lambda.*kB.*Tl./tstep);
xir0 = sqrt(2.*lambda.*kB.*Tr./tstep);

% pre-allocation
newa = zeros(N+2,1);
xl = zeros(N+2,1);
xr = zeros(N+2,1);

fcum = zeros(N+2, 1);
fcumsingle = 0;

j = zeros(nDataPoints,1);
jsingle = zeros(nDataPoints, 1);
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
    
        xr = circshift(x,-1);
        xr(1) = 0;

        diffl = (xl-x);
        diffr = (xr-x);
     
        % not sure if we're supposed to remove the linear term here
        fstep = k.*(diffl + diffr) + beta.*(diffl.*diffl.*diffl...
            + diffr.*diffr.*diffr);

        xiL = xil0*randn;
        xiR = xir0*randn;
        
        % adding xi and lambda terms to endpoints
        fstep(2) = fstep(2) + xiL - lambda.*v(2);
        fstep(N+1) = fstep(N+1) + xiR - lambda.*v(N+1);
    
        newa = fstep./m;
    
    
        % calculating v(t + delta_t)
        v = v + 0.5.*(a + newa).*tstep;
    
        a = newa;

        vl = circshift(v, 1);
        vl(N+2) = 0;
        sumvl = vl + v;

        % used to calculate moving average of j based on 
        % j_{n+1/2} = 1/2 * (v_n + v_{n+1}) * k * (x_{n+1} - x_n) from our
        % notes from last week, not sure if I interpreted it correctly or
        % which is better between:

        % sum all, and take the average at the end
        fcum = fcum + k/2.*sumvl.*diffl;

        % take the average at each iteration
        fcumsingle = fcumsingle + mean(fcum(3:N+1));
    
        %tcum = tcum + v.*v;
        tcuml = tcuml + v(2)*v(2);
        tcumr = tcumr + v(N+1)*v(N+1);
        
        % used to graph the temperature everywhere at the 20th (arbitrarily
        % small) iteration
        if i == 20
            tcumall = tcumall + v(2:N+1).*v(2:N+1);
        end
        
    end
    
    % calculating j at each data point
    j(i) = sum(fcum(3:N+1))./(nPoints.*(N-1));
    jsingle(i) = fcumsingle/(nPoints);

    m(1) = 0;
    m(N+2) = 0;
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
j = j/kB;

end