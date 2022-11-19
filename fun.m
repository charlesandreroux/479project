


function fun()

numPointsPerDeltaM = 6;
numDeltaMs = 5;
N = 1024;
Tl = 0.8;
Tr = 1.2;
lambda = 0.8;
kconstant = 1;
tstep = 0.05;
maxStep = 1e4;
nDataPoints = 100;
save = false;
plotl = false;

deltaMs = linspace(0, 0.8, numDeltaMs);
betas = logspace(-3,2,numPointsPerDeltaM);
js = zeros(numPointsPerDeltaM, numDeltaMs);
%tls = zeros(numPointsPerDeltaM, nDataPoints);
%trs = zeros(numPointsPerDeltaM, nDataPoints);
%tavgs = zeros(numPointsPerDeltaM, nDataPoints);


col = ["k", "r", "#EDB120", "b", "m"];

for k = 1:numDeltaMs    
    
    for i = 1:numPointsPerDeltaM
    
        [dat, j, m, time, T, Tl1, Tr1] = main(deltaMs(k), betas(i), plotl, 1);
    
    
        js(i,k) = mean(j(ceil(length(j)/2):length(j)));
    
    
    end
    
end
figure;
loglog(betas, js);
lbls = "{\Delta}m=" + deltaMs;
legend(lbls);

end