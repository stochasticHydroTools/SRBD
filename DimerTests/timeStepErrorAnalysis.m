close all; clear; clc;

%Will look at both weak and strong convergence of EM method. With CI.

N = 16; %Number of seeds
M = 7; % Number of time steps

for j = 1:M
   t(j) = 2^(-j);
end

ts = tinv([0.025 0.975], N-1); %The t value for 95% confidence

%  Read in all data. 
data_timeStep_1 = readSeedData(N,'relDist_dt1_seed');
data_timeStep_2 = readSeedData(N,'relDist_dt2_seed');
data_timeStep_3 = readSeedData(N,'relDist_dt3_seed');
data_timeStep_4 = readSeedData(N,'relDist_dt4_seed');
data_timeStep_5 = readSeedData(N,'relDist_dt5_seed');
data_timeStep_6 = readSeedData(N,'relDist_dt6_seed');
data_timeStep_7 = readSeedData(N,'relDist_dt7_seed');

allData = [data_timeStep_1;data_timeStep_2;data_timeStep_3;data_timeStep_4;data_timeStep_5;data_timeStep_6;data_timeStep_7];

means = zeros(M,1);
standardDevs = zeros(M,1);
ers = zeros(M,1);

for j = 1:M
   means(j) = mean(allData(j,:));
   standardDevs(j) = std(allData(j,:));
   ers(j) = ts(1)*standardDevs(j)/sqrt(N);
end

figure(2)
t = t'; 
loglog(t,means,'ok','Linewidth',1.25)
xlabel('\Delta t')
ylabel('E(\Delta t)')
hold on

p = polyfit(log(t),log(means),1);
p = [0.5,p(2)];
plot(t,exp(polyval(p,log(t))),'r','Linewidth',1.25)
errorbar(t,means,ers,'k-','Linewidth',1.5,'LineStyle','none')
hold off


% Function which makes a matrix for one time step where each column is a
% different seed, outputs mean error across seeds and standard dev. 
function timeStepData = readSeedData(N,fileRoot)
    close all;
    binWidth = 0.1;
    binEdges = -5:binWidth:5;
    
    timStepData = zeros(N,1);
        
    for seed = 5:20
        str = strcat(fileRoot,num2str(seed));
        str = strcat(str,'.txt');
        A = importdata(str);  %Create file name
        
        % Create histogram from a specific seed. 
        h1 = histogram(A(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges);
        
     
        timeStepData(seed - 4) = error(h1,A);
  
      
    end
    
    

end

%Error function for STRONG convergence.
function err = error(hist,TS)
    f = @(x) (1/sqrt(2*pi))*exp((-1*x.^2)./(2*1));
    means = hist.Values;
    BinCenters = (hist.BinEdges(2:end) + hist.BinEdges(1:end-1))/2;
    actual = f(BinCenters);
    %err = abs(std(TS) - 1);
    err = mean(abs(means-actual));   \
    
    % Error for WEAK convergence.
    %err = abs( mean(TS) );

end