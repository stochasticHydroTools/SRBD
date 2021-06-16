% Program to read in the relDist files 1-16 (corresponding to seeds 5-20)
% and create a histogram with error bars. The error bars will follow the
% formula (for each bin) of: (mean +- ts*s/4)
clear; clc; 
format long

N = 16;
ts = tinv([0.025 0.975], N-1);      %t distribution values for the 95% interval. Note this is -2.131, not 1.92


%Read in all data. (will probably take a while)
A1 = importdata('relDist1.txt')';
A2 = importdata('relDist2.txt')';
A3 = importdata('relDist3.txt')';
A4 = importdata('relDist4.txt')';
A5 = importdata('relDist5.txt')';
A6 = importdata('relDist6.txt')';
A7 = importdata('relDist7.txt')';
A8 = importdata('relDist8.txt')';
A9 = importdata('relDist9.txt')';
A10 = importdata('relDist10.txt')';
A11 = importdata('relDist11.txt')';
A12 = importdata('relDist12.txt')';
A13 = importdata('relDist13.txt')';
A14 = importdata('relDist14.txt')';
A15 = importdata('relDist15.txt')';
A16 = importdata('relDist16.txt')';


binWidth = 0.1;  %This follows the general rule of N^(-1/5), as N = 100,000
binEdges = -3.8:binWidth:3.8;    %-3.8 - 3.8 chosen so that I all my data

% Instantiate all histogram objects.
figure(1)
h1 = histogram(A1(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
hold on
h2 = histogram(A2(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h3 = histogram(A3(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h4 = histogram(A4(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h5 = histogram(A5(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h6 = histogram(A6(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h7 = histogram(A7(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h8 = histogram(A8(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h9 = histogram(A9(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h10 = histogram(A10(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h11 = histogram(A11(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h12 = histogram(A12(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h13 = histogram(A13(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h14 = histogram(A14(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h15 = histogram(A15(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 
h16 = histogram(A16(1000:end),'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges); 


%Each Column corresponds to a different bin. There are 76 bins total (by
%design)
myData = [h1.Values; h2.Values; h3.Values; h4.Values; h5.Values; h6.Values; h7.Values; h8.Values; h9.Values; h10.Values; h11.Values; h12.Values; h13.Values; h14.Values; h15.Values; h16.Values;];

numBins = length(h1.Values); %Should be equal to all other histograms, so we just take the first

means = zeros(numBins,1);
standardDevs = zeros(numBins,1);
CI = zeros(numBins,2);
ers = zeros(numBins,1);

%Get mean and STD for each bin.
for k = 1:numBins
    means(k) = mean(myData(:,k));
    standardDevs(k) = std(myData(:,k));
    CI(k,:) = [means(k) + standardDevs(k)*ts(1)/sqrt(N), means(k) + standardDevs(k)*ts(2)/sqrt(N)];
    ers(k) = standardDevs(k)*ts(2)/sqrt(N);     %May be 2 not 1
end

BinCenters = (h1.BinEdges(2:end) + h1.BinEdges(1:end-1))/2;

errorbar(BinCenters,means,ers,'k-','Linewidth',1.25)


figure(2)
plot(BinCenters,means,'Linewidth',1.25)

figure(3)
errorbar(BinCenters,means,ers,'k-','Linewidth',1.5)
hold on
hTotal = histogram([A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16],'Normalization','pdf','Binwidth',binWidth,'BinEdges',binEdges);
xlabel('Displacement r1 - r2')
ylabel('Probability')
x = linspace(-4,4,10001);
y = (1/sqrt(2*pi))*exp((-1*x.^2)./ (2*1));
plot(x,y,'r','Linewidth',1.5)
legend('Mean particle displacement','Histogram of Displacement','Gibbs-Boltzmann Distribution')










