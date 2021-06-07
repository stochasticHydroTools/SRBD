% Plots the boltzmann distribution, as of right now, the mean of this
% distribution is not centered where I think it should be (l0). 
clear
format long
A = importdata('data4.txt');
h = histogram(A,'Normalization','probability')
xlabel('Distance between particles r1, r2')
ylabel('Probability')

hold on 
k = 0.029; l0 = 1e-10; kb = 1.38e-23; T = 300;
x = linspace(0,2e-9,10001);
y = (1/(sqrt(2*pi*kb*T/k)))*exp((-k*(x-l0).^2)./ (2*kb*T));  %Will be the boltzmann distribution and SHOULD overlap (but doesn't)

%plot(x,y)







