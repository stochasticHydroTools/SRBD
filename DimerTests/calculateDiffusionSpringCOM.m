clear;
clc;
format long
%File to calculate the diffusion coefficient, in this case for two
%particles joined by a spring as well as histogram relative distance
%THe files I am currently reading in are in 1D and have mu = D = 1

d = 1;                          %Dimension
mu = 1.0;
k = 1.0;
l0 = 0;

A = importdata('diffusion.txt')';
tau = 1 / (k* mu);        
jmax = 100;
Diffusion = zeros(jmax,1); 
temp = zeros(jmax,1);

for j = 1:jmax
    sum = 0;
    for p = 1:length(A) - j
       sum = sum + norm(A(:,p+j) - A(:,p))^2;
    end
    temp(j) = sum ./ (length(A) - j);
    Diffusion(j) = temp(j) ./ (2*d*j*tau);  
end

B = importdata('relDist.txt');

C = zeros(length(B),1);

% If in 3D, will need the norm, if in 1D, then only histogram B and not C
for f = 1:length(B)
    C(f) = norm(B(f,:));    
end 
 

disp('As j increases, we are sampling less points so best approx for diffusion coeff is j = 1.')
disp('For j = 1 we have:  ')
disp(Diffusion(1))

kb = physconst('Boltzmann'); T = 300;
disp('Compare to the diffusive coefficient divided by two:')
disp(1 / 2)

figure(1)
%If 1D, histogram B, if 3D, histogram C (this corresponds to displacement
%vs distance)
h = histogram(B(1000:end),'Normalization','pdf','Binwidth',0.1); 
xlabel('Distance between particles r1, r2')
ylabel('Probability')

hold on
x = linspace(-5,5,10001);
y = (1/(sqrt(2*pi*1/k)))*exp((-k*(x-l0).^2)./ (2*1));  %Will be the boltzmann distribution ( D = 1)
plot(x,y,'Linewidth',2)

% To find confidence Intervals....
N = length(B);
s = std(B);     %Sample SD
SE = s / sqrt(N); %Standard Error
ts = tinv([0.025 0.975], N-1); %Returns t-values for this specific case with N-1 DOF. (95% CI)
CI = mean(B) + ts*SE;

% Graph the confidence bands
xline(CI(1), 'Linewidth',1.25)
xline(CI(2),'Linewidth',1.25)

figure(2)
t = 1:1:jmax;
t = t*tau;

plot(t,temp,'o','Linewidth',2)
xlabel('n*dt')
ylabel('<r_{cm}^2>')


 

