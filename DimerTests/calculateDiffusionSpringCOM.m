clear;
clc;
format long
%File to calculate the diffusion coefficient, in this case for two
%particles joined by a spring as well as histogram relative distance
%THe files I am currently reading in are in 1D and have mu = D = 1

d = 1;                          %Dimension
mu = 1.0;
k = 0.029;
l0 = 0;

A = importdata('diffusion.txt')';
tau = 1 / (k* mu);        
jmax = 100;
Diffusion = zeros(jmax,1);   

for j = 1:jmax
    sum = 0;
    for p = 1:length(A) - j
       sum = sum + norm(A(:,p+j) - A(:,p))^2;
    end
    Diffusion(j) = sum ./ (length(A) - j);
    Diffusion(j) = Diffusion(j) ./ (2*d*j*tau);  
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
h = histogram(B(1000:end),'Normalization','pdf'); 
xlabel('Distance between particles r1, r2')
ylabel('Probability')

hold on
x = linspace(-20,20,10001);
y = (1/(sqrt(2*pi*1/k)))*exp((-k*(x-l0).^2)./ (2*1));  %Will be the boltzmann distribution ( D = 1)
plot(x,y,'Linewidth',2)


    
 

