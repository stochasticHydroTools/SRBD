clear
clc
format long

%Initialize Parameters and allocate arrays
k = 1; mu = 1; d = 1; tau = 1 / (k* mu); 
N = 16;
jmax = 1000;
jump = 10;
Diffusion = zeros(jmax,N);   
temp = zeros(jmax,N);
t = 1:jump:jmax;
t = t*tau;


% Norm data from each file, and put the norm squared for each file in a
% seperate column. 
for k = 1:N                             % Iterate over all files
    A = importdata(append('diffusion', num2str(k) ,'.txt'))';
    counter = 1;
    for j = 1:jump:jmax
        sum = 0;
        for p = 1:length(A) - j         %Each p is a different t0. Want to average amongst these. 
           sum = sum + norm(A(:,p+j) - A(:,p))^2;
        end
        temp(counter,k) = sum ./ (length(A) - j);      %Average over all t0's.
        Diffusion(counter,k) = temp(counter,k) ./ (2*d*j*tau);  
        counter = counter + 1;
    end
end

means = zeros(jmax/jump,1);
standardDevs = zeros(jmax/jump,1);
ers = zeros(jmax/jump,1);
ts = tinv([0.025 0.975], N-1);          %This will tell me what the t value is for the 95% 

for m = 1:jmax/jump 
   means(m) = mean(temp(m,:));
   standardDevs(m) = std(temp(m,:));
   ers(m) = standardDevs(m)*ts(2)/sqrt(N);
end

errorbar(t,means,ers,'k-','Linewidth',1.5)
xlabel('n*dt')
ylabel('<r_{cm}^2>')

p = polyfit(t',means,1);   %Find slope
D_eff = p(1) / 2*d;

disp('The approximate effective diffusive coefficient is: ')
disp(D_eff)


