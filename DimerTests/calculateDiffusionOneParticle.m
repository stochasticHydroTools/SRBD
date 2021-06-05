%File to calculate the diffusion coefficient, in this case for a single
%brownian particle
% In the fortran code, I specified D = 9 (so i expect this to be returned)
% Also sampled 10,000 brownian iterates.

format long
A = readmatrix('data1.txt')';
sum = 0;
tau = 0.001;        %Arbitrarily chosen, but must be consistent with what is in fortran code
jmax = 100;
Diffusion = zeros(jmax,1);       

for j = 1:jmax
    sum = 0;
    for p = 1:length(A) - j
        sum = sum + norm(A(:,p+j) - A(:,p))^2;
    end
    Diffusion(j) = sum ./ (length(A) - j);
    Diffusion(j) = Diffusion(j) ./ (6*j*tau);
end
disp('As j increases, we are sampling less points (j = 100 means only 100 points sampled) so best approx for diffusion coeff is j = 1.')
disp('For j = 1 we have:  ')
disp(Diffusion(1))
disp('Approximately 9, as expected.')


 
 

