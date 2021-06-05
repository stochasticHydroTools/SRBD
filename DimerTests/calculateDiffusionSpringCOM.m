%File to calculate the diffusion coefficient, in this case for two
%particles joined by a spring
% In the fortran code, I specified D = 9 for both r1,r2 (these should
% really be specified by einsteins law, but i found D was too small, and I
% would just get 0 from this code, so I am trying to brute force 9 for the
% time being.
% Also sampled 10,000 brownian iterates.

format long
A = readmatrix('data3.txt')';
sum = 0;
tau = 6*pi*0.2*1 / 10;        %Arbitrarily chosen, but must be consistent with what is in fortran code
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
disp('Approximately 4.5, if this is true, must show it analytically.')

%Probably not correct, will have to review this, doesn't make sense why
%this would be as diffusive as the 2 particles no spring case.


 
 

