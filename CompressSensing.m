%initialize constants and variables 
rng(0);         %set RNG seed
N = 256;        %length of the signal 
P = 5;          %number of nonzero peaks 
K = 64;         %number of measurements to take (N < L)
x = zeros(N,1); %original signal (P-sparse)

%generate signal with P randomly spread values 
peaks = randperm(N); 
peaks = peaks(1:P);
x(peaks) = rand(1, P);

