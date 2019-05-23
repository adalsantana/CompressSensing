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
amp = 1.2*max(abs(x)); 
figure; 
subplot(3, 1, 1); 
plot(x); 
title('Original Signal');
xlim([1 N]);
ylim([-amp amp]);

%obtain K measurements 
A = rand(K, N);
y = A*x; 
subplot(3, 1, 2);
plot(y); 
title('K measured values');
xlim([1 K]);

%perform compressed sensing recovery 
x0 = A.'*y; 
xp = l1eq_pd(x0, A, [], y);
subplot(3, 1, 3);
plot(real(xp));
title('K Measured Values');
xlim([1 N]);
ylim([-amp amp]);
