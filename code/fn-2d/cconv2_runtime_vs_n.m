clear;clc;close all;
N= 40;
runtime= zeros(N,1);
for n = 1: N
    f=magic(n); h =magic(n)';
    tic
    w = cconv2(f,h);
    runtime(n)=toc;
    fprintf(1,'n: %d, time: %f\n',n,runtime(n));
end
plot(1:N,runtime);
xlabel('filter size')
ylabel('running time for cconv2')
