This program used to simulate Gibbs Sampling.

1. Use matlab to generate uniform distribution samples and normal distribution samples.
uniform : matlab_tools/generateRand.m 
normal  : matlab_tools/generateRandn.m


2. execute make in command line to generate GibbsSampling.execute

3. GibbsSampling.exe randU.dat randN.dat randGibbs.dat
randGibbs.dat is output file, contains Gibbs samples.

4. Use matlab to draw result, run matlab_tools/GibbsSamplingByC.m


Note: default sampling size is very big, run it will take a lot of time, you can 
dicrease it to speed up, but results will loss resolution.
dicrease it to speed up, but results will loss resolution.