N=10000000;
fid=fopen('randU.dat','w');
for i=1:N
    x=rand();
    fprintf(fid,"%.12f\n",x);
end
fclose(fid);