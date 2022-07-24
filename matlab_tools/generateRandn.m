N=10000000;
fid=fopen('randN.dat','w');
for i=1:N
    x=randn();
    fprintf(fid,"%.12f\n",x);
end
fclose(fid);