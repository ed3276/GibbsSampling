obj = main.o src/normpdf.o
GibbsSampling.exe : $(obj)
	gcc -o GibbsSampling $(obj)
main.o : main.c inc/myStatFunc.h
	gcc -c main.c -o main.o
src/normpdf.o : src/normpdf.c inc/myStatFunc.h
	gcc -c src/normpdf.c -o src/normpdf.o
clean :
	rm GibbsSampling $(obj)
