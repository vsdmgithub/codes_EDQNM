# MAKEFILE FOR EDQNM

# DEFINE VARIABLES
# ---------------------------start-----
prog=edqnm_code.f90
cc=gfortran
obj=edqnm_code.o
cc_lc= -I/usr/local/include
run=./ex
#----------------------------end-------


# MAKEFILE
# ---------------------------start----- 
ex:$(ob)
	$(cc) -c $(prog) 
	$(cc) $(obj) -o ex 
	$(run)
#----------------------------end-------

# CLEANING
# ---------------------------start----- 
clean:
	rm ex
	rm *.mod
	rm *.o
cl:
	rm *.mod
	rm *.o
#----------------------------end-------
