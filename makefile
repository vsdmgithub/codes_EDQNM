# MAKEFILE FOR EULER

# COMPILER
cc=gfortran

# PROGRAM
program=code_EDQNM.f90

# MODULES
timer_mod            =modules-secondary/system_timer.f90
constants_mod        =modules-secondary/system_constants.f90
auxilaries_mod       =modules-secondary/system_auxilaries.f90
basicvariables_mod   =modules-primary/system_basicvariables.f90
basicoutput_mod      =modules-primary/system_basicoutput.f90
initialcondition_mod =modules-primary/system_initialcondition.f90
basicfunctions_mod   =modules-primary/system_basicfunctions.f90
solvereqns_mod       =modules-secondary/system_solver_eqns.f90
solver_mod           =modules-secondary/system_solver.f90
advfunctions_mod     =modules-primary/system_advfunctions.f90
main_mod             =modules-primary/system_main.f90

# OBJECTS
obj=system_timer.o\
	system_constants.o\
	system_auxilaries.o\
	system_basicvariables.o\
	system_basicoutput.o\
	system_initialcondition.o\
	system_basicfunctions.o\
	system_solver_eqns.o\
	system_solver.o\
	system_advfunctions.o\
	system_main.o

# EXECUTABLE
run=./ex
# ex_DIM_GHD
# CLEAN COMMANDS
rmex=rm ex

mkcl=make cl
#----------------------------end-------

# MAKEFILE
# ---------------------------start-----
ex:$(ob)
	$(cc) -c $(timer_mod)
	$(cc) -c $(constants_mod)
	$(cc) -c $(auxilaries_mod)
	$(cc) -c $(basicvariables_mod)
	$(cc) -c $(basicoutput_mod)
	$(cc) -c $(initialcondition_mod)
	$(cc) -c $(basicfunctions_mod)
	$(cc) -c $(solvereqns_mod)
	$(cc) -c $(solver_mod)
	$(cc) -c $(advfunctions_mod)
	$(cc) -c $(main_mod)
	$(cc) $(program) $(obj) -o ex_P_XXX_YYY
	$(mkcl)
	# $(run)

#----------------------------end-------

# CLEANING
# ---------------------------start-----
clean:
	rm ex
	clear
cl:
	rm *.mod
	rm *.o
#----------------------------end-------
