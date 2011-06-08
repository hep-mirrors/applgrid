FC = gfortran
FFLAGS =  
QLDIR  = ql
FFDIR  = ff

# -----------------------------------------------------------------------------
# Specify the object files. 

TEST=test.o 


all: test

.PHONY: fflib qllib

.f.o:	
		$(FC) $(F77FLAGS) -c $<  -DREPLACEPATH=\'$(PWD)\'


# -----------------------------------------------------------------------------


test: qllib fflib $(TEST)
	$(FC) $(FFLAGS) -o $@  $(TEST) -L$(QLDIR) -lqcdloop  -L$(FFDIR) -lff 

fflib: 
	make -C $(FFDIR) 

qllib:
	make -C $(QLDIR) 


# -----------------------------------------------------------------------------

clean:
	- rm -f *.o *.s *.prj *~ core 

# -----------------------------------------------------------------------------

# DO NOT DELETE


