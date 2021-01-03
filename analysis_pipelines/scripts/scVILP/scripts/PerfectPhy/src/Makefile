EXECNAME=perfectphy

fast : 
	g++ -Wall -O3 -o $(EXECNAME) sprop4.cpp

default : 
	g++ -Wall -o $(EXECNAME) sprop4.cpp

debug : 
	g++ -Wall -g -D_GLIBCXX_DEBUG -o $(EXECNAME) sprop4.cpp

clean : 
	rm -f $(EXECNAME)

############### The targets below are for experimental use only! ############### 

runtestmintree :
	g++ -Wall -O3 -DRUNTESTMINTREE=true -o $(EXECNAME)_mintree sprop4.cpp

variety :
	g++ -Wall -O3 -DMISEXT=true -DMISLIMITK=true -DKMAXISBRANCHTYPE=3 -o kw_kmaxis.out sprop4.cpp
	g++ -Wall -O3 -DMISEXT=true -DMISLIMITK=true -DKMAXISBRANCHTYPE=2 -o kw_kmaxis_first.out sprop4.cpp
	g++ -Wall -O3 -DMISEXT=true -DMISLIMITK=true -DKMAXISBRANCHTYPE=1 -o kw_kmaxis_second.out sprop4.cpp
	g++ -Wall -O3 -DMISEXT=true -DMISLIMITK=false -o kw_tsumis.out sprop4.cpp
	g++ -Wall -O3 -DMISEXT=false -DMISLIMITK=false -o kw_brute.out sprop4.cpp
	g++ -Wall -O3 -DMISEXT=false -DMISLIMITK=true -o kw_kbrute.out sprop4.cpp

gcov :
	g++ -Wall -g -fprofile-arcs -ftest-coverage -o $(EXECNAME) sprop4.cpp
	./$(EXECNAME) $(ARGS)
	gcov sprop4.cpp

clean_gcov : 
	rm -f *.gcov *.gcda *.gcno

gprof :
	g++ -Wall -pg -o $(EXECNAME) sprop4.cpp
	./$(EXECNAME) $(ARGS)
	gprof $(EXECNAME) > gprof_output.txt

clean_gprof : 
	rm -f gmon.out gprof_output.txt
