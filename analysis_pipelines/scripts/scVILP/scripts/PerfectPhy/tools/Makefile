all : convert missingdata newicktodot

convert : convert.cpp deconvert.cpp
	g++ -Wall -O3 -o convert convert.cpp

missingdata : missingdata.cpp
	g++ -Wall -O3 -o missingdata missingdata.cpp

newicktodot : newicktodot.cpp
	g++ -Wall -O3 -o newicktodot newicktodot.cpp
