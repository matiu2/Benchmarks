all: var-st var-mt

var-st: var-st.cpp
	c++ -O6 -march=native -o $@ $<

var-mt: var-mt.cpp
	c++ -O6 -march=native -lpthread -o $@ $<
