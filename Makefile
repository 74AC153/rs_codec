default: poly_test

poly_test: poly_test.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

clean:
	-rm poly_test
