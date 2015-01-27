OUTPUT = poly_test gf_2_4_test poly_gf_2_4_test

default: ${OUTPUT}

poly_test: poly_test.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

gf_2_4_test: gf_2_4_test.cc gf_2_4.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

poly_gf_2_4_test: poly_gf_2_4_test.cc gf_2_4.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

clean:
	-rm -f ${OUTPUT}
