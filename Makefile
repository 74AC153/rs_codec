OUTPUT = poly_test gf_2_4_test poly_gf_2_4_test rs_codec corrupt

default: ${OUTPUT}

poly_test: poly_test.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

gf_2_4_test: gf_2_4_test.cc gf_2_4.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

poly_gf_2_4_test: poly_gf_2_4_test.cc gf_2_4.cc gf_2_8.cc gf_2_16.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

rs_codec: rs_codec.cc gf_2_8.cc
	g++ -g3 -Wall -Wextra -std=c++11 -o $@ $^

corrupt: corrupt.c
	gcc -g3 -Wall -Wextra -std=c99 -D_POSIX_C_SOURCE=2 -o $@ $^

clean:
	-rm -f ${OUTPUT}
