#include <stdio.h>
#include <assert.h>

#include "gf_2_4.h"

int main()
{
	printf("from index to poly to index:\n");
	for(unsigned i = 0; i < 16; i++)
		printf("%u -> %u -> %u\n",
		       i, lut_singleton.from_index(i),
		       lut_singleton.to_index(lut_singleton.from_index(i)));

	printf("\nlogarithms:\n");
	for(unsigned i = 0; i < 16; i++)
		printf("log(%u) = %u\n", i, lut_singleton.log(i));

	printf("\nexponentials:\n");
	for(unsigned i = 0; i <= 28; i++)
		printf("exp(%u) = %u\n", i, lut_singleton.exp(i));

	printf("\naddition:\n");
	for(unsigned i = 0; i < 16; i++) {
		for(unsigned j = 0; j < 16; j++) {
			gf_2_4 result = gf_2_4(i) + gf_2_4(j);
			printf("%3u ", (unsigned) result);
		}
		printf("\n");
	}

	printf("\nmultiplication:\n");
	for(unsigned i = 0; i < 16; i++) {
		for(unsigned j = 0; j < 16; j++) {
			gf_2_4 result = gf_2_4(i) * gf_2_4(j);
			printf("%3u ", (unsigned) result);
		}
		printf("\n");
	}

	printf("\ndivision:\n");
	for(unsigned i = 0; i < 16; i++) {
		for(unsigned j = 1; j < 16; j++) {
			gf_2_4 result = gf_2_4(i) / gf_2_4(j);
			printf("%3u ", (unsigned) result);
		}
		printf("\n");
	}

	// validate inverses
	for(unsigned i = 0; i < 16; i++) {
		for(unsigned j = 0; j < 16; j++) {
			gf_2_4 x(i);
			gf_2_4 y(j);

			// additive
			gf_2_4 z = x + y;
			gf_2_4 w = z - y;
			assert(w == x);

			if(j == 0)
				continue;

			// multiplicative
			z = x * y;
			w = z / y;
			assert(w == x);
		}
	}
	

	return 0;
}
