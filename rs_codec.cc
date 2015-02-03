#include <stdio.h>
#include <unistd.h>

#include "rs_codec.hpp"
#include "gf_2_8.h"

typedef gf_2_8 coeff_type;
int main(int argc, char *argv[])
{
	FILE *infile, *outfile;
	unsigned t = 0;
	char *endp;
	bool decode = false;

	int opt;
	while((opt = getopt(argc, argv, "di:o:t:")) != -1) {
		switch(opt) {
		case 'd':
			decode = true;
			break;
		case 'i':
			infile = fopen(optarg, "rb");
			if(! infile) {
				perror("fopen(infile)");
				return -1;
			}
			break;
		case 'o':
			outfile = fopen(optarg, "wb");
			if(! outfile) {
				perror("fopen(outfile)");
				return -1;
			}
			break;
		case 't':
			t = strtol(optarg, &endp, 0);
			if(*endp || t > 63) {
				fprintf(stderr, "t must be <= 63\n");
				return -1;
			}
		}
	}

	if(!infile) {
		fprintf(stderr, "-i <infile> required\n");
		return -1;
	}
	if(!outfile) {
		fprintf(stderr, "-o <outfile> required\n");
		return -1;
	}

	std::vector<coeff_type> inbuf;
	for(unsigned i = 0; i < coeff_type::order() - 1 - (2 * t); i++) {
		int ch = fgetc(infile);
		if(ch == EOF) break;
		inbuf.push_back(coeff_type((unsigned)ch));
	}

	polynomial<coeff_type> input(inbuf);
	polynomial<coeff_type> output;

	std::vector<coeff_type> generator_roots;
	for(unsigned i = 0; i < 2 * t; i++)
		generator_roots.push_back(coeff_type::exp(i));

	if(decode) {
		polynomial<coeff_type> syndrome;
		rs_calc_syndrome(&syndrome, input, generator_roots);

		polynomial<coeff_type> sigma;
		rs_berlekamp(&sigma, syndrome);

		std::vector<unsigned> err_locs;
		rs_chien_search(&err_locs, coeff_type::order()-1, sigma);

		polynomial<coeff_type> correction;
		rs_forney(&correction, sigma, syndrome, err_locs, generator_roots);

		output = (correction + input) >> (size_t)(2U * t);
	} else {
		rs_encode(&output, input, generator_roots);
	}

	for(unsigned i = 0; i < output.terms(); i++) {
		unsigned char byte = (unsigned) output[i];
		fwrite(&byte, 1, 1, outfile);
	}

	return 0;
}
