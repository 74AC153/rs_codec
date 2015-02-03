#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "soak_stream.h"

void swap(int *x, int *y)
{
	int temp = *x;
	*x = *y;
	*y = temp;
}

void shuffle(int *vals, int len)
{
	for(int i = 0; i < len; i++) {
		int pos = rand() % (len - i);
		swap(vals + pos, vals + len - i - 1);
	}
}


int main(int argc, char *argv[])
{
	unsigned int s = 0;
	int n = -1;
	char *endp;
	FILE *instream = NULL, *outstream = NULL;
	unsigned char *buf;
	size_t buflen;

	int opt;
	while((opt = getopt(argc, argv, "i:o:s:n:")) != -1) {
		switch(opt) {
		case 'i':
			instream = fopen(optarg, "r");
			if(instream == NULL) {
				perror("open(infile)");
				return -1;
			}
			break;
		case 'o':
			outstream = fopen(optarg, "w");
			if(outstream == NULL) {
				perror("fopen(outfile)");
				return -1;
			}
			break;
		case 's':
			s = strtol(optarg, &endp, 0);
			if(*endp) {
				fprintf(stderr, "could not parse s\n");
				return -1;
			}
			break;
		case 'n':
			n = strtol(optarg, &endp, 0);
			if(*endp) {
				fprintf(stderr, "could not parse n\n");
				return -1;
			}
			break;
		default:
			return -1;
		}
	}

	if(instream == NULL) {
		fprintf(stderr, "-i <infile> required\n");
		return -1;
	} else if(outstream == NULL) {
		fprintf(stderr, "-o <outfile> required\n");
		return -1;
	} else if(n < 0) {
		fprintf(stderr, "-n <badbytes> required\n");
		return -1;
	} else if(s == 0) {
		fprintf(stderr, "-s <rand-seed> required and != 0\n");
		return -1;
	}

	int rc = soak_FILE(instream, &buf, &buflen);
	if(rc) {
		return rc;
	}
	
	srand(s);

	int *positions = malloc(buflen * sizeof(*positions));
	for(int i = 0; i < (int) buflen; i++)
		positions[i] = i;

	shuffle(positions, buflen);

	for(unsigned i = 0; i < (unsigned) n; i++) {
		int pos = positions[i];
		unsigned char err_val;
		while(! (err_val = rand())); // must be nonzero corruption
		fprintf(stderr, "%u: %x\n", pos, err_val);

		buf[pos] ^= err_val;	
	}

	free(positions);

	rc = fwrite(buf, 1, buflen, outstream);
	if(rc != (int) buflen) {
		perror("fwrite()");
		return -1;
	}

	return 0;
}
