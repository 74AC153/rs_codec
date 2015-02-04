#if ! defined(SOAK_STREAM_H_DEFINED)
#define SOAK_STREAM_H_DEFINED

#include <stdio.h>
#include <stddef.h>

static int soak_FILE(FILE *stream, unsigned char **buf_out, size_t *len_out)
{
	int rc;
	size_t bytes_read = 0;
	size_t buflen = 4096;
	unsigned char *buf = NULL;

	buf = malloc(buflen);
	if(! buf) {
		return -1;
	}

	while((rc = fread(buf + bytes_read, 1, buflen - bytes_read, stream)) > 0) {
		bytes_read += rc;
		if(bytes_read == buflen) {
			unsigned char *newbuf = realloc(buf, buflen * 2);
			if(! newbuf) {
				goto cleanup;
			}
			buflen *= 2;
			buf = newbuf;
		}
	}
	if(rc < 0) {
		goto cleanup;
	}
	*buf_out = buf;
	*len_out = bytes_read;
	return 0;

cleanup:
	free(buf);
	return -2;
}

#endif
