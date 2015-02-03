#include "gf_2_8.h"

static unsigned poly_form_from_idx[256] = { 0 };
static unsigned log_poly_form[256] = { 0 };

gf_2_8_lut lut_2_8_singleton;

gf_2_8_lut::gf_2_8_lut()
{
	// field generator polynomial: x^8 + x^4 + x^3 + x^2 + 1 --> 0x11d
	poly_form_from_idx[0] = 0;
	poly_form_from_idx[1] = 1;
	for(unsigned i = 2; i < 256; i++) {
		poly_form_from_idx[i] = poly_form_from_idx[i-1] << 1;
		if(poly_form_from_idx[i] & 0x100)
			poly_form_from_idx[i] ^= 0x11d;
	}

	for(unsigned i = 1; i < 256; i++) {
		unsigned poly = poly_form_from_idx[i];
		log_poly_form[poly] = i - 1;
	}
}

gf_2_8_lut::~gf_2_8_lut()
{
}

unsigned gf_2_8_lut::from_index(unsigned index)
{
	return poly_form_from_idx[index];
}

unsigned gf_2_8_lut::to_index(unsigned rawval)
{
	if(rawval == 0)
		return 0;
	return log_poly_form[rawval] + 1;
}

unsigned gf_2_8_lut::log(unsigned rawval)
{
	return log_poly_form[rawval];
}

unsigned gf_2_8_lut::exp(int lograwval)
{
	if(lograwval < 0)
		lograwval += 255;
	if(lograwval >= 255)
		lograwval -= 255;
	return poly_form_from_idx[lograwval+1];
}
