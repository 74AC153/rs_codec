#include "gf_2_4.h"

static unsigned poly_form_from_idx[16] = { 0 };
static unsigned log_poly_form[16] = { 0 };

gf_2_4_lut lut_2_4_singleton;

gf_2_4_lut::gf_2_4_lut()
{
	// field generator polynomial: a^4 = a + 1
	poly_form_from_idx[0] = 0;
	poly_form_from_idx[1] = 1;
	for(unsigned i = 2; i < 16; i++) {
		poly_form_from_idx[i] = poly_form_from_idx[i-1] << 1;
		if(poly_form_from_idx[i] & 0x10)
			poly_form_from_idx[i] ^= 0x13;
	}

	for(unsigned i = 1; i < 16; i++) {
		unsigned poly = poly_form_from_idx[i];
		log_poly_form[poly] = i - 1;
	}
}

gf_2_4_lut::~gf_2_4_lut()
{
}

unsigned gf_2_4_lut::from_index(unsigned index)
{
	return poly_form_from_idx[index];
}

unsigned gf_2_4_lut::to_index(unsigned rawval)
{
	if(rawval == 0)
		return 0;
	return log_poly_form[rawval] + 1;
}

unsigned gf_2_4_lut::log(unsigned rawval)
{
	return log_poly_form[rawval];
}

unsigned gf_2_4_lut::exp(int lograwval)
{
	if(lograwval < 0)
		lograwval += 15;
	if(lograwval > 14)
		lograwval -= 15;
	return poly_form_from_idx[lograwval+1];
}
