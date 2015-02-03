#include "gf_2_16.h"

static unsigned poly_form_from_idx[65536] = { 0 };
static unsigned log_poly_form[65536] = { 0 };

gf_2_16_lut lut_2_16_singleton;

gf_2_16_lut::gf_2_16_lut()
{
	// FIXME:
	// field generator polynomial: x^16 + x^12 +x^3 +x + 1 --> 0x1100b
	poly_form_from_idx[0] = 0;
	poly_form_from_idx[1] = 1;
	for(unsigned i = 2; i < 65536; i++) {
		poly_form_from_idx[i] = poly_form_from_idx[i-1] << 1;
		if(poly_form_from_idx[i] & 0x10000)
			poly_form_from_idx[i] ^= 0x1100b;
	}

	for(unsigned i = 1; i < 65536; i++) {
		unsigned poly = poly_form_from_idx[i];
		log_poly_form[poly] = i - 1;
	}
}

gf_2_16_lut::~gf_2_16_lut()
{
}

unsigned gf_2_16_lut::from_index(unsigned index)
{
	return poly_form_from_idx[index];
}

unsigned gf_2_16_lut::to_index(unsigned rawval)
{
	if(rawval == 0)
		return 0;
	return log_poly_form[rawval] + 1;
}

unsigned gf_2_16_lut::log(unsigned rawval)
{
	return log_poly_form[rawval];
}

unsigned gf_2_16_lut::exp(int lograwval)
{
	if(lograwval < 0)
		lograwval += 65535;
	if(lograwval >= 65535)
		lograwval -= 65535;
	return poly_form_from_idx[lograwval+1];
}
