#include <stdio.h>

#include "polynomial.hpp"
#include "gf_2_4.h"

typedef gf_2_4 coeff_type;
typedef polynomial<coeff_type> poly;

void poly_print(poly &P)
{
	for(unsigned i = 0; i < P.terms(); i++)
		printf("%u ", (unsigned) P[i]);
	printf("( ");
	for(unsigned i = 0; i < P.terms(); i++) {
		if((unsigned) P[i] == 0) printf("- ");
		else printf("%u ", lut_singleton.log((unsigned) P[i]));
	}
	printf(")\n");
}

void rs_encode(
	poly *transmit,
	const poly &message,
	const std::vector<coeff_type> &generator_roots)
{
	poly generator = poly::from_roots(generator_roots);
	printf("generator:\n");
	poly_print(generator);

	poly msg_shift = message << (generator.terms()-1);
	printf("msg shift:\n");
	poly_print(msg_shift);

	poly remainder = msg_shift % generator;
	printf("remainder:\n");
	poly_print(remainder);

	// transmit = msg * x^(2t) % generator
	*transmit = msg_shift + remainder;
}

void calc_syndrome(
	poly *syndrome,
	const poly &msg,
	const std::vector<coeff_type> &generator_roots)
{
	// syndrome[i] = msg.eval(gen_root[i]) * x^i
	*syndrome = poly();
	for(unsigned i = 0; i < generator_roots.size(); i++) {
		coeff_type s_i = msg.evaluate(generator_roots[i]);
		*syndrome += poly(i, s_i);
	}
}

void berlekamp(
	poly *sigma_out,
	const poly &syndrome)
{
	unsigned K = 1;
	unsigned L = 0;
	poly sigma { 1 };
	poly sigma_new;
	poly C { 0, 1 };

	coeff_type e;

	while(K <= syndrome.terms()) {
		e = syndrome[K-1];
		for(unsigned i = 1; i <= L; i++)
			e = e + sigma[i] * syndrome[K - 1 - i];

		if(e != 0) {
			sigma_new = sigma + C * e;
			if(2 * L < K) {
				L = K - L;
				C = sigma / e;
			}
			sigma = sigma_new;
		}

		C *= poly(1, 1);
		K++;
	}

	*sigma_out = sigma;
}

void chien_search(
	std::vector<unsigned> *err_locs,
	unsigned len,
	const poly &sigma)
{
	//poly err_eval_poly;
	for(unsigned i = 0; i < len; i++) {
		coeff_type test_root = lut_singleton.exp(-(int)i);
		coeff_type val = sigma.evaluate(test_root);
		//err_eval_poly += poly(i, val);
		if(val == 0) {
			err_locs->push_back(i);
		}
	}
}

void forney(
	poly *correction,
	const poly &sigma,
	const poly &syndrome,
	const std::vector<unsigned> &err_locs,
	const std::vector<coeff_type> &generator_roots)
{
	// calculate derivative of sigma
	// "set even-powered coefficients to 0, then divide by x"
	poly sigma_deriv;
	for(unsigned i = 1; i < sigma.terms(); i += 2) {
		sigma_deriv += poly(i-1, sigma[i]);
	}
	printf("sigma_deriv:\n");
	poly_print(sigma_deriv);

	// calculate omega
	// omega = (syndrome * sigma) mod x^(2t)
	poly omega = (sigma * syndrome) % poly(generator_roots.size(), 1);
	printf("omega:\n");
	poly_print(omega);

	// forney algorithm
	// err_j = X_j^(1-b) omega(X_j^(-1)) / lambda_deriv(X_j^(-1))
	// X_j = a^(err_loc[j])
	// b = exponent for first generator root
	// NB: assume b is 0 for now (see FIXME below)
	*correction = poly();
	for(unsigned i = 0; i < err_locs.size(); i++) {
		int loc = (int) err_locs[i];
		coeff_type root = lut_singleton.exp(loc);
		coeff_type root_inv = lut_singleton.exp(-loc);
		// FIXME: fixup root if initial power of gen. poly is not 0
		coeff_type val = root * omega.evaluate(root_inv) / sigma_deriv.evaluate(root_inv);
		*correction += poly(loc, val);
	}
}

int main(void)
{
	std::vector<coeff_type> generator_roots {
		lut_singleton.exp(0),
		lut_singleton.exp(1),
		lut_singleton.exp(2),
		lut_singleton.exp(3),
		lut_singleton.exp(4),
		lut_singleton.exp(5),
	};

	//poly message({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1});
	poly message({9, 8, 7, 6, 5, 4, 3, 2, 1});
	printf("message:\n");
	poly_print(message);

	
	poly transmit;
	rs_encode(&transmit, message, generator_roots);
	printf("transmit:\n");
	poly_print(transmit);

	poly syndrome_good;
	calc_syndrome(&syndrome_good, transmit, generator_roots);
	printf("good syndrome:\n");
	poly_print(syndrome_good);
	

	//poly errors_actual({0, 0, 2, 0, 0, 0, 0, 0, 0, 13});
	poly errors_actual({0, 0, 1, 0, 15, 0, 0, 0, 0, 1});
	printf("errors_actual:\n");
	poly_print(errors_actual);

	poly received = transmit + errors_actual;
	printf("received:\n");
	poly_print(received);


	poly syndrome_bad;
	calc_syndrome(&syndrome_bad, received, generator_roots);
	printf("bad syndrome:\n");
	poly_print(syndrome_bad);

	poly sigma;
	berlekamp(&sigma, syndrome_bad);
	printf("sigma:\n");
	poly_print(sigma);

	std::vector<unsigned> err_locs;
	chien_search(&err_locs, coeff_type::order()-1, sigma);
	printf("err locs:\n");
	for(unsigned i = 0; i < err_locs.size(); i++)
		printf("%u ", err_locs[i]);
	printf("\n");

	poly correction;
	forney(&correction, sigma, syndrome_bad, err_locs, generator_roots);
	printf("correction:\n");
	poly_print(correction);

	return 0;
}
