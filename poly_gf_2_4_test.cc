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

void calc_syndrome(
	poly *syndrome,
	const poly &msg,
	const std::vector<coeff_type> &generator_roots)
{
	*syndrome = poly();
	for(unsigned i = 0; i < generator_roots.size(); i++) {
		coeff_type s_i = msg.evaluate(generator_roots[i]);
		*syndrome += poly(i, s_i);
	}
}

void berlekamp(poly *sigma_out, const poly &syndrome)
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

		C *= poly({ 0, 1 });
		K++;
	}

	*sigma_out = sigma;
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
	
	poly generator = poly::from_roots(generator_roots);
	printf("generator:\n");
	poly_print(generator);

	//poly message({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1});
	poly message({9, 8, 7, 6, 5, 4, 3, 2, 1});
	printf("message:\n");
	poly_print(message);

	poly msg_shift = message << (generator.terms()-1);
	printf("msg shift:\n");
	poly_print(msg_shift);

	poly remainder = msg_shift % generator;
	printf("remainder:\n");
	poly_print(remainder);

	poly transmit = msg_shift + remainder;
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

	// error loc: chien search
	poly err_eval_poly;
	std::vector<unsigned> err_locs;
	for(unsigned i = 0; i < coeff_type::order()-1; i++) {
		coeff_type test_root = lut_singleton.exp(-(int)i);
		coeff_type val = sigma.evaluate(test_root);
		err_eval_poly += poly(i, val);
		if(val == 0) {
			err_locs.push_back(i);
		}
	}
	printf("err_eval_poly:\n");
	poly_print(err_eval_poly);
	printf("err locs:\n");
	for(unsigned i = 0; i < err_locs.size(); i++)
		printf("%u ", err_locs[i]);
	printf("\n");


	// calculate error value polynomial

	// calculate derivative of sigma
	poly sigma_deriv;
	for(unsigned i = 1; i < sigma.terms(); i += 2) {
		sigma_deriv += poly(i-1, sigma[i]);
	}
	printf("sigma_deriv:\n");
	poly_print(sigma_deriv);

	// calculate omega
	poly omega = sigma * syndrome_bad % poly(generator_roots.size(), 1);
	printf("omega:\n");
	poly_print(omega);

	// forney algorithm
	poly correction;
	for(unsigned i = 0; i < err_locs.size(); i++) {
		int loc = (int) err_locs[i];
		coeff_type root = lut_singleton.exp(loc);
		coeff_type root_inv = lut_singleton.exp(-loc);
		// FIXME: fixup root if initial power of gen. poly is not 0
		coeff_type val = root * omega.evaluate(root_inv) / sigma_deriv.evaluate(root_inv);
		correction += poly(loc, val);
	}
	printf("correction:\n");
	poly_print(correction);

	return 0;
}

#if 0
	poly R(generator_roots.size(), 1);
	poly S = syndrome_bad;
	poly A = poly(0, 1);
	poly B;
	int i = 0;
	while(S.terms() >= generator_roots.size() / 2) {
		poly Q = R / S;
		poly S_ = R % S;
		poly A_ = Q * A + B;
		poly R_ = S;
		poly B_ = A;
		i++;
		if(S.terms() >= generator_roots.size() / 2)
			break;
		S = S_;
		A = A_;
		R = R_;
		B = B_;
	}
	poly lambda = A / poly(0, A.term(0));
	poly omega = S / poly(0, A.term(0));

#elif 0
	poly lambda, F, omega;
	// lambda(x) * S(x) + F(x) * x^2t = omega(x)
	poly::polynomial_ext_euclid(
		&omega,
		&F,
		&lambda,
		poly(generator_roots.size(), 1),
		syndrome_bad);

	printf("lambda/sigma:\n");
	poly_print(lambda);

	printf("omega/A:\n");
	poly_print(omega);

	printf("F/B:\n");
	poly_print(F);

	poly sigma_r = lambda.reversed();
	printf("sigma_r:\n");
	poly_print(sigma_r);
#endif

#if 0
	// error values
	std::vector<coeff_type> err_vals;
	for(unsigned i = 0; i < err_locs.size(); i++) {
		coeff_type zi(lut_singleton.exp(err_locs[i]));

		coeff_type zi_inv = coeff_type(1) / zi;
		coeff_type numer = omega.evaluate(zi_inv);

		coeff_type prod = 1;
		for(unsigned j = 0; j < err_locs.size(); j++) {
			if(j == i)
				continue;
			coeff_type zj(lut_singleton.exp(err_locs[j]));
			coeff_type t = zj / zi;
			prod *= t + coeff_type(1);
		}
		coeff_type denom = zi * prod;

		coeff_type yi = numer / denom;
		err_vals.push_back(yi);
	}

	poly errors_est;
	for(unsigned i = 0; i < err_locs.size(); i++) {
		errors_est += poly(err_locs[i], err_vals[i]);
	}

	printf("errors estimated:\n");
	poly_print(errors_est);

	poly corrected = received + errors_est;
	printf("corrected:\n");
	poly_print(corrected);
	printf("transmit compare:\n");
	poly_print(transmit);
	
#endif
