#include <stdio.h>

#include "polynomial.hpp"
#include "gf_2_4.h"

typedef polynomial<gf_2_4> poly;

void poly_print(poly &P)
{
	for(unsigned i = 0; i < P.terms(); i++)
		printf("%u ", (unsigned) P.term(i));
	printf("( ");
	for(unsigned i = 0; i < P.terms(); i++) {
		if((unsigned) P.term(i) == 0) printf("- ");
		else printf("%u ", lut_singleton.log((unsigned) P.term(i)));
	}
	printf(")\n");
}

void berlekamp(poly *sigma_out, const poly &S)
{
	unsigned K = 1;
	unsigned L = 0;
	poly lambda { 1 };
	poly lambda_new;
	poly C { 0, 1 };

	gf_2_4 e;

	while(K <= S.terms()) {
		e = S.term(K-1);
		for(unsigned i = 1; i <= L; i++)
			e = e + lambda.term(i) * S.term(K - 1 - i);

		if(e != 0) {
			lambda_new = lambda + C * e;
			if(2 * L < K) {
				L = K - L;
				C = lambda / e;
			}
		}

		C *= poly({ 0, 1 });
		lambda = lambda_new;
		K++;
	}

	*sigma_out = lambda;
}

int main(void)
{
	std::vector<gf_2_4> generator_roots {
		lut_singleton.exp(0),
		lut_singleton.exp(1),
		lut_singleton.exp(2),
		lut_singleton.exp(3),
		//lut_singleton.exp(4),
		//lut_singleton.exp(5),
	};
	
	poly generator({1});
	for(unsigned i = 0; i < generator_roots.size(); i++) {
		generator *= poly({ generator_roots[i], 1 });
	}

	printf("generator:\n");
	poly_print(generator);

	poly message({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1});
	//poly message({9, 8, 7, 6, 5, 4, 3, 2, 1});
	printf("message:\n");
	poly_print(message);

	poly msg_shift = message.shifted(generator.terms()-1);
	printf("msg shift:\n");
	poly_print(msg_shift);

	poly remainder = msg_shift % generator;
	printf("remainder:\n");
	poly_print(remainder);

	poly transmit = message.shifted(generator.terms()-1) + remainder;
	printf("transmit:\n");
	poly_print(transmit);

	poly syndrome_good;
	for(unsigned i = 0; i < generator_roots.size(); i++) {
		gf_2_4 s_i = transmit.evaluate(generator_roots[i]);
		syndrome_good += poly({s_i}).shifted(i);
	}
	printf("good syndrome:\n");
	poly_print(syndrome_good);
	
	poly errors_actual({0, 0, 2, 0, 0, 0, 0, 0, 0, 13});
	printf("errors_actual:\n");
	poly_print(errors_actual);

	poly received = transmit + errors_actual;
	printf("received:\n");
	poly_print(received);


	poly syndrome_bad;
	for(unsigned i = 0; i < generator_roots.size(); i++) {
		gf_2_4 s_i = received.evaluate(generator_roots[i]);
		syndrome_bad += poly({s_i}).shifted(i);
	}
	printf("bad syndrome:\n");
	poly_print(syndrome_bad);

	poly sigma;
	berlekamp(&sigma, syndrome_bad);
	printf("sigma:\n");
	poly_print(sigma);

	// calculate omega
	poly omega = sigma * syndrome_bad % poly(generator_roots.size(), 1);
	printf("omega:\n");
	poly_print(omega);

	// error loc: chien search
	poly err_eval_poly;
	std::vector<unsigned> err_locs;
	for(unsigned i = 1; i < gf_2_4::order(); i++) {
		gf_2_4 test_root = lut_singleton.exp(-(int)i);
		gf_2_4 val = sigma.evaluate(test_root);
		err_eval_poly += poly(i, val);
		if(val == 0) {
			err_locs.push_back(i % (gf_2_4::order() - 1));
		}
	}
	printf("err_eval_poly:\n");
	poly_print(err_eval_poly);
	printf("err locs:\n");
	for(unsigned i = 0; i < err_locs.size(); i++)
		printf("%u ", err_locs[i]);
	printf("\n");

	// calculate derivative of sigma
	poly sigma_deriv;
	for(unsigned i = 1; i < sigma.terms(); i += 2) {
		sigma_deriv += poly(i-1, sigma.term(i));
	}
	printf("sigma_deriv:\n");
	poly_print(sigma_deriv);

	// calculate error value polynomial
	poly correction;
	for(unsigned i = 0; i < err_locs.size(); i++) {
		int loc = (int) err_locs[i];
		gf_2_4 root = lut_singleton.exp(loc);
		gf_2_4 root_inv = lut_singleton.exp(-loc);
		// FIXME: fixup root if initial root power of gen. poly is not 0
		gf_2_4 val = root * omega.evaluate(root_inv) / sigma_deriv.evaluate(root_inv);
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
	std::vector<gf_2_4> err_vals;
	for(unsigned i = 0; i < err_locs.size(); i++) {
		gf_2_4 zi(lut_singleton.exp(err_locs[i]));

		gf_2_4 zi_inv = gf_2_4(1) / zi;
		gf_2_4 numer = omega.evaluate(zi_inv);

		gf_2_4 prod = 1;
		for(unsigned j = 0; j < err_locs.size(); j++) {
			if(j == i)
				continue;
			gf_2_4 zj(lut_singleton.exp(err_locs[j]));
			gf_2_4 t = zj / zi;
			prod *= t + gf_2_4(1);
		}
		gf_2_4 denom = zi * prod;

		gf_2_4 yi = numer / denom;
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
