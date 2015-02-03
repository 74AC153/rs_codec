#if ! defined(RS_CODEC_H_INCLUDED)
#define RS_CODEC_H_INCLUDED

#include "polynomial.hpp"

template <class coeff_type>
void rs_encode(
	polynomial<coeff_type> *transmit,
	const polynomial<coeff_type> &message,
	const std::vector<coeff_type> &generator_roots)
{
	polynomial<coeff_type> generator = polynomial<coeff_type>::from_roots(generator_roots);
	polynomial<coeff_type> msg_shift = message << (generator.terms()-1);
	polynomial<coeff_type> remainder = msg_shift % generator;
	// transmit = msg * x^(2t) % generator
	*transmit = msg_shift + remainder;
}

template <class coeff_type>
void rs_calc_syndrome(
	polynomial<coeff_type> *syndrome,
	const polynomial<coeff_type> &msg,
	const std::vector<coeff_type> &generator_roots)
{
	// syndrome[i] = msg.eval(gen_root[i]) * x^i
	*syndrome = polynomial<coeff_type>();
	for(unsigned i = 0; i < generator_roots.size(); i++) {
		coeff_type s_i = msg.evaluate(generator_roots[i]);
		*syndrome += polynomial<coeff_type>(i, s_i);
	}
}

template <class coeff_type>
void rs_berlekamp(
	polynomial<coeff_type> *sigma_out,
	const polynomial<coeff_type> &syndrome)
{
	unsigned K = 1;
	unsigned L = 0;
	polynomial<coeff_type> sigma { 1 };
	polynomial<coeff_type> sigma_new;
	polynomial<coeff_type> C { 0, 1 };

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

		C *= polynomial<coeff_type>(1, 1);
		K++;
	}

	*sigma_out = sigma;
}

template <class coeff_type>
void rs_chien_search(
	std::vector<unsigned> *err_locs,
	unsigned len,
	const polynomial<coeff_type> &sigma)
{
	//poly err_eval_poly;
	for(unsigned i = 0; i < len; i++) {
		coeff_type test_root = coeff_type::exp(-(int)i);
		coeff_type val = sigma.evaluate(test_root);
		//err_eval_poly += poly(i, val);
		if(val == 0) {
			err_locs->push_back(i);
		}
	}
}

template <class coeff_type>
void rs_forney(
	polynomial<coeff_type> *correction,
	const polynomial<coeff_type> &sigma,
	const polynomial<coeff_type> &syndrome,
	const std::vector<unsigned> &err_locs,
	const std::vector<coeff_type> &generator_roots)
{
	// calculate derivative of sigma
	// "set even-powered coefficients to 0, then divide by x"
	polynomial<coeff_type> sigma_deriv;
	for(unsigned i = 1; i < sigma.terms(); i += 2) {
		sigma_deriv += polynomial<coeff_type>(i-1, sigma[i]);
	}

	// calculate omega
	// omega = (syndrome * sigma) mod x^(2t)
	polynomial<coeff_type> omega =
		(sigma * syndrome) % polynomial<coeff_type>(generator_roots.size(), 1);

	// forney algorithm
	// err_j = X_j^(1-b) omega(X_j^(-1)) / lambda_deriv(X_j^(-1))
	// X_j = a^(err_loc[j])
	// b = exponent for first generator root
	// NB: assume b is 0 for now (see FIXME below)
	*correction = polynomial<coeff_type>();
	for(unsigned i = 0; i < err_locs.size(); i++) {
		int loc = (int) err_locs[i];
		coeff_type root = coeff_type::exp(loc);
		coeff_type root_inv = coeff_type::exp(-loc);
		// FIXME: fixup root if initial power of gen. poly is not 0
		coeff_type val = root * omega.evaluate(root_inv) / sigma_deriv.evaluate(root_inv);
		*correction += polynomial<coeff_type>(loc, val);
	}
}

#endif
