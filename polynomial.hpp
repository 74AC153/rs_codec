#if ! defined(POLYNOMIAL_HPP_INCLUDED)
#define POLYNOMIAL_HPP_INCLUDED

#include <algorithm>
#include <vector>

template <class coeff>
class polynomial {
public:
	polynomial() {}
	polynomial(std::initializer_list<coeff> c) { p = std::vector<coeff>(c); }
	polynomial(const std::vector<coeff> &x) { p = x; }
	polynomial(const polynomial &x) { p = x.p; }
	polynomial(size_t order, coeff val) { p.resize(order+1); p[order] = val; }
	~polynomial() {}

	static polynomial from_roots(const std::vector<coeff> &roots)
	{
		polynomial ret(0, 1);
		for(unsigned i = 0; i < roots.size(); i++)
			ret *= polynomial({ roots[i], 1 });
		return ret;
	}

	std::vector<coeff> rawdata()
	{ return p; }

	polynomial &operator=(const polynomial &rhs)
	{ p = rhs.p; return *this; }

	bool operator==(const polynomial &rhs)
	{
		if(p.size() != rhs.p.size())
			return false;
		for(size_t i = 0; i < p.size() ; i++)
			if(p[i] != rhs.p[i])
				return false;
		return true;
	}

	bool operator!=(const polynomial &rhs)
	{ return ! (*this == rhs); }

	operator bool()
	{ return p.size() != 0; }

	polynomial operator+(const polynomial &rhs) const
	{ polynomial ret(*this); ret += rhs; return ret; }

	polynomial &operator+=(const polynomial &rhs)
	{
		if(p.size() < rhs.p.size())
			p.resize(rhs.p.size());
		for(size_t i = 0; i < rhs.p.size();  i++)
			p[i] += rhs.p[i];
		trim(&p);
		return *this;
	}

	polynomial operator-(const polynomial &rhs) const
	{ polynomial ret(*this); ret -= rhs; return ret; }

	polynomial &operator-=(const polynomial &rhs)
	{
		if(p.size() < rhs.p.size())
			p.resize(rhs.p.size());
		for(size_t i = 0; i < rhs.p.size();  i++)
			p[i] -= rhs.p[i];
		trim(&p);
		return *this;
	}

	polynomial operator-() const
	{
		std::vector<coeff> ret(p.size());
		for(size_t i = 0; i < p.size(); i++)
			ret[i] = -p[i];
		trim(&ret);
		return polynomial(ret);
	}

	polynomial operator*(const polynomial &rhs) const
	{ polynomial ret(*this); ret *= rhs; return ret; }

	polynomial &operator*=(const polynomial &rhs)
	{
		size_t i, j;
	
		std::vector<coeff> ret(p.size() + rhs.p.size() - 1);
	
		for(i = 0; i < p.size(); i++)
			for(j = 0; j < rhs.p.size(); j++)
				ret.at(i+j) += p.at(i) * rhs.p.at(j);

		trim(&ret);
		p = ret;
		return *this;
	}

	polynomial operator/(const polynomial &rhs) const
	{ polynomial ret(*this); ret /= rhs; return ret; }

	polynomial &operator/=(const polynomial &rhs)
	{
		polynomial R;
		div_quot_rem(this, &R, *this, rhs);
		return *this;
	}

	polynomial operator%(const polynomial &rhs) const
	{ polynomial ret(*this); ret %= rhs; return ret; }

	polynomial &operator%=(const polynomial &rhs)
	{
		polynomial Q;
		div_quot_rem(&Q, this, *this, rhs);
		return *this;
	}

	polynomial operator*(const coeff &rhs) const
	{ polynomial ret(*this); ret *= rhs; return ret; }

	polynomial &operator*=(const coeff &rhs)
	{ *this *= polynomial(0, rhs); return *this; }

	polynomial operator<<(size_t rhs) const
	{ polynomial ret(*this); ret <<= rhs; return ret; }

	polynomial &operator<<=(size_t rhs)
	{ *this *= polynomial(rhs, 1); return *this; }

	polynomial operator>>(size_t rhs) const
	{ polynomial ret(*this); ret >>= rhs; return ret; }

	polynomial &operator>>=(size_t rhs)
	{ *this /= polynomial(rhs, 1); return *this; }

	polynomial operator/(const coeff &rhs) const
	{ polynomial ret(*this); ret /= polynomial(0, rhs); return ret; }

	polynomial &operator/=(const coeff &rhs)
	{ *this /= polynomial(rhs); return *this; }

	coeff evaluate(const coeff &x) const
	{
		coeff ret = 0;
		for(unsigned i = p.size(); i; i--) {
			ret = ret * x + p[i-1];
		}
		return ret;
	}

	size_t terms() const
	{ return p.size(); }

	const coeff &operator[](size_t rhs) const
	{ return p[rhs]; }

	coeff &operator[](size_t rhs)
	{ return p[rhs]; }

	void plus_s_Bx_x_pow_n(const coeff &s,
	                       const polynomial<coeff> &B,
	                       unsigned n)
	{
		unsigned B_shifted_len = B.terms() + n;
		unsigned len = p.size() > B_shifted_len ? p.size() : B_shifted_len;
		p.resize(len);
		for(unsigned i = len; i > n; i--) {
			if(i <= B_shifted_len)
				p[i-1] += s * B[i - 1 - n];
		}
		trim(&p);
	}

	static void div_quot_rem(polynomial<coeff> *Q, 
	                         polynomial<coeff> *R,
	                         const polynomial<coeff> &A, 
	                         const polynomial<coeff> &B)
	{
		//  Initialize: _R(x) = A(x)
		//              align = n - j
		//
		//  at each step, calculate (until n < j):
		//
		//                Q_m x^(n-j) +...
		//              -------------------
		//  B_j(x) +... ) R_n x +...
		//              -   T(x)
		//              ----------
		//                       R_n-1 x +...
		//
		// Q_m = R_n / B_j
		// T(x) = B(x) * Q_m * x^(n-j)
		// R(x) = R(x) - T(x)

		*R = A;
		int B_shift = R->terms() - B.terms();
		Q->p.resize(0);
		if(B_shift < 0)
			return;
		Q->p.resize(B_shift+1);
		while(1) {
			B_shift = R->terms() - B.terms();
			if(B_shift < 0)
				break;
			coeff R_lead = R->terms() ? (*R)[R->terms()-1] : coeff(0);
			coeff B_lead = B.terms() ? B[B.terms()-1] : coeff(0);
			coeff Q_m(R_lead / B_lead);
			R->plus_s_Bx_x_pow_n(-Q_m, B, B_shift);
			(*Q)[B_shift] = Q_m;
		}
		trim(&Q->p);
	}

private:
	std::vector<coeff> p;
	static void trim(std::vector<coeff> *v)
	{
		size_t i;
		for(i = v->size(); i && (v->at(i-1) == 0); i--);
		v->resize(i);
	}
};

#endif
