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

	static void div_quot_rem(polynomial<coeff> *Q, 
	                         polynomial<coeff> *R,
	                         const polynomial<coeff> &A, 
	                         const polynomial<coeff> &B)
	{
		//  Initialize: _R(x) = A(x)
		//              align = n - j
		//
		//  at each step, calculate:
		//
		//                Q_m x^align +...
		//              -------------------
		//  B_j(x) +... ) R_n x +...
		//              -   T(x)
		//              ----------
		//                       R_n-1 x +...
		//
		//  do:
		//    Q_m = -R_n / B_j
		//    T(x) = (B(x) * x^align) * Q_m
		//    R(x) = R(x) + T(x)
		//    align = align - 1
		//  while align >= 0

		polynomial _R(A.p);
		std::vector<coeff> _Q;
		int align = A.terms() - B.terms();
		
		if(align < 0) {
			// special case if order(A) < order(B)
			*Q = polynomial(_Q);
			*R = _R;
			return;
		}

		coeff Q_m;
		polynomial L, T;
		do {
			coeff _R_lead = _R.terms() ? _R[_R.terms()-1] : coeff(0);
			coeff B_lead = B.terms() ? B[B.terms()-1] : coeff(0);
			Q_m = _R_lead / B_lead;
			T = (B << align) * Q_m;
			if(_R.terms() < T.terms())
				_Q.push_back(0);
			else {
				_Q.push_back(Q_m);
				_R -= T;
			}
		} while(--align >= 0);

		std::reverse(_Q.begin(), _Q.end());
		*Q = polynomial(_Q);
		*R = _R;
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
