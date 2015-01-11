template <class coeff>
class polynomial {
public:
	polynomial() {}
	polynomial(std::initializer_list<coeff> c) { p = std::vector<coeff>(c); }
	polynomial(const std::vector<coeff> &x) { p = x; }
	polynomial(const polynomial &x) { p = x.p; }
	polynomial(size_t order, coeff val) { p.resize(order+1); p[order] = val; }
	~polynomial() {}

	polynomial &operator=(const polynomial &rhs)
	{
		p = rhs.p;
		return *this;
	}

	polynomial &operator=(std::initializer_list<coeff> c)
	{
		p = std::vector<coeff>(c);
		return *this;
	}

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
	{
		return ! (*this == rhs);
	}

	operator bool()
	{
		return p.size() != 0;
	}

	polynomial operator+(const polynomial &rhs) const
	{
		polynomial ret(*this);
		ret += rhs;
		return ret;
	}

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
	{
		polynomial ret(*this);
		ret -= rhs;
		return ret;
	}

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
	{
		polynomial ret(*this);
		ret *= rhs;
		return ret;
	}

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
	{
		polynomial Q, R;
		div_quot_rem(&Q, &R, *this, rhs);
		return Q;
	}

	polynomial &operator/=(const polynomial &rhs)
	{
		polynomial R;
		div_quot_rem(this, &R, *this, rhs);
	}

	polynomial operator%(const polynomial &rhs) const
	{
		polynomial Q, R;
		div_quot_rem(&Q, &R, *this, rhs);
		return R;
	}

	polynomial &operator%=(const polynomial &rhs)
	{
		polynomial Q;
		div_quot_rem(&Q, this, *this, rhs);
	}

	coeff evaluate(const coeff &x)
	{
		coeff ret = 0;
		for(unsigned i = p.size(); i; i--) {
			ret = ret * x + p[i-1];
		}
		return ret;
	}

	size_t terms() const
	{ return p.size(); }

	coeff lead_term() const
	{ return p.size() ? p[p.size()-1] : 0; }

	const coeff &term(size_t degree) const
	{ return p[degree]; }

	coeff &term(size_t degree)
	{ return p[degree]; }

	static void div_quot_rem(polynomial<coeff> *Q, 
	                         polynomial<coeff> *R,
	                         const polynomial<coeff> &A, 
	                         const polynomial<coeff> &B)
	{

		//  Initialize: R(x) = A(x)
		//             align = n - j
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
		//    L(x) = x^align
		//    Q_m = -R_n / B_j
		//    T(x) = (B(x) * L(X)) * Q_m
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
			L = polynomial(align, 1);
			Q_m = _R.lead_term() / B.lead_term();
			T = B * L * polynomial({Q_m});
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

	// given A(x), B(x), A > B
	// return G, X, Y, s.t.:
	// A(x) * X(x) + B(x) * Y(x) = G(x)
	static void polynomial_ext_euclid(
		polynomial<coeff> *G,
		polynomial<coeff> *X,
		polynomial<coeff> *Y,
		const polynomial<coeff> &A,
		const polynomial<coeff> &B)
	{
		polynomial<coeff> S2({1}), T2({0}), R2(A);
		polynomial<coeff> S1({0}), T1({1}), R1(B);
		polynomial<coeff> S0, T0, R0, Q;

		polynomial<coeff> checkS, checkT, checkR;

		while(1) {
			Q = R2 / R1;
			S0 = S2 - Q * S1;
			T0 = T2 - Q * T1;
			R0 = R2 - Q * R1;

			checkS = S0 * A;
			checkT = T0 * B;
			checkR = checkS + checkT;

			if(! R0)
				break;

			R2 = R1; R1 = R0;
			S2 = S1; S1 = S0;
			T2 = T1; T1 = T0;
		}

		*G = R1;
		*X = S1;
		*Y = T1;
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
