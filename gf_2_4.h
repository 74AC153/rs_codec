#if ! defined(GF_2_4_H_INCLUDED)
#define GF_2_4_H_INCLUDED

class gf_2_4_lut
{
public:
	gf_2_4_lut();
	~gf_2_4_lut();

	unsigned from_index(unsigned index);
	unsigned to_index(unsigned rawval);

	unsigned log(unsigned rawval);
	unsigned exp(int lograwval);

private:
};

extern gf_2_4_lut lut_singleton;

class gf_2_4
{
public:
	gf_2_4() : rawval(0) { }
	gf_2_4(unsigned _rawval) : rawval(_rawval) {}
	gf_2_4(const gf_2_4 &x) : rawval(x.rawval) {}
	~gf_2_4() {}

	operator unsigned()
	{ return rawval; }

	bool operator==(const gf_2_4 &rhs)
	{ return rawval == rhs.rawval; }

	bool operator!=(const gf_2_4 &rhs)
	{ return rawval != rhs.rawval; }

	gf_2_4 &operator=(const gf_2_4 &rhs)
	{ rawval = rhs.rawval; return *this; };

	gf_2_4 &operator=(unsigned rhs_rawval)
	{ rawval = rhs_rawval; return *this; };

	gf_2_4 operator+(const gf_2_4 &rhs) const
	{ gf_2_4 ret = *this; ret += rhs; return ret; }

	gf_2_4 &operator+=(const gf_2_4 &rhs)
	{ rawval ^= rhs.rawval; return *this; };

	gf_2_4 operator-(const gf_2_4 &rhs) const
	{ gf_2_4 ret = *this; ret -= rhs; return ret; }

	gf_2_4 &operator-=(const gf_2_4 &rhs)
	{ rawval ^= rhs.rawval; return *this; };

	gf_2_4 operator*(const gf_2_4 &rhs) const
	{ gf_2_4 ret = *this; ret *= rhs; return ret; };

	gf_2_4 &operator*=(const gf_2_4 &rhs)
	{
		if(rawval == 0 || rhs.rawval == 0)
			rawval = 0;
		else {
			unsigned lg_this = lut_singleton.log(rawval);
			unsigned lg_rhs = lut_singleton.log(rhs.rawval);
			rawval = lut_singleton.exp(lg_this + lg_rhs);
		}
		return *this;
	};

	gf_2_4 operator/(const gf_2_4 &rhs) const
	{ gf_2_4 ret = *this; ret /= rhs; return ret; };

	gf_2_4 &operator/=(const gf_2_4 &rhs)
	{
		if(rawval != 0) {
			unsigned lg_this = lut_singleton.log(rawval);
			unsigned lg_rhs = lut_singleton.log(rhs.rawval);
			rawval = lut_singleton.exp(lg_this - lg_rhs);
		}
		return *this;
	};

private:
	unsigned rawval;
};

#endif