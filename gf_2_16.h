#if ! defined(GF_2_16_H_INCLUDED)
#define GF_2_16_H_INCLUDED

class gf_2_16_lut
{
public:
	gf_2_16_lut();
	~gf_2_16_lut();

	unsigned from_index(unsigned index); // returns rawval
	unsigned to_index(unsigned rawval); // returns index

	unsigned log(unsigned rawval); // returns lograwval
	unsigned exp(int lograwval); // returns rawval

private:
};

extern gf_2_16_lut lut_2_16_singleton;

class gf_2_16
{
public:
	gf_2_16() : rawval(0) { }
	gf_2_16(unsigned _rawval) : rawval(_rawval) {}
	gf_2_16(const gf_2_16 &x) : rawval(x.rawval) {}
	~gf_2_16() {}

	static gf_2_16 exp(unsigned pow)
	{ return gf_2_16(lut_2_16_singleton.exp(pow)); }

	static unsigned log(unsigned val)
	{ return lut_2_16_singleton.log(val); }

	operator unsigned()
	{ return rawval; }

	bool operator==(const gf_2_16 &rhs)
	{ return rawval == rhs.rawval; }

	bool operator==(int rhs)
	{ return rawval == (unsigned) rhs; }

	bool operator==(unsigned rhs)
	{ return rawval == rhs; }

	bool operator!=(const gf_2_16 &rhs)
	{ return rawval != rhs.rawval; }

	bool operator!=(int rhs)
	{ return rawval != (unsigned) rhs; }

	bool operator!=(unsigned rhs)
	{ return rawval != rhs; }

	gf_2_16 &operator=(const gf_2_16 &rhs)
	{ rawval = rhs.rawval; return *this; };

	gf_2_16 &operator=(unsigned rhs_rawval)
	{ rawval = rhs_rawval; return *this; };

	gf_2_16 operator+(const gf_2_16 &rhs) const
	{ gf_2_16 ret = *this; ret += rhs; return ret; }

	gf_2_16 &operator+=(const gf_2_16 &rhs)
	{ rawval ^= rhs.rawval; return *this; };

	gf_2_16 operator-(const gf_2_16 &rhs) const
	{ gf_2_16 ret = *this; ret -= rhs; return ret; }

	gf_2_16 &operator-=(const gf_2_16 &rhs)
	{ rawval ^= rhs.rawval; return *this; };

	gf_2_16 operator*(const gf_2_16 &rhs) const
	{ gf_2_16 ret = *this; ret *= rhs; return ret; };

	gf_2_16 &operator*=(const gf_2_16 &rhs)
	{
		if(rawval == 0 || rhs.rawval == 0)
			rawval = 0;
		else {
			unsigned lg_this = lut_2_16_singleton.log(rawval);
			unsigned lg_rhs = lut_2_16_singleton.log(rhs.rawval);
			rawval = lut_2_16_singleton.exp(lg_this + lg_rhs);
		}
		return *this;
	};

	gf_2_16 operator/(const gf_2_16 &rhs) const
	{ gf_2_16 ret = *this; ret /= rhs; return ret; };

	gf_2_16 &operator/=(const gf_2_16 &rhs)
	{
		if(rawval != 0) {
			unsigned lg_this = lut_2_16_singleton.log(rawval);
			unsigned lg_rhs = lut_2_16_singleton.log(rhs.rawval);
			rawval = lut_2_16_singleton.exp(lg_this - lg_rhs);
		}
		return *this;
	};

	static unsigned order()
	{ return 65536; }

private:
	unsigned rawval;
};

#endif
