#include <iostream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <assert.h>

#include "polynomial.hpp"

class doublecoeff {
public:
	doublecoeff() { val = 0.0; };
	doublecoeff(double x) { val = (double)x; }
	doublecoeff(const doublecoeff &x) { val = x.val; }
	~doublecoeff() {};
	
	doublecoeff &operator=(const doublecoeff &rhs)
	{ val = rhs.val; return *this; };

	bool operator==(const doublecoeff &rhs)
	{ return val == rhs.val; }
	bool operator!=(const doublecoeff &rhs)
	{ return val != rhs.val; }

	doublecoeff operator+(const doublecoeff &rhs) const
	{ return doublecoeff(val + rhs.val); }
	doublecoeff operator-(const doublecoeff &rhs) const
	{ return doublecoeff(val - rhs.val); }
	doublecoeff operator*(const doublecoeff &rhs) const
	{ return doublecoeff(val * rhs.val); }
	doublecoeff operator/(const doublecoeff &rhs) const
	{ return doublecoeff(val / rhs.val); }

	doublecoeff operator-() const
	{ return doublecoeff(-val); }

	doublecoeff operator+=(const doublecoeff &rhs)
	{ val += rhs.val; return *this; }
	doublecoeff operator-=(const doublecoeff &rhs)
	{ val -= rhs.val; return *this; }
	doublecoeff operator*=(const doublecoeff &rhs)
	{ val *= rhs.val; return *this; }
	doublecoeff operator/=(const doublecoeff &rhs)
	{ val /= rhs.val; return *this; }

	double val;
};

void poly_print(const polynomial<doublecoeff> &P)
{
	if(P.terms() == 0)
		printf("%g ", 0.0);
	else for(unsigned i = 0; i < P.terms(); i++)
		printf("%g ", P.term(i).val);
}

struct { polynomial<doublecoeff> A, B, Z; } add_tests[] = {
	{ { }, { }, { } },
	{ { 1.0 }, { }, { 1.0 } },
	{ { }, { 1.0 }, { 1.0 } },
	{ { 1.0, 2.0 }, { 1.0, 1.0 }, { 2.0, 3.0 } },
	{ { 1.0, 2.0, 3.0 }, { 1.0, 1.0 }, { 2.0, 3.0, 3.0 } },
};

struct { polynomial<doublecoeff> A, B, Z; } sub_tests[] = {
	{ { }, { }, { } },
	{ { 1.0 }, { }, { 1.0 } },
	{ { }, { 1.0 }, { -1.0 } },
	{ { 1.0, 2.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } },
};

struct { polynomial<doublecoeff> A, Z; } neg_tests[] = {
	{ { }, { } },
	{ { 1.0 }, { -1.0 } },
};

struct { polynomial<doublecoeff> A, B, Z; } mul_tests[] = {
	{ { 1.0, 2.0 }, { 0.0, 1.0 }, { 0.0, 1.0, 2.0 } },
	{ { 1.0, 2.0 }, { 1.0, 1.0 }, { 1.0, 3.0, 2.0 } },
};

struct { polynomial<doublecoeff> A, B, Z; } div_tests[] = {
	{ { }, { 1.0 }, { } },
	{ { 0.0, 1.0, 2.0 }, { 1.0, 2.0 }, { 0.0, 1.0 } },
	{ { 1.0, 3.0, 2.0 }, { 1.0, 2.0 }, { 1.0, 1.0 } },
	{ { -3.0, 10.0, -5.0, 3.0 }, { 1.0, 3.0 }, { 4.0, -2.0, 1.0 } },
	{ { 15.0, 0.0, -9.0, 2.0 }, { -5.0, 2.0 }, { -5.0, -2.0, 1.0 } },
	{ { 1.0, 2.0, 0.0, 3.0, 4.0 }, { 2.0, 1.0, 1.0 }, { -7.0, -1.0, 4.0 } },
};

struct { polynomial<doublecoeff> A, B, Z; } mod_tests[] = {
	{ { }, { 1.0 }, { } },
	{ { 0.0, 1.0, 2.0 }, { 1.0, 2.0 }, { } },
	{ { 1.0, 3.0, 2.0 }, { 1.0, 2.0 }, { } },
	{ { -3.0, 10.0, -5.0, 3.0 }, { 1.0, 3.0 }, { -7.0 } },
	{ { 15.0, 0.0, -9.0, 2.0 }, { -5.0, 2.0 }, { -10.0 } },
	{ { 1.0, 2.0, 0.0, 3.0, 4.0 }, { 2.0, 1.0, 1.0 }, { 15.0, 11.0 } },
};

struct { polynomial<doublecoeff> A; doublecoeff x, z; } eval_tests[] = {
	{ { 0.0, 1.0, 2.0 }, 3.0, 21.0 },
};

struct { polynomial<doublecoeff> A, B, X, Y, G; } euclid_tests[] = {
	{ { 2, 3, 1 }, { 1, 2, 1 }, { 1 }, { -1 }, { 1, 1 } },
	{ { 2, 5, 4, 1 }, { 1, 3, 3, 1 }, { 1 }, { -1 }, { 1, 2, 1 } },
};

int main(void)
{
#define DO_OP1_1(OP, OPNAME, A, Z) \
do { \
		polynomial<doublecoeff> test = OP A; \
		if(test == Z) printf("[OK]   "); \
		else printf("[FAIL] "); \
		printf(OPNAME " "); poly_print(A); printf("= "); poly_print(test); \
		printf("(=? "); poly_print(Z); printf(")\n"); \
} while(0)

#define DO_OP2_1(OP, OPNAME, A, B, Z) \
do { \
		polynomial<doublecoeff> test = A OP B; \
		if(test == Z)  printf("[OK]   "); \
		else  printf("[FAIL] "); \
		poly_print(A); printf(OPNAME " "); poly_print(B); \
		printf("= "); poly_print(test); \
		printf("(=? "); poly_print(Z); printf(")\n"); \
} while(0)

	for(unsigned i = 0; i < sizeof(add_tests) / sizeof(add_tests[0]); i++) {
		DO_OP2_1(+, "+", add_tests[i].A, add_tests[i].B, add_tests[i].Z);
	}

	for(unsigned i = 0; i < sizeof(sub_tests) / sizeof(sub_tests[0]); i++) {
		DO_OP2_1(-, "-", sub_tests[i].A, sub_tests[i].B, sub_tests[i].Z);
	}

	for(unsigned i = 0; i < sizeof(neg_tests) / sizeof(neg_tests[0]); i++) {
		DO_OP1_1(-, "-", neg_tests[i].A, neg_tests[i].Z);
	}

	for(unsigned i = 0; i < sizeof(mul_tests) / sizeof(mul_tests[0]); i++) {
		DO_OP2_1(*, "*", mul_tests[i].A, mul_tests[i].B, mul_tests[i].Z);
	}

	for(unsigned i = 0; i < sizeof(div_tests) / sizeof(div_tests[0]); i++) {
		DO_OP2_1(/, "/", div_tests[i].A, div_tests[i].B, div_tests[i].Z);
	}

	for(unsigned i = 0; i < sizeof(mod_tests) / sizeof(mod_tests[0]); i++) {
		DO_OP2_1(%, "%%", mod_tests[i].A, mod_tests[i].B, mod_tests[i].Z);
	}

	for(unsigned i = 0; i < sizeof(eval_tests) / sizeof(eval_tests[0]); i++) {
		doublecoeff test = eval_tests[i].A.evaluate(eval_tests[i].x);
		if(test == eval_tests[i].z)  printf("[OK]   ");
		else  printf("[FAIL] ");
		poly_print(eval_tests[i].A); printf("@ %f", eval_tests[i].x.val);
		printf(" = %g ", test.val);
		printf("(=? %f)\n", eval_tests[i].z.val);
	}

	for(unsigned i = 0; i < sizeof(euclid_tests) / sizeof(euclid_tests[0]); i++) {
		polynomial<doublecoeff> Gtest;
		polynomial<doublecoeff> Xtest;
		polynomial<doublecoeff> Ytest;

		polynomial<doublecoeff>::polynomial_ext_euclid(
			&Gtest, &Xtest, &Ytest, euclid_tests[i].A, euclid_tests[i].B);

		if(Gtest != euclid_tests[i].G ||
		   Xtest != euclid_tests[i].X ||
		   Ytest != euclid_tests[i].Y)
			printf("[FAIL] ");
		else
			printf("[OK]   ");

		printf("ext_euc( ");
		poly_print(euclid_tests[i].A);
		printf(", ");
		poly_print(euclid_tests[i].B);
		printf(") -> ");
		poly_print(Gtest);
		printf("(=? "); poly_print(euclid_tests[i].G); printf(") ");
		printf("X = "); poly_print(Xtest);
		printf("(=? "); poly_print(euclid_tests[i].X); printf(") ");
		printf(", Y = "); poly_print(Ytest);
		printf("(=? "); poly_print(euclid_tests[i].Y); printf(") ");
		printf("\n");
	}

	return 0;
}
