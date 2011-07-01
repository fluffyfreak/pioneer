#ifndef _VECTOR2_H
#define _VECTOR2_H

#include <math.h>
#include <stdio.h>

// Need this pragma due to operator[] implementation.
#pragma pack(4)

template <typename T>
class vector2 {
public:
	T x,y;

	// Constructor definitions are outside class declaration to enforce that
	// only float and double versions are possible.
	vector2();
	vector2(const vector2<float > &v);
	vector2(const vector2<double> &v);
	vector2(T val);
	vector2(T _x, T _y);
	vector2(const float  vals[2]);
	vector2(const double vals[2]);

	const T& operator[](const size_t i) const { return (const_cast<const T *>(&x))[i]; }
	T& operator[](const size_t i) { return (&x)[i]; }

	vector2 operator+(const vector2 &a) const { return vector2 (a.x+x, a.y+y); }
	vector2 &operator+=(const vector2 &a) { x+=a.x; y+=a.y; return *this; }
	vector2 &operator-=(const vector2 &a) { x-=a.x; y-=a.y; return *this; }
	vector2 &operator*=(const float a) { x*=a; y*=a; return *this; }
	vector2 &operator*=(const double a) { x*=a; y*=a; return *this; }
	vector2 &operator/=(const float a) { const T inva = T(1.0/a); x*=inva; y*=inva; return *this; }
	vector2 &operator/=(const double a) { const T inva = T(1.0/a); x*=inva; y*=inva; return *this; }
	vector2 operator-(const vector2 &a) const { return vector2(x-a.x, y-a.y); }
	vector2 operator-() const { return vector2(-x, -y); }
	bool operator==(const vector2 &a) const { return ((a.x==x)&&(a.y==y)); }
	bool operator!=(const vector2 &a) const { return ((a.x!=x)||(a.y!=y)); }

	friend vector2 operator*(const vector2 &a, const float  scalar) { return vector2(T(a.x*scalar), T(a.y*scalar)); }
	friend vector2 operator*(const vector2 &a, const double scalar) { return vector2(T(a.x*scalar), T(a.y*scalar)); }
	friend vector2 operator*(const float  scalar, const vector2 &a) { return a*scalar; }
	friend vector2 operator*(const double scalar, const vector2 &a) { return a*scalar; }
	friend vector2 operator/(const vector2 &a, const float  scalar) { const T inv = 1.0/scalar; return vector2(a.x*inv, a.y*inv); }
	friend vector2 operator/(const vector2 &a, const double scalar) { const T inv = 1.0/scalar; return vector2(a.x*inv, a.y*inv); }

	T Dot(const vector2 &b) const { return x*b.x + y*b.y; }
	T Length() const { return sqrt (x*x + y*y); }
	T LengthSqr() const { return x*x + y*y; }
	vector2 Normalized() const { const T l = 1.0f / sqrt(x*x + y*y); return vector2(x*l, y*l); }
	vector2 NormalizedSafe() const {
		T l = sqrt(x*x + y*y);
		if (l==0.0) return vector2(1,0);
		return vector2(x/l, y/l);
	}

	void Print() const { printf("v(%f,%f)\n", x, y); }
};

// These are here in this manner to enforce that only float and double versions are possible.
template<> inline vector2<float >::vector2() {}
template<> inline vector2<double>::vector2() {}
template<> inline vector2<float >::vector2(const vector2<float > &v): x(v.x), y(v.y) {}
template<> inline vector2<float >::vector2(const vector2<double> &v): x(float(v.x)), y(float(v.y)) {}
template<> inline vector2<double>::vector2(const vector2<float > &v): x(v.x), y(v.y) {}
template<> inline vector2<double>::vector2(const vector2<double> &v): x(v.x), y(v.y) {}
template<> inline vector2<float >::vector2(float  val): x(val), y(val) {}
template<> inline vector2<double>::vector2(double val): x(val), y(val) {}
template<> inline vector2<float >::vector2(float  _x, float  _y): x(_x), y(_y) {}
template<> inline vector2<double>::vector2(double _x, double _y): x(_x), y(_y) {}
template<> inline vector2<float >::vector2(const float  vals[2]): x(vals[0]), y(vals[1]) {}
template<> inline vector2<float >::vector2(const double vals[2]): x(float(vals[0])), y(float(vals[1])) {}
template<> inline vector2<double>::vector2(const float  vals[2]): x(vals[0]), y(vals[1]) {}
template<> inline vector2<double>::vector2(const double vals[2]): x(vals[0]), y(vals[1]) {}

#pragma pack()

typedef vector2<float > vector2f;
typedef vector2<double> vector2d;

#endif /* _VECTOR2_H */
