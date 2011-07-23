#ifndef _VECTOR4_H
#define _VECTOR4_H

#include <math.h>
#include <stdio.h>

// Need this pragma due to operator[] implementation.
#pragma pack(4)

template <typename T>
class vector4 {
public:
	T x,y,z,w;

	// Constructor definitions are outside class declaration to enforce that
	// only float and double versions are possible.
	vector4();
	vector4(const vector4<float > &v);
	vector4(const vector4<double> &v);
	vector4(T val);
	vector4(T _x, T _y, T _z, T _w);
	vector4(const float  vals[4]);
	vector4(const double vals[4]);

	const T& operator[](const size_t i) const { return (const_cast<const T *>(&x))[i]; }
	T& operator[](const size_t i) { return (&x)[i]; }

	vector4 operator+(const vector4 &a) const { return vector4 (a.x+x, a.y+y, a.z+z); }
	vector4 &operator+=(const vector4 &a) { x+=a.x; y+=a.y; z+=a.z; return *this; }
	vector4 &operator-=(const vector4 &a) { x-=a.x; y-=a.y; z-=a.z; return *this; }
	vector4 &operator*=(const float a) { x*=a; y*=a; z*=a; return *this; }
	vector4 &operator*=(const double a) { x*=a; y*=a; z*=a; return *this; }
	vector4 &operator/=(const float a) { const T inva = T(1.0/a); x*=inva; y*=inva; z*=inva; w*=inva; return *this; }
	vector4 &operator/=(const double a) { const T inva = T(1.0/a); x*=inva; y*=inva; z*=inva; w*=inva; return *this; }
	vector4 operator-(const vector4 &a) const { return vector4(x-a.x, y-a.y, z-a.z); }
	vector4 operator-() const { return vector4(-x, -y, -z); }
	bool operator==(const vector4 &a) const { return ((a.x==x)&&(a.y==y)&&(a.z==z)&&(a.w==w)); }
	bool operator!=(const vector4 &a) const { return ((a.x!=x)||(a.y!=y)||(a.z!=z)||(a.w!=w)); }

	friend vector4 operator*(const vector4 &a, const float  scalar) { return vector4(T(a.x*scalar), T(a.y*scalar), T(a.z*scalar), T(a.w*scalar)); }
	friend vector4 operator*(const vector4 &a, const double scalar) { return vector4(T(a.x*scalar), T(a.y*scalar), T(a.z*scalar), T(a.w*scalar)); }
	friend vector4 operator*(const float  scalar, const vector4 &a) { return a*scalar; }
	friend vector4 operator*(const double scalar, const vector4 &a) { return a*scalar; }
	friend vector4 operator/(const vector4 &a, const float  scalar) { const T inv = 1.0/scalar; return vector4(a.x*inv, a.y*inv, a.z*inv, a.w*inv); }
	friend vector4 operator/(const vector4 &a, const double scalar) { const T inv = 1.0/scalar; return vector4(a.x*inv, a.y*inv, a.z*inv, a.w*inv); }

	//vector4 Cross(const vector4 &b) const { return vector4 (y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
	T Dot(const vector4 &b) const { return x*b.x + y*b.y + z*b.z + w*b.w; }
	T Length() const { return sqrt (x*x + y*y + z*z + w*w); }
	T LengthSqr() const { return x*x + y*y + z*z + w*w; }
	vector4 Normalized() const { const T l = 1.0f / sqrt(x*x + y*y + z*z + w*w); return vector4(x*l, y*l, z*l, w*l); }
	vector4 NormalizedSafe() const {
		T l = sqrt(x*x + y*y + z*z + w*w);
		if (l==0.0) return vector4(1,0,0,1);
		return vector4(x/l, y/l, z/l, w/l);
	}

	void Print() const { printf("v(%f,%f,%f,%f)\n", x, y, z,w); }

	/* Rotate this vector about point o, in axis defined by v. */
	void ArbRotateAroundPoint(const vector4 &o, const vector4 &__v, T ang) {
		vector4 t;
		T a = o.x;
		T b = o.y;
		T c = o.z;
		T u = __v.x;
		T v = __v.y;
		T w = __v.z;
		T cos_a = cos (ang);
		T sin_a = sin (ang);
		T inv_poo = 1.0f/(u*u+v*v+w*w);
		t.x = a*(v*v+w*w)+u*(-b*v-c*w+u*x+v*y+w*z)+(-a*(v*v+w*w)+u*(b*v+c*w-v*y-w*z)+(v*v+w*w)*x)*cos_a+
			sqrtf (u*u+v*v+w*w)*(-c*v+b*w-w*y+v*z)*sin_a;
		t.x *= inv_poo;
		t.y = b*(u*u+w*w)+v*(-a*u-c*w+u*x+v*y+w*z)+(-b*(u*u+w*w)+v*(a*u+c*w-u*x-w*z)+(u*u+w*w)*y)*cos_a+
			sqrtf (u*u+v*v+w*w)*(-c*u-a*w+w*x-u*z)*sin_a;
		t.y *= inv_poo;
		t.z = c*(u*u+v*v)+w*(-a*u+b*v+u*x+v*y+w*z)+(-c*(u*u+v*v)+w*(a*u+b*v-u*x-v*y)+(u*u+v*v)*z)*cos_a+
			sqrtf (u*u+v*v+w*w)*(-b*u+a*v-v*x+u*y)*sin_a;
		t.z *= inv_poo;
		*this = t;
	}

	/* Rotate this vector about origin, in axis defined by v. */
	void ArbRotate(const vector4 &__v, T ang) {
		vector4 t;
		T u = __v.x;
		T v = __v.y;
		T w = __v.z;
		T cos_a = cos(ang);
		T sin_a = sin(ang);
		T inv_poo = 1.0f/(u*u+v*v+w*w);
		t.x = u*(u*x+v*y+w*z)+(u*(-v*y-w*z)+(v*v+w*w)*x)*cos_a+
			sqrtf (u*u+v*v+w*w)*(-w*y+v*z)*sin_a;
		t.x *= inv_poo;
		t.y = v*(u*x+v*y+w*z)+(v*(-u*x-w*z)+(u*u+w*w)*y)*cos_a+
			sqrtf (u*u+v*v+w*w)*(w*x-u*z)*sin_a;
		t.y *= inv_poo;
		t.z = w*(u*x+v*y+w*z)+(w*(-u*x-v*y)+(u*u+v*v)*z)*cos_a+
			sqrtf (u*u+v*v+w*w)*(-v*x+u*y)*sin_a;
		t.z *= inv_poo;
		*this = t;
	}
};

// These are here in this manner to enforce that only float and double versions are possible.
template<> inline vector4<float >::vector4() {}
template<> inline vector4<double>::vector4() {}
template<> inline vector4<float >::vector4(const vector4<float > &v): x(v.x), y(v.y), z(v.z) {}
template<> inline vector4<float >::vector4(const vector4<double> &v): x(float(v.x)), y(float(v.y)), z(float(v.z)) {}
template<> inline vector4<double>::vector4(const vector4<float > &v): x(v.x), y(v.y), z(v.z) {}
template<> inline vector4<double>::vector4(const vector4<double> &v): x(v.x), y(v.y), z(v.z) {}
template<> inline vector4<float >::vector4(float  val): x(val), y(val), z(val) {}
template<> inline vector4<double>::vector4(double val): x(val), y(val), z(val) {}
template<> inline vector4<float >::vector4(float  _x, float  _y, float  _z, float  _w): x(_x), y(_y), z(_z), w(_w) {}
template<> inline vector4<double>::vector4(double _x, double _y, double _z, double _w): x(_x), y(_y), z(_z), w(_w) {}
template<> inline vector4<float >::vector4(const float  vals[4]): x(vals[0]), y(vals[1]), z(vals[2]) {}
template<> inline vector4<float >::vector4(const double vals[4]): x(float(vals[0])), y(float(vals[1])), z(float(vals[2])) {}
template<> inline vector4<double>::vector4(const float  vals[4]): x(vals[0]), y(vals[1]), z(vals[2]) {}
template<> inline vector4<double>::vector4(const double vals[4]): x(vals[0]), y(vals[1]), z(vals[2]) {}

#pragma pack()

typedef vector4<float > vector4f;
typedef vector4<double> vector4d;

#endif /* _VECTOR4_H */
