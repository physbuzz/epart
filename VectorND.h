#ifndef VECTORND_H
#define VECTORND_H
#include <iostream>
#include <math.h>
#include <array>

/*
template<typename Float,int DIM>
class VectorND {
public:
    std::array<Float,DIM> x;
    VectorND();
    VectorND(std::array<Float,DIM> y);

    //vector addition and subtraction
    VectorND operator+(const VectorND& b) const;
    VectorND& operator+= (const VectorND& b);
    VectorND operator- (const VectorND& b) const;
    VectorND& operator-= (const VectorND& b);
    VectorND operator* (Float b) const;
    VectorND& operator*= (Float b);
    VectorND operator/ (Float b) const;
    VectorND& operator/= (Float b);
    //return the length squared of the vector
    Float length2() const;
    //return the length of the vector
    Float length() const;
    //calculate length and mutate to (x/length,y/length).
    void normalize();
    //return normalized version of the vector.
    VectorND normalized() const;
    //return the product of all elements of the vector
    Float product() const;
    //clamp each component of the vector to between mmin and mmax, inclusive.
    void clampCube(VectorND<Float,DIM> mmin,VectorND<Float,DIM> mmax);
};

template<typename Float,int DIM>
inline std::ostream& operator<<(std::ostream& out,const VectorND<Float,DIM>& c);
template<typename Float,int DIM>
inline VectorND<Float,DIM> operator* (Float b,const VectorND<Float,DIM>& c);
//Linearly interpolate between vec1 and vec2.
template<typename Float,int DIM>
inline VectorND<Float,DIM> lerp(Float t,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2);
//cubic interpolation between v1 and v2 with control points v0 v3.
template<typename Float,int DIM>
inline VectorND<Float,DIM> cerp(Float t,const VectorND<Float,DIM>& vec0,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2,const VectorND<Float,DIM>& vec3);
*/

template<typename Float,int DIM>
class VectorND {
    static_assert(0<DIM,"VectorND DIM must be positive");
public:
    std::array<Float,DIM> x;
    VectorND() : x{} {}
    VectorND(Float arg) : x{} {
        for(int i=0;i<DIM;i++){
            x[i]=arg;
        }
    }
    VectorND(std::array<Float,DIM> y) : x(y) {}

    Float& operator[](size_t arg){
        return x[arg];
    }
    const Float& operator[](size_t arg) const{
        return x[arg];
    }

    bool operator==(const VectorND& b) const {
        bool ret=true;
        for(int i=0;i<DIM;i++){
            ret=ret&&(x[i]==b[i]);
        }
        return ret;
    }

    //vector addition and subtraction
    VectorND operator+(const VectorND& b) const {
        VectorND c;
        for(int i=0;i<DIM;i++)
            c.x[i]=x[i]+b.x[i];
        return c;
    }
    inline VectorND& operator+= (const VectorND& b){
        for(int i=0;i<DIM;i++)
            x[i]+=b.x[i];
        return *this;
    }
    VectorND operator- (const VectorND& b) const{
        VectorND c;
        for(int i=0;i<DIM;i++)
            c.x[i]=x[i]-b.x[i];
        return c;
    }
    VectorND& operator-= (const VectorND& b){
        for(int i=0;i<DIM;i++)
            x[i]-=b.x[i];
        return *this;
    }

    //multiplication/division on righthand side by scalar.
    VectorND operator* (Float b) const {
        VectorND c;
        for(int i=0;i<DIM;i++)
            c.x[i]=x[i]*b;
        return c;
    }
    VectorND& operator*= (Float b) {
        for(int i=0;i<DIM;i++)
            x[i]*=b;
        return *this;
    }
    VectorND operator/ (Float b) const{
        VectorND c;
        for(int i=0;i<DIM;i++)
            c.x[i]=x[i]/b;
        return c;
    }

    VectorND& operator/= (Float b){
        for(int i=0;i<DIM;i++)
            x[i]/=b;
        return *this;
    }
    //return the length squared of the vector
    Float length2() const {
        Float l2=Float(0);
        for(int i=0;i<DIM;i++){
            l2+=x[i]*x[i];
        }
        return l2;
    }

    //return the length squared of the vector
    Float dot(const VectorND& arg) const {
        Float ret=Float(0);
        for(int i=0;i<DIM;i++){
            ret+=x[i]*arg.x[i];
        }
        return ret;
    }

    //return the length of the vector
    Float length() const {
        return std::sqrt(length2());
    }

    //calculate length and mutate to (x/length,y/length).
    void normalize() {
        *this/=length();
    }

    Float max() const {
        Float m=x[0];
        for(int i=1;i<DIM;i++)
            m=(x[i]>m)?x[i]:m;
        return m;
    }
    Float min() const {
        Float m=x[0];
        for(int i=1;i<DIM;i++)
            m=(x[i]<m)?x[i]:m;
        return m;
    }

    //return normalized version of the vector.
    VectorND normalized() const{
        return VectorND(*this)/length();
    }

    Float product() const {
        Float prod=Float(1);
        for(int i=0;i<DIM;i++){
            prod*=x[i];
        }
        return prod;
    }

    void clampCube(VectorND<Float,DIM> mmin,VectorND<Float,DIM> mmax) {
        for(int i=0;i<DIM;i++){
            if(x[i]<mmin[i])
                x[i]=mmin[i];
            if(x[i]>mmax[i])
                x[i]=mmax[i];
        }
    }
};

template<typename Float,int DIM>
inline std::ostream& operator<<(std::ostream& out,const VectorND<Float,DIM>& c) {
    out<<'(';
    for(int i=0;i<DIM-1;i++)
        out<<c.x[i]<<',';
    out<<c.x[DIM-1]<<")";
    return out;
}
//Vector left multiplication
template<typename Float,int DIM>
inline VectorND<Float,DIM> operator* (Float b,const VectorND<Float,DIM>& c){
    VectorND<Float,DIM> d;
    for(int i=0;i<DIM;i++){
        d.x[i]=b*c.x[i];
    }
    return d;
}

//Linearly interpolate between vec1 and vec2.
template<typename Float,int DIM>
inline VectorND<Float,DIM> lerp(Float t,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2) {
    return (vec2-vec1)*t+vec1;
}

//cubic interpolation between v1 and v2 with control points v0 v3.
//This particular form has the cool property that cerp'(1,v0,v1,v2,v3) == cerp'(0,v1,v2,v3,v4), where
//prime denotes differentiation with respect to the first parameter t. So cubic interpolation is great 
//for making long continuous paths through a bunch of points.
template<typename Float,int DIM>
inline VectorND<Float,DIM> cerp(Float t,const VectorND<Float,DIM>& vec0,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2,const VectorND<Float,DIM>& vec3) {
    Float t2=t*t;
    VectorND<Float,DIM> a0=vec3-vec2-vec0+vec1;
    VectorND<Float,DIM> a1=vec0-vec1-a0;
    VectorND<Float,DIM> a2=vec2-vec0;
    return a0*t*t2+a1*t2+a2*t+vec1;
}



/*
//Basic mutable const-correct vector class with lots of operator overloading!

template<typename Float,int DIM>
class VectorND {
    static_assert(0<=DIM);
public:
    std::array<Float,DIM> x;
    VectorND();
    VectorND(Float y[DIM]);

    //multiplication/division on righthand side by scalar.
    //Division by zero throws an exception as expected.
    VectorND operator*  (Float b) const;
    VectorND& operator*=    (Float b);
    VectorND operator/  (Float b) const;
    VectorND& operator/=    (Float b);

    //calculate length and mutate to (x/length,y/length).
    void normalize();

    //return normalized version of the vector.
    VectorND normalized() const;

    //return the length squared of the vector
    Float length2() const;

    //return the length of the vector
    Float length() const;
};

//Vector left multiplication
template<typename Float,int DIM>
inline VectorND<Float,DIM> operator* (Float b,const VectorND<Float,DIM>& c);

//Writing a vector to an ostream as "(x,y)"
template<typename Float,int DIM>
inline std::ostream& operator<<(std::ostream& out,const VectorND<Float,DIM>& c);

//Linearly interpolate between vec1 and vec2.
template<typename Float,int DIM>
inline VectorND<Float,DIM> lerp(Float t,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2);

//cubic interpolation between v1 and v2 with control points v0 v3.
template<typename Float,int DIM>
inline VectorND<Float,DIM> cerp(Float t,const VectorND<Float,DIM>& vec0,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2,const VectorND<Float,DIM>& vec3);

template<typename Float,int DIM>
inline VectorND<Float,DIM>::VectorND<Float,DIM>() : {
    for(int i=0;i<DIM;i++)
        x[i]=Float(0);
}
template<typename Float,int DIM>
inline VectorND<Float,DIM>::VectorND<Float,DIM>(std::array<Float,DIM> y) x:y {}
template<typename Float,int DIM>
inline VectorND<Float,DIM>::VectorND<Float,DIM>() : x(a),y(b) { }
template<typename Float,int DIM>
inline VectorND<Float,DIM> VectorND<Float,DIM>::operator+ (const VectorND<Float,DIM>& b) const {
    return VectorND<Float,DIM>(x+b.x,y+b.y);
}
template<typename Float,int DIM>
inline VectorND<Float,DIM>& VectorND<Float,DIM>::operator+= (const VectorND<Float,DIM>& b) {
    x=x+b.x;
    y=y+b.y;
    return *this;
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> VectorND<Float,DIM>::operator-(const VectorND<Float,DIM>& b) const {
    return VectorND<Float,DIM>(x-b.x,y-b.y);
}
template<typename Float,int DIM>
inline VectorND<Float,DIM>& VectorND<Float,DIM>::operator-=(const VectorND<Float,DIM>& b) {
    x=x-b.x;
    y=y-b.y;
    return *this;
}
template<typename Float,int DIM>
inline Float VectorND<Float,DIM>::operator* (const VectorND<Float,DIM>& b) const {
    return (x*b.x+y*b.y);
}
template<typename Float,int DIM>
inline Float VectorND<Float,DIM>::operator^(const VectorND<Float,DIM>& b) const {
    return (x*b.y-y*b.x);
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> VectorND<Float,DIM>::operator* (Float b) const {
    return VectorND<Float,DIM>(x*b,y*b);
}
template<typename Float,int DIM>
inline VectorND<Float,DIM>& VectorND<Float,DIM>::operator*= (Float b)
{
    x=x*b;
    y=y*b;
    return *this;
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> VectorND<Float,DIM>::operator/ (Float b) const {
    return VectorND<Float,DIM>(x/b,y/b);
}
template<typename Float,int DIM>
inline VectorND<Float,DIM>& VectorND<Float,DIM>::operator/= (Float b) {
    x=x/b;
    y=y/b;
    return *this;
}
template<typename Float,int DIM>
inline void VectorND<Float,DIM>::normalize() {
    Float dist = sqrt(x*x+y*y);
    x/=dist;
    y/=dist;
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> VectorND<Float,DIM>::normalized() const {
    Float dist = sqrt(x*x+y*y);
    return VectorND<Float,DIM>(x/dist,y/dist);
}
template<typename Float,int DIM>
inline Float VectorND<Float,DIM>::length() const {
    return sqrt(x*x+y*y);
}
template<typename Float,int DIM>
inline Float VectorND<Float,DIM>::length2() const {
    return x*x+y*y;
}

template<typename Float,int DIM>
inline std::ostream& operator<<(std::ostream& out,const VectorND<Float,DIM>& c) {
    out << '(' << c.x << ',' << c.y << ')';
    return out;
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> operator* (Float b,const VectorND<Float,DIM>& c) {
    VectorND<Float,DIM> temp;
    temp.x=c.x*b;
    temp.y=c.y*b;
    return temp;
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> lerp(Float t,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2) {
    return (vec2-vec1)*t+vec1;
}
template<typename Float,int DIM>
inline VectorND<Float,DIM> cerp(Float t,const VectorND<Float,DIM>& vec0,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2,const VectorND<Float,DIM>& vec3) {
    Float t2=t*t;
    VectorND<Float,DIM> a0=vec3-vec2-vec0+vec1;
    VectorND<Float,DIM> a1=vec0-vec1-a0;
    VectorND<Float,DIM> a2=vec2-vec0;
    return a0*t*t2+a1*t2+a2*t+vec1;
}
*/
#endif
