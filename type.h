#ifndef ADAMS_TYPE_H
#define ADAMS_TYPE_H

#ifdef USE_LONG_DOUBLE
typedef long double DOUBLE;

#define fmod fmodl
#define sqrt sqrtl
#define pow powl
#else
typedef double DOUBLE;
#endif


#endif //ADAMS_TYPE_H
