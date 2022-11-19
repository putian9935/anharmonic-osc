#include "gmpxx.h"

#if defined(O0)
int const m = 0;
mpq_class hermite_coeff[] = {1};
#elif defined(O1)
int const m = 1;
mpq_class hermite_coeff[] = {0, 2};
#elif defined(O2)
int const m = 2;
mpq_class hermite_coeff[] = {-2, 0, 4};
#elif defined(O3)
int const m = 3;
mpq_class hermite_coeff[] = {0, -12, 0, 8};
#elif defined(O4)
int const m = 4;
mpq_class hermite_coeff[] = {12, 0, -48, 0, 16};
#elif defined(O5)
int const m = 5;
mpq_class hermite_coeff[] = {0, 120, 0, -160, 0, 32};
#elif defined(O6)
int const m = 6;
mpq_class hermite_coeff[] = {-120, 0, 720, 0, -480, 0, 64};
#elif defined(O7)
int const m = 7;
mpq_class hermite_coeff[] = {0, -1680, 0, 3360, 0, -1344, 0, 128};
#elif defined(O8)
int const m = 8;
mpq_class hermite_coeff[] = {1680, 0, -13440, 0, 13440, 0, -3584, 0, 256};
#elif defined(O9)
int const m = 9;
mpq_class hermite_coeff[] = {0, 30240, 0, -80640, 0, 48384, 0, -9216, 0, 512};
#endif