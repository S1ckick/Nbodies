#define NUMBER_DOUBLE_DOUBLE 1
//#define NUMBER_DOUBLE 1

#define SAVE_STEPS 1
#define SAVE_DIFF 1

//#define SAVE_INV 1

#define TAYLOR
#ifdef TAYLOR
  #define TAYLOR_TYPE ABMD_DOUBLE
#endif


const int YEARS = 40;

enum class Planets{
  SUN, MERCURY, VENUS, EARTH, 
  MARS, JUPITER, SATURN, URANUS, 
  NEPTUNE, PLUTO, MOON, CERES_1,
  PALLAS_1, VESTA_4, IRIS_7, BAMBERGA_324
};

const int DIFF_PLANET = int(Planets::MOON);

#ifdef NUMBER_DOUBLE_DOUBLE
  #include <qd/dd_real.h>
  #include <qd/fpu.h>

  using main_type = dd_real;
  using helper_type = double;
#else
  using main_type = long double;
  using helper_type = long double;
#endif

const int moonNum = int(Planets::MOON);
const int earthNum = int(Planets::EARTH);
const int barrier = int(Planets::BAMBERGA_324) + 1;