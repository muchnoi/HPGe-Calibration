#include <memory>

#include "Math/SpecFunc.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/Interpolator.h"

namespace RM = ROOT::Math;

namespace {
    void dummy_gsl_handler (const char * , const char * , int , int ) { }
}

double airy_ai_int(const double *xx, const double *)
{
  static const double xmin  = -50.0;
  static const double xmax  =  6.00;
  static const double xstep =  0.05;

  static std::auto_ptr<RM::Interpolator> spline;
  if( !spline.get() ) {
    const RM::WrappedFunction<> wrapped(&RM::airy_Ai);

    RM::IntegratorOneDim engine(wrapped, RM::IntegrationOneDim::kADAPTIVE);

    const size_t nx = static_cast<size_t>((xmax - xmin)/xstep+1e-8);
    std::vector<double> x(nx), y(nx);

    for( size_t i = 0; i<nx; ++i ) {
      x[i] = xmin + i*xstep;
      y[i] = engine.Integral(0.0, x[i]);
    }
    spline.reset(new RM::Interpolator(x, y, RM::Interpolation::kAKIMA));
  }

  if( xx[0]      <xmin ) return 1;
  if( xx[0]+xstep>xmax ) return 0;

  return 1.0/3.0 - spline->Eval(xx[0]);
}
