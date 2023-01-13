#ifndef INCLUDED_CHEMBOX
#define INCLUDED_CHEMBOX 1
#include "draw/graphics.h"
#include "sci/sim.h"

namespace chembox
{
  class CuboidBox : public simtools::PBCCuboidBox
  {
    graphics::Frame myFrame;
    double timestep=1;
  public:
    CuboidBox(int w,int b,int h,int wW,int wH,const std::vector<simtools::Particle>& parts) : simtools::PBCCuboidBox(Vector(w,b,h),parts), myFrame(wW,wH,"System",true)
    {
    }
    
  };
}
#endif
