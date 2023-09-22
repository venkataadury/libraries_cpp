#include <cmath>
#include "Eigen337/Eigen/Dense" //Can be removed
#include "Eigen337/Eigen/Core" //Can be removed (included in Atom.hpp)
namespace phases
{
    static constexpr double TWOPI=2*M_PI;
    class Planet
    {
    public:
        double position; //angle
        const double vel; //angle per second

        Planet(float v,float p=0.0f) : vel(v) {position=p;}

        //time in seconds
        inline void step(long time)
        {
            position+=vel*time;
            while(position>=TWOPI) position-=TWOPI;
            while(position<0) position+=TWOPI;
        }
        inline double getPosition() const {return position;}
        inline double getVelocity() const {return vel;}

        inline int asPhase(int phases) const
        {
            const double divs=TWOPI/phases;
            return (int)(std::round(position/divs))%phases;
        }
    };

    class Ellipse //Center is always at the origin
    {
        Eigen::Vector3d a,b;
        double c;
    public:
        Ellipse() {}
        Ellipse(const Eigen::Vector3d& v1, const Eigen::Vector3d v2)
        {
            a=v1;
            b=v2;
            c=sqrt(v1.norm()-v2.norm());
        }

        inline double at(double theta) const {return a*cos(theta)+b*sin(theta);}
        inline double operator[](double theta) const {return this->at(theta);}
        inline Eigen::Vector3d getCorrectionVector() const {return -c*(a/a.norm());}
    };

    class CelestialBody
    {
    protected:
        Eigen::Vector3d pos;
        CelestialBody() {}
    };

    class OrbitingBody : public CelestialBody
    {
    public:
        CelestialBody* relative=nullptr; //This is the focal body (at the FOCUS, not centre)
        Ellipse path;
        Eigen::Vector3d correction; //Correction from focal body to get centre of ellipse
        double vext; // Velocity at extreme (at r=a)
    public:
        OrbitingBody(CelestialBody* around, Eigen::Vector3d a, Eigen::Vector3d b, double vel) : CelestialBody()
        {
            relative=around;
            path=Ellipse(a,b);
        }
    };
}
