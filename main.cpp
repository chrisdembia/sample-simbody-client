// A mass-spring-damper-actuator system.

#include "Simbody.h"

// #define RIGID_CONTACT

using namespace SimTK;
using std::cout;
using std::endl;

class LinearActuator : public Force::Custom::Implementation
{
public:
    LinearActuator(const MobilizedBody::Slider& mass) : m_mass(mass) {}
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const
        override
    {
        if (0.5 < state.getTime() && state.getTime() < 1.0)
        {
            m_mass.applyOneMobilityForce(state, 0, -300, mobilityForces);
        }
    }

    Real calcPotentialEnergy(const State& state) const override { return 0; }

private:
    const MobilizedBody::Slider& m_mass;
};

int main() {

    // Parameters.
    Real g = 9.8;
    Real m = 1.0;
    Real k = 1000.0;
    Real d = 10.0;
    Real L = 1.0;
    Real F = 100;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Matter.
    Body::Rigid massless(MassProperties(0.001, Vec3(0), Inertia(0)));
    Body::Rigid body(MassProperties(m, Vec3(0), Inertia(0)));

    Transform jointTransform(Rotation(Pi / 2, ZAxis));
    MobilizedBody::Slider platform(matter.Ground(), jointTransform,
                                   massless,        jointTransform);
    MobilizedBody::Slider mass(platform, jointTransform,
                               body,     jointTransform);

    // Forces.
    Force::Gravity gravity(forces, matter, -YAxis, g);
#ifndef RIGID_CONTACT
    Force::MobilityLinearStop stop(forces, platform, MobilizerQIndex(0),
            1000, 10, 0);
#endif
    Force::MobilityLinearSpring spring(forces, mass, MobilizerQIndex(0), k, L);
    Force::MobilityLinearDamper damper(forces, mass, MobilizerUIndex(0), d);
    Force::Custom actuator(forces, new LinearActuator(mass));

#ifdef RIGID_CONTACT
    // Constraints.
    HardStopLower* stop = new HardStopLower(platform, MobilizerQIndex(0),
            0.0, 0.0);
    matter.adoptUnilateralContact(stop);
#endif

    // Visualization.
    Visualizer viz(system);
    viz.setShowSimTime(true);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));

    // Initialize the system and state.
    State state = system.realizeTopology();
    // platform.setLength(state, 1.0);
    Real restLength = L - m * g / k;
    mass.setLength(state, restLength);

#ifdef RIGID_CONTACT
    SemiExplicitEulerTimeStepper ts(system);
    // TODO ts.setImpulseSolverType(SemiExplicitEulerTimeStepper::PLUS);
#else
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
#endif
    ts.initialize(state);
    ts.stepTo(20.0);

    return 0;
}

