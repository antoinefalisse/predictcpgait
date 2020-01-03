#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
#include <OpenSim/Simulation/Model/MovingPathPoint.h>
#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
#include <recorder.hpp>

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;
using namespace OpenSim;

// Declare inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 21;        // degrees of freedom
constexpr int NX = ndof*2;      // states
constexpr int NU = ndof;        // controls
constexpr int NR = ndof+6+6+8;  // residual torques + GRFs + MOMs.

/// Value
template<typename T>
T value(const Recorder& e) { return e; }

template<>
double value(const Recorder& e) { return e.getValue(); }

/* Function F, using templated type T
F(x,u) -> (tau)
*/

template<typename T>
int F_generic(const T** arg, T** res) {

    // OpenSim model
    /// bodies
    OpenSim::Model* model;
    OpenSim::Body* pelvis;
    OpenSim::Body* femur_r;
    OpenSim::Body* femur_l;
    OpenSim::Body* tibia_r;
    OpenSim::Body* tibia_l;
    OpenSim::Body* talus_r;
    OpenSim::Body* talus_l;
    OpenSim::Body* calcn_r;
    OpenSim::Body* calcn_l;
    OpenSim::Body* toes_r;
    OpenSim::Body* toes_l;
    OpenSim::Body* torso;
    /// joints
    OpenSim::CustomJoint* ground_pelvis;
    OpenSim::CustomJoint* hip_r;
    OpenSim::CustomJoint* hip_l;
    OpenSim::CustomJoint* knee_r;
    OpenSim::CustomJoint* knee_l;
    OpenSim::CustomJoint* ankle_r;
    OpenSim::CustomJoint* ankle_l;
    OpenSim::CustomJoint* subtalar_r;
    OpenSim::CustomJoint* subtalar_l;
    OpenSim::WeldJoint* mtp_r;
    OpenSim::WeldJoint* mtp_l;
    OpenSim::CustomJoint* back;
    /// contact elements
    OpenSim::HuntCrossleyForce_smooth* HC_1_r;
    OpenSim::HuntCrossleyForce_smooth* HC_2_r;
    OpenSim::HuntCrossleyForce_smooth* HC_3_r;
    OpenSim::HuntCrossleyForce_smooth* HC_4_r;
    OpenSim::HuntCrossleyForce_smooth* HC_1_l;
    OpenSim::HuntCrossleyForce_smooth* HC_2_l;
    OpenSim::HuntCrossleyForce_smooth* HC_3_l;
    OpenSim::HuntCrossleyForce_smooth* HC_4_l;
    /// Model
    model = new OpenSim::Model();
    /// Bodies - Definition
    pelvis = new OpenSim::Body("pelvis", 3.69727, Vec3(-0.0737526, -0.0226889, 0), Inertia(0.0238289107, 0.0270444653, 0.029685404, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 4.64393, Vec3(0, -0.159112, 0), Inertia(0.078825664, 0.0161678053, 0.078825664, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 4.64393, Vec3(0, -0.161537, 0), Inertia(0.078825664, 0.0161678053, 0.078825664, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 1.43323, Vec3(0, -0.165359, 0), Inertia(0.0132395587, 0.002275956, 0.0137832813, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 1.43323, Vec3(0, -0.162479, 0), Inertia(0.0132395587, 0.002275956, 0.0137832813, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.01986, Vec3(0.00464602, 0.00194288, 0), Inertia(1.324e-006, 4.413e-007, 1.324e-006, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.01986, Vec3(0.00459638, 0.00192212, 0), Inertia(1.324e-006, 4.413e-007, 1.324e-006, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 0.39058, Vec3(0.0860486, 0.0131778, 0), Inertia(0.000401172, 0.0010212453, 0.0011840973, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 0.39058, Vec3(0.085129, 0.013037, 0), Inertia(0.000401172, 0.0010212453, 0.0011840973, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.04303, Vec3(0.0259333, -0.0021963, -0.00887), Inertia(4.41333e-005, 0.00011254, 0.000129752, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.04303, Vec3(0.0256561, -0.00217283, 0.008775), Inertia(4.41333e-005, 0.00011254, 0.000129752, 0, 0, 0));
    torso = new OpenSim::Body("torso", 16.25541, Vec3(0.0123281, 0.218025, 0), Inertia(0.215542616846467, 0.135378536674077, 0.291668449684595, 0, 0, 0));
    /// Joints - Transforms
    // Ground-Pelvis
    SpatialTransform st_ground_pelvis;
    st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("lower_torso_RX", 1, 1));
    st_ground_pelvis[0].setFunction(new LinearFunction());
    st_ground_pelvis[0].setAxis(Vec3(1, 0, 0));
    st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("lower_torso_RY", 1, 1));
    st_ground_pelvis[1].setFunction(new LinearFunction());
    st_ground_pelvis[1].setAxis(Vec3(0, 1, 0));
    st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("lower_torso_RZ", 1, 1));
    st_ground_pelvis[2].setFunction(new LinearFunction());
    st_ground_pelvis[2].setAxis(Vec3(0, 0, 1));
    st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("lower_torso_TX", 1, 1));
    st_ground_pelvis[3].setFunction(new LinearFunction());
    st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
    st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("lower_torso_TY", 1, 1));
    st_ground_pelvis[4].setFunction(new LinearFunction());
    st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));
    st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("lower_torso_TZ", 1, 1));
    st_ground_pelvis[5].setFunction(new LinearFunction());
    st_ground_pelvis[5].setAxis(Vec3(0, 0, 1));
    // Hip_l
    SpatialTransform st_hip_l;
    st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
    st_hip_l[0].setFunction(new LinearFunction());
    st_hip_l[0].setAxis(Vec3(0, 0, 1));
    st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
    st_hip_l[1].setFunction(new LinearFunction());
    st_hip_l[1].setAxis(Vec3(-1, 0, 0));
    st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
    st_hip_l[2].setFunction(new LinearFunction());
    st_hip_l[2].setAxis(Vec3(0, -1, 0));
    // Hip_r
    SpatialTransform st_hip_r;
    st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
    st_hip_r[0].setFunction(new LinearFunction());
    st_hip_r[0].setAxis(Vec3(0, 0, 1));
    st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
    st_hip_r[1].setFunction(new LinearFunction());
    st_hip_r[1].setAxis(Vec3(1, 0, 0));
    st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
    st_hip_r[2].setFunction(new LinearFunction());
    st_hip_r[2].setAxis(Vec3(0, 1, 0));
    // Knee_l
    SpatialTransform st_knee_l;
    st_knee_l[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    st_knee_l[2].setFunction(new LinearFunction());
    st_knee_l[2].setAxis(Vec3(0, 0, -1));
    // Knee_r
    SpatialTransform st_knee_r;
    st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[2].setFunction(new LinearFunction());
    st_knee_r[2].setAxis(Vec3(0, 0, -1));
    // Ankle_l
    SpatialTransform st_ankle_l;
    st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
    st_ankle_l[0].setFunction(new LinearFunction());
    st_ankle_l[0].setAxis(Vec3(0.104529047112, 0.173649078266, 0.979244441356));
    // Ankle_r
    SpatialTransform st_ankle_r;
    st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
    st_ankle_r[0].setFunction(new LinearFunction());
    st_ankle_r[0].setAxis(Vec3(-0.104529047112, -0.173649078266, 0.979244441356));
    // Subtalar_l
    SpatialTransform st_subtalar_l;
    st_subtalar_l[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_l", 1, 1));
    st_subtalar_l[0].setFunction(new LinearFunction());
    st_subtalar_l[0].setAxis(Vec3(-0.787180020856, -0.604747016023, -0.120949003205));
    // Subtalar_r
    SpatialTransform st_subtalar_r;
    st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
    st_subtalar_r[0].setFunction(new LinearFunction());
    st_subtalar_r[0].setAxis(Vec3(0.787180020856, 0.604747016023, -0.120949003205));
    // Back
    SpatialTransform st_back;
    st_back[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_extension", 1, 1));
    st_back[0].setFunction(new LinearFunction());
    st_back[0].setAxis(Vec3(0, 0, 1));
    st_back[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_bending", 1, 1));
    st_back[1].setFunction(new LinearFunction());
    st_back[1].setAxis(Vec3(1, 0, 0));
    st_back[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_rotation", 1, 1));
    st_back[2].setFunction(new LinearFunction());
    st_back[2].setAxis(Vec3(0, 1, 0));
    /// Joints - Definition
    ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.0384699, -0.0708953, -0.0646186), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.0443968, -0.0673054, 0.0671603), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *femur_l, Vec3(0.00000009, -0.393218, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0.105096325674, 0.048462535877, -0.005110024376), st_knee_l);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(0.0000001, -0.394719, -0.00000002), Vec3(0), *tibia_r, Vec3(0), Vec3(-0.064206849712, -0.022848068757, -0.001468892896), st_knee_r);
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.386280874342227, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.386028856205262, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
    subtalar_l = new CustomJoint("subtalar_l", *talus_l, Vec3(-0.041198, -0.035436, -0.00669), Vec3(0), *calcn_l, Vec3(0), Vec3(0),st_subtalar_l);
    subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.040757, -0.035058, 0.006619), Vec3(0), *calcn_r, Vec3(0), Vec3(0),st_subtalar_r);
    mtp_l = new WeldJoint("mtp_l", *calcn_l, Vec3(0.149349, -0.001689, -0.000912), Vec3(0), *toes_l, Vec3(0), Vec3(0));
    mtp_r = new WeldJoint("mtp_r", *calcn_r, Vec3(0.147753, -0.001671, 0.000903), Vec3(0), *toes_r, Vec3(0), Vec3(0));
    back = new CustomJoint("back", *pelvis, Vec3(-0.072018, 0.029113, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
    /// bodies and joints
    model->addBody(pelvis);		    model->addJoint(ground_pelvis);
    model->addBody(femur_l);		model->addJoint(hip_l);
    model->addBody(femur_r);		model->addJoint(hip_r);
    model->addBody(tibia_l);		model->addJoint(knee_l);
    model->addBody(tibia_r);		model->addJoint(knee_r);
    model->addBody(talus_l);		model->addJoint(ankle_l);
    model->addBody(talus_r);		model->addJoint(ankle_r);
    model->addBody(calcn_l);		model->addJoint(subtalar_l);
    model->addBody(calcn_r);		model->addJoint(subtalar_r);
    model->addBody(toes_l);		    model->addJoint(mtp_l);
    model->addBody(toes_r);		    model->addJoint(mtp_r);
    model->addBody(torso);          model->addJoint(back);
    /// Contact elements
    /// Parameters
    osim_double_adouble radiusSphere = 0.032;
    osim_double_adouble stiffness = 1000000;
    osim_double_adouble dissipation = 2.0;
    osim_double_adouble staticFriction = 0.8;
    osim_double_adouble dynamicFriction = 0.8;
    osim_double_adouble viscousFriction = 0.5;
    osim_double_adouble transitionVelocity = 0.2;
    Vec3 normal = Vec3(0, 1, 0);
    osim_double_adouble offset = 0;
    Vec3 locSphere_s1_r(-0.023098,   -0.021859, -0.0029017);
    Vec3 locSphere_s2_r(0.16997,     -0.021859, -0.0038755);
    Vec3 locSphere_s3_r(0.10809,     -0.021859,  0.055221);
    Vec3 locSphere_s4_r(0.0912,      -0.021859,  0.051336);
    Vec3 locSphere_s1_l(-0.023098,   -0.021859,  0.0029017);
    Vec3 locSphere_s2_l(0.16997,     -0.021859,  0.0038755);
    Vec3 locSphere_s3_l(0.10809,     -0.021859, -0.055221);
    Vec3 locSphere_s4_l(0.0912,      -0.021859, -0.051336);
    /// Left foot contact shere specifications
    HC_1_l = new HuntCrossleyForce_smooth("sphere_1_l", "calcn_l", locSphere_s1_l, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_2_l = new HuntCrossleyForce_smooth("sphere_2_l", "calcn_l", locSphere_s2_l, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_3_l = new HuntCrossleyForce_smooth("sphere_3_l", "calcn_l", locSphere_s3_l, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_4_l = new HuntCrossleyForce_smooth("sphere_4_l", "calcn_l", locSphere_s4_l, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    /// Add left foot contact spheres to model
    model->addComponent(HC_1_l);
    HC_1_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_2_l);
    HC_2_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_3_l);
    HC_3_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_4_l);
    HC_4_l->connectSocket_body_sphere(*calcn_l);
    /// Right foot contact shere specifications
    HC_1_r = new HuntCrossleyForce_smooth("sphere_1_r", "calcn_r", locSphere_s1_r, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_2_r = new HuntCrossleyForce_smooth("sphere_2_r", "calcn_r", locSphere_s2_r, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_3_r = new HuntCrossleyForce_smooth("sphere_3_r", "calcn_r", locSphere_s3_r, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_4_r = new HuntCrossleyForce_smooth("sphere_4_r", "calcn_r", locSphere_s4_r, radiusSphere,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    /// Add right foot contact spheres to model
    model->addComponent(HC_1_r);
    HC_1_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_2_r);
    HC_2_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_3_r);
    HC_3_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_4_r);
    HC_4_r->connectSocket_body_sphere(*calcn_r);
    /// Initialize  system and state.
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);

    // States and controls
    T ua[NU]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];    // states
    for (int i = 0; i < 12; ++i) ua[i] = u[i];      // controls
    ua[12] = u[18]; // 12 Simbody (lumbar-ext) is 18 OpenSim
    ua[13] = u[19]; // 13 Simbody (lumbar-bend) is 19 OpenSim
    ua[14] = u[20]; // 14 Simbody (lumbar-rot) is 20 OpenSim
    ua[15] = u[12]; // 15 Simbody (knee-angle-l) is 12 OpenSim
    ua[16] = u[13]; // 16 Simbody (knee-angle-r) is 13 OpenSim
    ua[17] = u[14]; // 17 Simbody (ankle-angle-l) is 14 OpenSim
    ua[18] = u[15]; // 18 Simbody (ankle-angle-l) is 15 OpenSim
    ua[19] = u[16]; // 19 Simbody (subtalar-angle-l) is 16 OpenSim
    ua[20] = u[17]; // 20 Simbody (subtalar-angle-l) is 17 OpenSim

    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    // Residual forces
    /// appliedMobilityForces (# mobilities)
    Vector appliedMobilityForces(ndof);
    appliedMobilityForces.setToZero();
    /// appliedBodyForces (# bodies + ground)
    Vector_<SpatialVec> appliedBodyForces;
    int nbodies = model->getBodySet().getSize() + 1; // including ground
    appliedBodyForces.resize(nbodies);
    appliedBodyForces.setToZero();
    /// Gravity
    Vec3 gravity(0);
    gravity[1] = -9.81;
    /// Weight
    for (int i = 0; i < model->getBodySet().getSize(); ++i) {
        model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
            model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
    }
    /// Add contact forces to appliedBodyForces
    /// Right foot
    Array<osim_double_adouble> Force_values_1_r = HC_1_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_2_r = HC_2_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_3_r = HC_3_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_4_r = HC_4_r->getRecordValues(*state);
    SpatialVec GRF_1_r;
    GRF_1_r[0] = Vec3(Force_values_1_r[9], Force_values_1_r[10], Force_values_1_r[11]);
    GRF_1_r[1] = Vec3(Force_values_1_r[6], Force_values_1_r[7], Force_values_1_r[8]);
    SpatialVec GRF_2_r;
    GRF_2_r[0] = Vec3(Force_values_2_r[9], Force_values_2_r[10], Force_values_2_r[11]);
    GRF_2_r[1] = Vec3(Force_values_2_r[6], Force_values_2_r[7], Force_values_2_r[8]);
    SpatialVec GRF_3_r;
    GRF_3_r[0] = Vec3(Force_values_3_r[9], Force_values_3_r[10], Force_values_3_r[11]);
    GRF_3_r[1] = Vec3(Force_values_3_r[6], Force_values_3_r[7], Force_values_3_r[8]);
    SpatialVec GRF_4_r;
    GRF_4_r[0] = Vec3(Force_values_4_r[9], Force_values_4_r[10], Force_values_4_r[11]);
    GRF_4_r[1] = Vec3(Force_values_4_r[6], Force_values_4_r[7], Force_values_4_r[8]);
    int ncalcn_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
    appliedBodyForces[ncalcn_r] = appliedBodyForces[ncalcn_r] + GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r;
    /// Left foot
    Array<osim_double_adouble> Force_values_1_l = HC_1_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_2_l = HC_2_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_3_l = HC_3_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_4_l = HC_4_l->getRecordValues(*state);
    SpatialVec GRF_1_l;
    GRF_1_l[0] = Vec3(Force_values_1_l[9], Force_values_1_l[10], Force_values_1_l[11]);
    GRF_1_l[1] = Vec3(Force_values_1_l[6], Force_values_1_l[7], Force_values_1_l[8]);
    SpatialVec GRF_2_l;
    GRF_2_l[0] = Vec3(Force_values_2_l[9], Force_values_2_l[10], Force_values_2_l[11]);
    GRF_2_l[1] = Vec3(Force_values_2_l[6], Force_values_2_l[7], Force_values_2_l[8]);
    SpatialVec GRF_3_l;
    GRF_3_l[0] = Vec3(Force_values_3_l[9], Force_values_3_l[10], Force_values_3_l[11]);
    GRF_3_l[1] = Vec3(Force_values_3_l[6], Force_values_3_l[7], Force_values_3_l[8]);
    SpatialVec GRF_4_l;
    GRF_4_l[0] = Vec3(Force_values_4_l[9], Force_values_4_l[10], Force_values_4_l[11]);
    GRF_4_l[1] = Vec3(Force_values_4_l[6], Force_values_4_l[7], Force_values_4_l[8]);
    int ncalcn_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
    appliedBodyForces[ncalcn_l] = appliedBodyForces[ncalcn_l] + GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l;
    // Ground reaction forces
    SpatialVec GRF_r = GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r;
    SpatialVec GRF_l = GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l;

    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) {
        knownUdot[i] = ua[i];
    }
    // Residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Calculate contact torques about the ground origin
    // Step: calculate contact point positions in body frames
    /// sphere 1 left
    Vec3 pos_InGround_HC_s1_l = calcn_l->findStationLocationInGround(*state, locSphere_s1_l);
    Vec3 contactPointpos_InGround_HC_s1_l = pos_InGround_HC_s1_l - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s1_l_adj = contactPointpos_InGround_HC_s1_l - 0.5*contactPointpos_InGround_HC_s1_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s1_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_l_adj, *calcn_l);
    /// sphere 2 left
    Vec3 pos_InGround_HC_s2_l = calcn_l->findStationLocationInGround(*state, locSphere_s2_l);
    Vec3 contactPointpos_InGround_HC_s2_l = pos_InGround_HC_s2_l - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s2_l_adj = contactPointpos_InGround_HC_s2_l - 0.5*contactPointpos_InGround_HC_s2_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s2_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_l_adj, *calcn_l);
    /// sphere 3 left
    Vec3 pos_InGround_HC_s3_l = calcn_l->findStationLocationInGround(*state, locSphere_s3_l);
    Vec3 contactPointpos_InGround_HC_s3_l = pos_InGround_HC_s3_l - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s3_l_adj = contactPointpos_InGround_HC_s3_l - 0.5*contactPointpos_InGround_HC_s3_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s3_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_l_adj, *calcn_l);
    /// sphere 4 left
    Vec3 pos_InGround_HC_s4_l = calcn_l->findStationLocationInGround(*state, locSphere_s4_l);
    Vec3 contactPointpos_InGround_HC_s4_l = pos_InGround_HC_s4_l - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s4_l_adj = contactPointpos_InGround_HC_s4_l - 0.5*contactPointpos_InGround_HC_s4_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s4_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_l_adj, *calcn_l);
    /// sphere 1 right
    Vec3 pos_InGround_HC_s1_r = calcn_r->findStationLocationInGround(*state, locSphere_s1_r);
    Vec3 contactPointpos_InGround_HC_s1_r = pos_InGround_HC_s1_r - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s1_r_adj = contactPointpos_InGround_HC_s1_r - 0.5*contactPointpos_InGround_HC_s1_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s1_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_r_adj, *calcn_r);
    /// sphere 2 right
    Vec3 pos_InGround_HC_s2_r = calcn_r->findStationLocationInGround(*state, locSphere_s2_r);
    Vec3 contactPointpos_InGround_HC_s2_r = pos_InGround_HC_s2_r - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s2_r_adj = contactPointpos_InGround_HC_s2_r - 0.5*contactPointpos_InGround_HC_s2_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s2_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_r_adj, *calcn_r);
    /// sphere 3 right
    Vec3 pos_InGround_HC_s3_r = calcn_r->findStationLocationInGround(*state, locSphere_s3_r);
    Vec3 contactPointpos_InGround_HC_s3_r = pos_InGround_HC_s3_r - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s3_r_adj = contactPointpos_InGround_HC_s3_r - 0.5*contactPointpos_InGround_HC_s3_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s3_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_r_adj, *calcn_r);
    /// sphere 4 right
    Vec3 pos_InGround_HC_s4_r = calcn_r->findStationLocationInGround(*state, locSphere_s4_r);
    Vec3 contactPointpos_InGround_HC_s4_r = pos_InGround_HC_s4_r - radiusSphere*normal;
    Vec3 contactPointpos_InGround_HC_s4_r_adj = contactPointpos_InGround_HC_s4_r - 0.5*contactPointpos_InGround_HC_s4_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s4_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_r_adj, *calcn_r);
    // Contact forces
    Vec3 AppliedPointForce_s1_r = GRF_1_r[1];
    Vec3 AppliedPointForce_s2_r = GRF_2_r[1];
    Vec3 AppliedPointForce_s3_r = GRF_3_r[1];
    Vec3 AppliedPointForce_s4_r = GRF_4_r[1];
    Vec3 AppliedPointForce_s1_l = GRF_1_l[1];
    Vec3 AppliedPointForce_s2_l = GRF_2_l[1];
    Vec3 AppliedPointForce_s3_l = GRF_3_l[1];
    Vec3 AppliedPointForce_s4_l = GRF_4_l[1];
    /// Get transforms
    SimTK::Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);
    /// Calculate torques
    Vec3 AppliedPointTorque_s1_l, AppliedPointTorque_s2_l, AppliedPointTorque_s3_l, AppliedPointTorque_s4_l;
    Vec3 AppliedPointTorque_s1_r, AppliedPointTorque_s2_r, AppliedPointTorque_s3_r, AppliedPointTorque_s4_r;
    AppliedPointTorque_s1_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s1_l) % AppliedPointForce_s1_l;
    AppliedPointTorque_s2_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s2_l) % AppliedPointForce_s2_l;
    AppliedPointTorque_s3_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s3_l) % AppliedPointForce_s3_l;
    AppliedPointTorque_s4_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s4_l) % AppliedPointForce_s4_l;
    AppliedPointTorque_s1_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s1_r) % AppliedPointForce_s1_r;
    AppliedPointTorque_s2_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s2_r) % AppliedPointForce_s2_r;
    AppliedPointTorque_s3_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s3_r) % AppliedPointForce_s3_r;
    AppliedPointTorque_s4_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s4_r) % AppliedPointForce_s4_r;
    /// Contact torques
    Vec3 MOM_l, MOM_r;
    MOM_l = AppliedPointTorque_s1_l + AppliedPointTorque_s2_l + AppliedPointTorque_s3_l + AppliedPointTorque_s4_l;
    MOM_r = AppliedPointTorque_s1_r + AppliedPointTorque_s2_r + AppliedPointTorque_s3_r + AppliedPointTorque_s4_r;

    // Extract several joint origins to set constraints in problem
    Vec3 calcn_or_l  = calcn_l->getPositionInGround(*state);
    Vec3 calcn_or_r  = calcn_r->getPositionInGround(*state);
    Vec3 tibia_or_l  = tibia_l->getPositionInGround(*state);
    Vec3 tibia_or_r  = tibia_r->getPositionInGround(*state);

     //CasADi may not always request all outputs
     //if res[i] is a null pointer, this means that output i is not required
    int nc = 3; // # components in Vec3
    if (res[0]) {
        for (int i = 0; i < 12; ++i) {
            res[0][i] = value<T>(residualMobilityForces[i]); // residual torques
        }
        // order adjusted since order Simbody is different than order OpenSim
        res[0][12] = value<T>(residualMobilityForces[15]);
        res[0][13] = value<T>(residualMobilityForces[16]);
        res[0][14] = value<T>(residualMobilityForces[17]);
        res[0][15] = value<T>(residualMobilityForces[18]);
        res[0][16] = value<T>(residualMobilityForces[19]);
        res[0][17] = value<T>(residualMobilityForces[20]);
        res[0][18] = value<T>(residualMobilityForces[12]);
        res[0][19] = value<T>(residualMobilityForces[13]);
        res[0][20] = value<T>(residualMobilityForces[14]);
        for (int i = 0; i < nc; ++i) {
            res[0][i + 21] = value<T>(GRF_r[1][i]); // GRF_r
        }
        for (int i = 0; i < nc; ++i) {
            res[0][i + 21 + nc] = value<T>(GRF_l[1][i]); // GRF_l
        }
        for (int i = 0; i < nc; ++i) {
            res[0][i + 21 + nc + nc] = value<T>(MOM_r[i]); // MOM_r
        }
        for (int i = 0; i < nc; ++i) {
            res[0][i + 21 + nc + nc + nc] = value<T>(MOM_l[i]); // MOM_l
        }
        res[0][21+4*nc] = value<T>(calcn_or_r[0]);  /// calcn_or_r_x
        res[0][22+4*nc] = value<T>(calcn_or_r[2]);  /// calcn_or_r_z
        res[0][23+4*nc] = value<T>(calcn_or_l[0]);  /// calcn_or_l_x
        res[0][24+4*nc] = value<T>(calcn_or_l[2]);  /// calcn_or_l_x
        res[0][25+4*nc] = value<T>(tibia_or_r[0]);  /// tibia_or_r_x
        res[0][26+4*nc] = value<T>(tibia_or_r[2]);  /// tibia_or_r_z
        res[0][27+4*nc] = value<T>(tibia_or_l[0]);  /// tibia_or_l_x
        res[0][28+4*nc] = value<T>(tibia_or_l[2]);  /// tibia_or_l_z
    }
    return 0;
}

int main() {

    Recorder x[NX];
    Recorder u[NU];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    for (int i = 0; i < NU; ++i) u[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x,u };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}