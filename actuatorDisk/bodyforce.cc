#line 1 "bodyforce.loci"
//#############################################################################
//#
//# Copyright 2019, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <iostream>
#include <Loci.h>
#include <flowTypes.h>
#include <flowPsiIO.h>
#line 1 "flowPsi.lh"
//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci>
#include <flowTypes.h>

#line 1 "FVM.lh"
//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

// $type pos store<Loci::vector3d<Loci::real_t> > 
// $type cl Map
// $type cr Map
// $type ci Map
// $type ref Map
// $type pmap Map
// $type face2node multiMap

// $type upper multiMap
// $type lower multiMap
// $type boundary_map multiMap

// $type cellcenter store<Loci::vector3d<Loci::real_t> > 
// $type facecenter store<Loci::vector3d<Loci::real_t> > 
// $type area store<Loci::Area> 
// $type vol store<Loci::real_t> 
// $type grid_vol param<Loci::real_t> 

// $type mn store<Loci::vector3d<Loci::real_t> > 
// $type ln store<Loci::vector3d<Loci::real_t> > 

// $type grads(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type grads_f(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv_f(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d_f(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type limiters(X0) store<Loci::real_t> 
// $type limiterv(X0) storeVec<Loci::real_t> 
// $type limiterv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type lefts(X0) store<Loci::real_t> 
// $type rights(X0) store<Loci::real_t> 
// $type leftsP(X0,X1) store<Loci::real_t> 
// $type rightsP(X0,X1) store<Loci::real_t> 
// $type leftvM(X0) storeVec<Loci::real_t> 
// $type rightvM(X0) storeVec<Loci::real_t> 
// $type leftv3d(X0) store<Loci::vector3d<Loci::real_t> > 
// $type rightv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type cell2node(X0) store<float> 
// $type cell2node_v(X0) storeVec<float> 
// $type cell2node_v3d(X0) store<Loci::vector3d<float> > 
// $type cell2nodeMax(X0) store<float> 
// $type cell2nodeMin(X0) store<float> 
// $type cell2nodeMaxMag(X0) store<float> 
// $type cell2nodeMaxv3d(X0) store<vector3d<float> > 

// $type BC_options store<Loci::options_list> 

// $type integrateSurface(X0) store<Loci::real_t> 
// $type integrateFlux(X0) store<Loci::real_t> 

// $type petscScalarSolve(X0) store<Loci::real_t> 
// $type petscBlockedSolve(X0) storeVec<Loci::real_t> 
// $type petscBlockedSSolve(X0) storeVec<Loci::real_t> 

// $type L2Norm(X0) param<Loci::real_t> 
// $type L1Norm(X0) param<Loci::real_t> 
// $type LinfNorm(X0) param<Loci::real_t> 

// $type fileNumber(X0) store<int> 
#line 24 "flowPsi.lh"


// $type plot_output param<flowPsi::list_input> 
// $type plot_output_exclusive param<flowPsi::list_input> 

// $type modelName param<std::string> 
// Transform for periodic bc's (should this be in FVM.lh?)
// $type periodicTransform store<Loci::rigid_transform> 


// EoS settings
// $type molecularMass param<flowPsi::real> 
// $type Rtilde param<flowPsi::real> 
// $type gamma param<flowPsi::real> 
// $type Cp param<flowPsi::real> 

// $type p0 param<PressureValue> 

// Initial Conditions
// $type temperature_ic store<flowPsi::real> 
// $type gagePressure_ic store<flowPsi::real> 
// $type u_ic store<flowPsi::vect3d> 
// Previous timestep values for BDF2
// $type temperature_ico store<flowPsi::real> 
// $type gagePressure_ico store<flowPsi::real> 
// $type u_ico store<flowPsi::vect3d> 

// Ambient Pressure used to define gage pressure
// $type Pambient param<flowPsi::real> 
// $type temperature store<flowPsi::real> 
// cell density
// $type rho store<flowPsi::real> 
// cell pressure
// $type pressure store<flowPsi::real> 
// $type e_internal store<flowPsi::real> 
// cell speed of sound

// $type soundSpeed store<flowPsi::real> 
// $type gagePressure store<flowPsi::real> 

// fluid velocity at a cell
// $type u store<flowPsi::vect3d> 
// grid velocity at face
// $type us_n store<flowPsi::real> 
// $type us store<flowPsi::vect3d> 
// $type gcl_sum store<flowPsi::real> 

// $type mu store<flowPsi::real>  // fluid viscosity cell
// $type mu_f store<flowPsi::real>  // fluid viscosity face
// $type mu(X0) store<flowPsi::real> 
// $type tmu store<flowPsi::real>  // turbulent viscosity cell
// $type tmu_f store<flowPsi::real>  // turbulent viscosity face
// $type muTotal store<flowPsi::real> 
// $type muTotal_f store<flowPsi::real> 

// $type kconduct store<flowPsi::real>  // fluid conductivity cell
// $type kconduct_f store<flowPsi::real>  // fluid conductivity face
// $type kconduct(X0) store<flowPsi::real> 
// $type turbulentPrandtlNumber param<flowPsi::real>  // Turbulent Prandtl number
// Contributions from fluxes to cell
// $type src store<Loci::Array<flowPsi::real,5> > 
// Right hand side of Newton method
// $type rhs store<Loci::Array<flowPsi::real,5> > 
// Inviscid Flux
// $type iflux store<Loci::Array<flowPsi::real,5> > 
// Viscous Flux
// $type vflux store<Loci::Array<flowPsi::real,4> > 

// Low Speed Preconditioning Work
// $type Mref store<flowPsi::real>  // Reference Mach number used for preconditioning
// $type PLimPC param<flowPsi::real> // Pressure jump limit on preconditioning
// $type Minf param<flowPsi::real>  // Free stream Mach number
// $type etaT param<flowPsi::real>  // Stagnation point limit factor
// $type Eta_p store<flowPsi::real>  // Eta factor used in preconditioning
// $type Eta_pf store<flowPsi::real>  // Eta factor used in preconditioning
// $type Eta_bc store<flowPsi::real>  // used for boundary conditions, only includes viscous terms
// Jacobian Assembly
// $type fjp storeMat<flowPsi::real_fj> 
// $type fjm storeMat<flowPsi::real_fj> 
// $type srcJ storeMat<flowPsi::real_fj> 

// Fluid Matrix, lower, upper, and diagonal components
// $type fluid_U storeMat<flowPsi::real_fj> 
// $type fluid_L storeMat<flowPsi::real_fj> 
// $type fluid_D storeMat<flowPsi::real_fj> 
// $type fluid_B storeVec<flowPsi::real_fj> 
// $type fluid_D_inv storeMat<flowPsi::real_fj> 
// $type fluid_pivot storeVec<pivot_type> 
// Face values
// $type gagePressure_f store<flowPsi::real> 
// $type pressure_f store<flowPsi::real> 
// $type temperature_f store<flowPsi::real> 
// $type rho_f store<flowPsi::real> 
// $type soundSpeed_f store<flowPsi::real> 
// $type u_f store<flowPsi::vect3d> 

// $type bc_total_area store<flowPsi::real> 
// $type bc_total_force store<flowPsi::real> 
// $type bc_average_pressure store<flowPsi::real> 

// $type flow_direction store<flowPsi::vect3d> 
// $type rigid_u store<flowPsi::vect3d> 
// $type p0Ref store<flowPsi::real> 
// $type T0Ref store<flowPsi::real> 
// $type massFluxRef store<flowPsi::real> 
// $type temperatureRef store<flowPsi::real> 
// $type gagePressureRef store<flowPsi::real> 
// $type uRef store<flowPsi::vect3d> 
// Solution Clipping
// $type TclipMin param<TemperatureValue> 
// $type TclipMax param<TemperatureValue> 
// $type PclipMin param<PressureValue> 
// Fluid Stress
// $type tau store<flowPsi::symmetricTensor> 
// Boundary Condition Parameters
// $type Uwall_BC store<flowPsi::vect3d> 
// $type T_BC store<flowPsi::real> 
// $type T0_BC store<flowPsi::real> 
// $type p0_BC store<flowPsi::real> 
// $type rho_BC store<flowPsi::real> 
// $type Twall_BC store<flowPsi::real> 
// $type qwall_BC store<flowPsi::real> 
// $type massFlux_BC store<flowPsi::real> 
// $type mdot_BC store<flowPsi::real> 
// $type swirlAngle_BC store<flowPsi::real> 
// $type swirlCenter_BC store<flowPsi::vect3d> 
// $type momentCenter_BC store<flowPsi::vect3d> 
// $type M_BC store<flowPsi::vect3d> 
// $type u_BC store<flowPsi::vect3d> 
// $type p_BC store<flowPsi::real> 
// $type pMean_BC store<flowPsi::real> 
// $type BC_inflow_type store<int> 
// $type temperatureRef_BC store<flowPsi::real> 
// $type gagePressureRef_BC store<flowPsi::real> 
// $type uRef_BC store<flowPsi::vect3d> 

// $type angVel_BC store<flowPsi::vect3d> 
// $type rotAxis_BC store<flowPsi::vect3d> 
// $type rotCenter_BC store<flowPsi::vect3d> 
// $type swirlAxis_BC store<flowPsi::vect3d> 
// $type flowDir_BC store<flowPsi::vect3d> 
// $type boundary_area store<flowPsi::real> 

// $type Twall store<flowPsi::real> 
// $type Twall_prescribed store<flowPsi::real> 
// $type T_interface store<flowPsi::real> 
// $type qwall store<flowPsi::real> 
// $type qwall_prescribed store<flowPsi::real> 
// $type wallVelocity store<flowPsi::vect3d> 
// $type heatf store<flowPsi::real> 
// $type wall_stress store<flowPsi::vect3d> 
// $type yplus_w store<flowPsi::real> 

// $type stime param<flowPsi::real>  // simulation time
// $type ncycle param<int>  // simulation iteration number
// $type newtonMaxIter param<int> 
// $type newtonFinished param<bool> 
// $type lastNewton param<bool> 
// $type urelax param<flowPsi::real> 
// $type dtmax param<flowPsi::TimeValue> 
// $type timeStepSize param<flowPsi::real> 
// $type cflmax param<flowPsi::real> 
// $type cflpdt store<flowPsi::real>  // Multiply this by dt to get cfl
// $type cfl param<std::pair<flowPsi::real,flowPsi::real> > 

// $type dtau store<flowPsi::real> 
// $type iblank store<flowPsi::byte_t> 

// General Turbulence Model
// Map from cell to nearest viscous wall
// $type min_cell2noslip Map
// $type multi_scale param<std::string> 

// $type divu store<flowPsi::real> 

// vorticity
// $type vortMag store<flowPsi::real> 
// $type vort store<flowPsi::vect3d> 
// $type vort_f store<flowPsi::vect3d> 

// strain rate
// $type strainRate store<flowPsi::real> 

// helicity
// $type helicity store<flowPsi::real> 
// $type helicity_f store<flowPsi::real> 

// distance to nearest viscous wall
// $type dist_noslip store<flowPsi::real> 
// $type dist_noslip_f store<flowPsi::real> 

// $type initialConditions param<Loci::options_list> 
// $type initialConditionRegions param<Loci::options_list> 
// $type icRegionInfo blackbox<std::vector<flowPsi::ICstate_info> > 
// $type icRegionId store<int> 

// $type print_freq param<int> 
// $type plot_freq param<int> 
// $type plot_modulo param<int> 
// $type do_plot param<bool> 
// $type plot_postfix param<std::string> 
// $type do_boundary_plot param<bool> 

// $type plot_status param<int> 
// $type restart_freq param<int> 
// $type restart_modulo param<int> 
// $type restart_directory param<std::string> 
// $type restart_postfix param<std::string> 
// $type do_restart param<bool> 
// $type icfile param<std::string> 

// $type maximumRunTime param<flowPsi::TimeValue> 


// $type grid_vol_iblank param<flowPsi::real> 
// $type integratedOutputFileManager blackbox<flowPsi::integratedFileDBManager> 
// $type do_output_integrate param<bool> 

// $type localPBias store<realF> 
// $type localPBias_next store<realF> 
// $type timeSteppingMode param<int> 

// $type gravityAccel param<flowPsi::vect3d> 

// Solver averaging types
// $type meanFreq param<int> 
// $type meanCount param<flowPsi::real> 
// $type meanCountReset param<flowPsi::real> 
// $type meanCount_ic param<pair<int,flowPsi::real> > 
// $type scalarMean(X0) store<flowPsi::real> 
// $type scalarMean_f(X0) store<flowPsi::real> 
// $type scalarMeanAll_f(X0) store<flowPsi::real> 
// $type scalarMeanFaceAll_X_ic store<flowPsi::real> 
// $type scalarM2(X0) store<flowPsi::real> 
// $type scalarVariance(X0) store<flowPsi::real> 
// $type scalarM2_f(X0) store<flowPsi::real> 
// $type scalarVariance_f(X0) store<flowPsi::real> 
// $type vect3dMean(X0) store<flowPsi::vect3d> 
// $type vect3dMean_f(X0) store<flowPsi::vect3d> 
// $type vect3dM2(X0) store<flowPsi::vect3d> 
// $type vect3dVariance(X0) store<flowPsi::vect3d> 
// $type vect3dM2_f(X0) store<flowPsi::vect3d> 
// $type vect3dVariance_f(X0) store<flowPsi::vect3d> 
// $type vect3dMX2(X0) store<flowPsi::vect3d> 
// $type vect3dCoVariance(X0) store<flowPsi::vect3d> 
// $type vect3dMX2_f(X0) store<flowPsi::vect3d> 
// $type vect3dCoVariance_f(X0) store<flowPsi::vect3d> 
// $type scalarFavreMean(X0) store<flowPsi::real> 
// $type scalarFavreMeanBase(X0) store<flowPsi::real> 
// $type scalarFavreMean_f(X0) store<flowPsi::real> 
// $type scalarFavreMeanBase_f(X0) store<flowPsi::real> 
// $type scalarFavreVar(X0) store<flowPsi::real> 
// $type scalarFavreVarBase(X0) store<flowPsi::real> 
// $type scalarFavreVar_f(X0) store<flowPsi::real> 
// $type scalarFavreVarBase_f(X0) store<flowPsi::real> 
// $type vect3dFavreMean(X0) store<flowPsi::vect3d> 
// $type vect3dFavreMeanBase(X0) store<flowPsi::vect3d> 
// $type vect3dFavreMean_f(X0) store<flowPsi::vect3d> 
// $type vect3dFavreMeanBase_f(X0) store<flowPsi::vect3d> 
// $type vect3dFavreVar(X0) store<flowPsi::vect3d> 
// $type vect3dFavreVarBase(X0) store<flowPsi::vect3d> 
// $type vect3dFavreVar_f(X0) store<flowPsi::vect3d> 
// $type vect3dFavreVarBase_f(X0) store<flowPsi::vect3d> 
// $type vect3dFavreCoVar(X0) store<flowPsi::vect3d> 
// $type vect3dFavreCoVarBase(X0) store<flowPsi::vect3d> 
// $type vect3dFavreCoVar_f(X0) store<flowPsi::vect3d> 
// $type vect3dFavreCoVarBase_f(X0) store<flowPsi::vect3d> 

// $type scalarTransportP(X0,X1) store<flowPsi::real> 
// $type scalarTransport(X0,X1) store<flowPsi::real> 

// $type timeStepSteadyState Constraint
// $type timeStepAccurate Constraint

// $type timeStepSchemeBDF2 Constraint
// $type timeStepSchemeKEC Constraint

// $type LaminarSimulation Constraint
// $type TurbulentSimulation Constraint
// $type ViscousSimulation Constraint


// $type fluidLinearSolverSGS Constraint
// $type fluidLinearSolverFSGS Constraint
// $type fluidLinearSolverLSGS Constraint
// $type fluidLinearSolverPETSC Constraint

// $type wallLaw_BC Constraint
// $type viscousWall_BC Constraint
// $type symmetry_BC Constraint
// $type impermeable_BC Constraint
// $type reflecting_BC Constraint
// $type interface_BC Constraint
// $type turboInterface_BC Constraint
// $type fixedMass_BC Constraint
// $type fixedMassOutflow_BC Constraint
// $type isentropicInflow_BC Constraint
// $type extrapolate_BC Constraint
// $type supersonicOutflow_BC Constraint
// $type outflow_BC Constraint
// $type supersonicInflow_BC Constraint
// $type farfield_BC Constraint
// $type inflow_BC Constraint

// $type plotFreq_BCoption Constraint
// $type qwall_BCoption Constraint
// $type prescribed_qwall_BCoption Constraint
// $type Twall_BCoption Constraint
// $type adiabatic_BCoption Constraint
// $type Uwall_BCoption Constraint
// $type stationary_BCoption Constraint
// $type absoluteFrame_BCoption Constraint
// $type normal_BCoption Constraint
// $type angVel_BCoption Constraint
// $type rotAxis_BCoption Constraint
// $type rotCenter_BCoption Constraint
// $type rotSpeed_BCoption Constraint
// $type flowDir_BCoption Constraint
// $type prescribed_BCoption Constraint
// $type T_BCoption Constraint
// $type T0_BCoption Constraint
// $type p_BCoption Constraint
// $type pMean_BCoption Constraint
// $type p0_BCoption Constraint
// $type rho_BCoption Constraint
// $type massFlux_BCoption Constraint
// $type mdot_BCoption Constraint
// $type swirlAngle_BCoption Constraint
// $type swirlCenter_BCoption Constraint
// $type swirlAxis_BCoption Constraint
// $type momentCenter_BCoption Constraint
// $type M_BCoption Constraint
// $type u_BCoption Constraint
// $type v_BCoption Constraint
// $type name_BCoption Constraint

// $type AllViscousBCs store<bool> 
// $type AllWallBCs store<bool> 

// $type componentMotion param<options_list> 
// $type componentHierarchy param<options_list> 
// $type componentGeometry param<options_list> 
// $type componentPriority param<options_list> 
// $type componentNodes_X store<bool> 

// $type movingMesh store<bool> 
// $type coriolis param<coriolis_options> 
// $type gauss_seidel_iter param<int> 
  





#line 25 "bodyforce.loci"

#line 1 "AD.lh"
//#############################################################################
//#
//# Copyright 2019, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

namespace flowPsi {
  struct actuatorDisk {
    real thrust ;
    real torque ;
    vect3d center ;
    vect3d axis ;
    real hub_radius ;
    real rotor_radius ;
    real thickness ;
    int type ;
  } ;
  typedef std::vector<actuatorDisk> actuatorDiskList ;
}

// $type actuatorDiskRegions param<Loci::options_list> 
// $type ActuatorDiskList blackbox<flowPsi::actuatorDiskList> 
// $type ActuatorID store<int> 
// $type ActuatorWeightThrust store<flowPsi::real> 
// $type ActuatorWeightTorque store<flowPsi::real> 
// $type ActuatorThrustFactor param<std::vector<flowPsi::real> > 
// $type ActuatorTorqueFactor param<std::vector<flowPsi::real> > 
// $type bodyForce store<flowPsi::vect3d> 
// $type bodyForce_f store<flowPsi::vect3d> 

#line 26 "bodyforce.loci"


using std::cerr ;
using std::endl ;

namespace flowPsi {


  namespace {class file_bodyforce000_1559584371m785 : public Loci::optional_rule {
#line 34 "bodyforce.loci"
    Loci::param<Loci::options_list>  L_actuatorDiskRegions_ ; 
#line 34 "bodyforce.loci"
public:
#line 34 "bodyforce.loci"
    file_bodyforce000_1559584371m785() {
#line 34 "bodyforce.loci"
       name_store("actuatorDiskRegions",L_actuatorDiskRegions_) ;
#line 34 "bodyforce.loci"
       output("actuatorDiskRegions") ;
#line 34 "bodyforce.loci"
    }
#line 34 "bodyforce.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 34 "bodyforce.loci"
Loci::register_rule<file_bodyforce000_1559584371m785> register_file_bodyforce000_1559584371m785 ;
#line 34 "bodyforce.loci"
}
#line 34 "bodyforce.loci"


  namespace {class file_bodyforce001_1559584371m786 : public Loci::singleton_rule {
#line 36 "bodyforce.loci"
    Loci::const_param<Loci::options_list>  L_actuatorDiskRegions_ ; 
#line 36 "bodyforce.loci"
    Loci::blackbox<flowPsi::actuatorDiskList>  L_ActuatorDiskList_ ; 
#line 36 "bodyforce.loci"
public:
#line 36 "bodyforce.loci"
    file_bodyforce001_1559584371m786() {
#line 36 "bodyforce.loci"
       name_store("actuatorDiskRegions",L_actuatorDiskRegions_) ;
#line 36 "bodyforce.loci"
       name_store("ActuatorDiskList",L_ActuatorDiskList_) ;
#line 36 "bodyforce.loci"
       input("actuatorDiskRegions") ;
#line 36 "bodyforce.loci"
       output("ActuatorDiskList") ;
#line 36 "bodyforce.loci"
    }
#line 36 "bodyforce.loci"
    void compute(const Loci::sequence &seq) { 
    std::vector<actuatorDisk> adl ;
    const options_list &ol = (*L_actuatorDiskRegions_);

    Loci::options_list::option_namelist nl = ol.getOptionNameList() ;
    Loci::options_list::option_namelist::const_iterator ii ;
    for(ii=nl.begin();ii!=nl.end();++ii) {
      string n = *ii ;
      Loci::option_value_type vt = ol.getOptionValueType(n) ;
      Loci::option_values ov = ol.getOption(n) ;
      if(vt == Loci::FUNCTION) {
	string name ;
	ov.get_value(name) ;
	options_list::arg_list value_list ;
	ov.get_value(value_list) ;
	options_list of ;
	of.Input(value_list) ;

	actuatorDisk ad ;
	ad.type = 0 ;
	if(name == "uniform") {
	  ad.type = 0 ;
	} else if(name == "goldstein") {
	  ad.type = 1 ;
	} else {
	  cerr << "only uniform or goldstein actuator disks types supported"
	       << endl ;
	} 
	
        ad.axis= vect3d(1,0,0) ;
        if(of.optionExists("axis")){
	  of.getOptionUnits("axis","",ad.axis) ;
	  ad.axis = ad.axis/(norm(ad.axis)+1e-30) ;
        } else {
	  cerr << "warning, actuator disk missing definition of 'axis'"
	       << endl ;
	}
       
        ad.center = vect3d(0,0,0) ;
        if(of.optionExists("center")) {
	  of.getOptionUnits("center","m",ad.center) ;
        } else {
	  cerr << "warning actuator disk missing definition of 'center'"
	       << endl ;
	}

        ad.hub_radius = 0.0 ;
        if(of.optionExists("hub_radius")) {
	  of.getOptionUnits("hub_radius","m",ad.hub_radius) ;
        } else {
	  cerr << "warning actuator disk missing definition of 'hub_radius'"
	       << endl ;
	}
	
        ad.thrust = 0.0 ;
        if(of.optionExists("thrust")) {
	  of.getOptionUnits("thrust","N",ad.thrust) ;
        } else {
	  cerr << "warning actuator disk missing definition of 'thrust'"
	       << endl ;
	}
	

        ad.torque = 0.0 ;
        if(of.optionExists("torque")) {
	  of.getOptionUnits("torque","N*m",ad.torque) ;
        } else {
	  cerr << "warning actuator disk missing definition of 'torque'"
	       << endl ;
        }

        ad.rotor_radius = 1.0 ;
        if(of.optionExists("rotor_radius")) {
	  of.getOptionUnits("rotor_radius","m",ad.rotor_radius) ;
        } else {
	  cerr << "warning actuator disk missing definition of 'rotor_radius'"
	       << endl ;
        }
 
        ad.thickness = 0.1 ;
        if(of.optionExists("thickness")) {
	  of.getOptionUnits("thickness","m",ad.thickness) ;
        } else {
	  cerr << "warning actuator disk missing definition of 'thickness'"
	       << endl ;
        }
	
	adl.push_back(ad) ;
      } else {
	cerr << "error in actuator disk specification" << endl ;
      }
    }
    (*L_ActuatorDiskList_)= adl ;
  }} ;
#line 129 "bodyforce.loci"
Loci::register_rule<file_bodyforce001_1559584371m786> register_file_bodyforce001_1559584371m786 ;
#line 129 "bodyforce.loci"
}
#line 129 "bodyforce.loci"


  namespace {class file_bodyforce002_1559584371m786 : public Loci::pointwise_rule {
#line 132 "bodyforce.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 132 "bodyforce.loci"
    Loci::const_blackbox<flowPsi::actuatorDiskList>  L_ActuatorDiskList_ ; 
#line 132 "bodyforce.loci"
    Loci::store<int>  L_ActuatorID_ ; 
#line 132 "bodyforce.loci"
    Loci::store<flowPsi::real>  L_ActuatorWeightThrust_ ; 
#line 132 "bodyforce.loci"
    Loci::store<flowPsi::real>  L_ActuatorWeightTorque_ ; 
#line 132 "bodyforce.loci"
public:
#line 132 "bodyforce.loci"
    file_bodyforce002_1559584371m786() {
#line 132 "bodyforce.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 132 "bodyforce.loci"
       name_store("ActuatorDiskList",L_ActuatorDiskList_) ;
#line 132 "bodyforce.loci"
       name_store("ActuatorID",L_ActuatorID_) ;
#line 132 "bodyforce.loci"
       name_store("ActuatorWeightThrust",L_ActuatorWeightThrust_) ;
#line 132 "bodyforce.loci"
       name_store("ActuatorWeightTorque",L_ActuatorWeightTorque_) ;
#line 132 "bodyforce.loci"
       input("cellcenter,ActuatorDiskList") ;
#line 132 "bodyforce.loci"
       output("ActuatorID") ;
#line 132 "bodyforce.loci"
       output("ActuatorWeightThrust") ;
#line 132 "bodyforce.loci"
       output("ActuatorWeightTorque") ;
#line 132 "bodyforce.loci"
    }
#line 132 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 133 "bodyforce.loci"
    int id = -1 ;
    for (size_t i =0;i <L_ActuatorDiskList_[_e_].size ();++i ) {
      vect3d pt = L_cellcenter_[_e_]- L_ActuatorDiskList_[_e_][i ].center ;
      real dist = dot (pt ,L_ActuatorDiskList_[_e_][i ].axis ) ;
      if (dist < 0 || dist > L_ActuatorDiskList_[_e_][i ].thickness )
	continue ;
      real r = norm (pt -dist *L_ActuatorDiskList_[_e_][i ].axis ) ;
      if ( r < L_ActuatorDiskList_[_e_][i ].hub_radius ||
	  r > L_ActuatorDiskList_[_e_][i ].rotor_radius )
	continue ;
      if (id != -1) {
	cerr << "overlapping actuator disks!" << endl ;
      }
      id = i ;
      L_ActuatorWeightThrust_[_e_]= 1 ;
      L_ActuatorWeightTorque_[_e_]= 1 ;
      if (L_ActuatorDiskList_[_e_][i ].type == 1) {
	real rp = r /L_ActuatorDiskList_[_e_][i ].rotor_radius ;
	real rhp = L_ActuatorDiskList_[_e_][i ].hub_radius /L_ActuatorDiskList_[_e_][i ].rotor_radius ;
	real rst = (rp -rhp )/(1.-rhp ) ;
	L_ActuatorWeightThrust_[_e_]= rst *sqrt (1.-rst ) ;
	L_ActuatorWeightTorque_[_e_]= L_ActuatorWeightThrust_[_e_]/(rst *(1-rhp )+rhp ) ;
      }
    }
    L_ActuatorID_[_e_]= id ;
  }    void compute(const Loci::sequence &seq) { 
#line 158 "bodyforce.loci"
      do_loop(seq,this) ;
#line 158 "bodyforce.loci"
    }
#line 158 "bodyforce.loci"
} ;
#line 158 "bodyforce.loci"
Loci::register_rule<file_bodyforce002_1559584371m786> register_file_bodyforce002_1559584371m786 ;
#line 158 "bodyforce.loci"
}
#line 158 "bodyforce.loci"


  template <class T> struct vectorSumOp {
    void operator()(T &lhs, const T &rhs) {
      const int sz = min(lhs.size(),rhs.size()) ;
      for(int i=0;i<sz;++i) {
	lhs[i] += rhs[i] ;
      }
    }
  } ;

  namespace {class file_bodyforce003_1559584371m787 : public Loci::unit_rule {
#line 169 "bodyforce.loci"
    Loci::const_blackbox<flowPsi::actuatorDiskList>  L_ActuatorDiskList_ ; 
#line 169 "bodyforce.loci"
    Loci::param<std::vector<flowPsi::real> >  L_ActuatorThrustFactor_ ; 
#line 169 "bodyforce.loci"
public:
#line 169 "bodyforce.loci"
    file_bodyforce003_1559584371m787() {
#line 169 "bodyforce.loci"
       name_store("ActuatorDiskList",L_ActuatorDiskList_) ;
#line 169 "bodyforce.loci"
       name_store("ActuatorThrustFactor",L_ActuatorThrustFactor_) ;
#line 169 "bodyforce.loci"
       input("ActuatorDiskList") ;
#line 169 "bodyforce.loci"
       output("ActuatorThrustFactor") ;
#line 169 "bodyforce.loci"
    }
#line 169 "bodyforce.loci"
    void compute(const Loci::sequence &seq) { 
    std::vector<real> tmp((*L_ActuatorDiskList_).size(),0) ;
    (*L_ActuatorThrustFactor_).swap(tmp) ;
  }} ;
#line 172 "bodyforce.loci"
Loci::register_rule<file_bodyforce003_1559584371m787> register_file_bodyforce003_1559584371m787 ;
#line 172 "bodyforce.loci"
}
#line 172 "bodyforce.loci"


  namespace {class file_bodyforce004_1559584371m787 : public Loci::apply_rule< param<std::vector<flowPsi::real> > ,vectorSumOp<std::vector<flowPsi::real> > >  {
#line 174 "bodyforce.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 174 "bodyforce.loci"
    Loci::const_store<int>  L_ActuatorID_ ; 
#line 174 "bodyforce.loci"
    Loci::const_store<flowPsi::real>  L_ActuatorWeightThrust_ ; 
#line 174 "bodyforce.loci"
    Loci::param<std::vector<flowPsi::real> >  L_ActuatorThrustFactor_ ; 
#line 174 "bodyforce.loci"
public:
#line 174 "bodyforce.loci"
    file_bodyforce004_1559584371m787() {
#line 174 "bodyforce.loci"
       name_store("vol",L_vol_) ;
#line 174 "bodyforce.loci"
       name_store("ActuatorID",L_ActuatorID_) ;
#line 174 "bodyforce.loci"
       name_store("ActuatorWeightThrust",L_ActuatorWeightThrust_) ;
#line 174 "bodyforce.loci"
       name_store("ActuatorThrustFactor",L_ActuatorThrustFactor_) ;
#line 174 "bodyforce.loci"
       input("ActuatorWeightThrust,ActuatorID,vol") ;
#line 174 "bodyforce.loci"
       output("ActuatorThrustFactor") ;
#line 174 "bodyforce.loci"
    }
#line 174 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 175 "bodyforce.loci"
    int id = L_ActuatorID_[_e_];
    if (id >=0) {
      L_ActuatorThrustFactor_[_e_][id ]+= L_vol_[_e_]*L_ActuatorWeightThrust_[_e_];
    }
  }    void compute(const Loci::sequence &seq) { 
#line 179 "bodyforce.loci"
      do_loop(seq,this) ;
#line 179 "bodyforce.loci"
    }
#line 179 "bodyforce.loci"
} ;
#line 179 "bodyforce.loci"
Loci::register_rule<file_bodyforce004_1559584371m787> register_file_bodyforce004_1559584371m787 ;
#line 179 "bodyforce.loci"
}
#line 179 "bodyforce.loci"

  
  namespace {class file_bodyforce005_1559584371m788 : public Loci::unit_rule {
#line 181 "bodyforce.loci"
    Loci::const_blackbox<flowPsi::actuatorDiskList>  L_ActuatorDiskList_ ; 
#line 181 "bodyforce.loci"
    Loci::param<std::vector<flowPsi::real> >  L_ActuatorTorqueFactor_ ; 
#line 181 "bodyforce.loci"
public:
#line 181 "bodyforce.loci"
    file_bodyforce005_1559584371m788() {
#line 181 "bodyforce.loci"
       name_store("ActuatorDiskList",L_ActuatorDiskList_) ;
#line 181 "bodyforce.loci"
       name_store("ActuatorTorqueFactor",L_ActuatorTorqueFactor_) ;
#line 181 "bodyforce.loci"
       input("ActuatorDiskList") ;
#line 181 "bodyforce.loci"
       output("ActuatorTorqueFactor") ;
#line 181 "bodyforce.loci"
    }
#line 181 "bodyforce.loci"
    void compute(const Loci::sequence &seq) { 
    std::vector<real> tmp((*L_ActuatorDiskList_).size(),0) ;
    (*L_ActuatorTorqueFactor_).swap(tmp) ;
  }} ;
#line 184 "bodyforce.loci"
Loci::register_rule<file_bodyforce005_1559584371m788> register_file_bodyforce005_1559584371m788 ;
#line 184 "bodyforce.loci"
}
#line 184 "bodyforce.loci"
  

  namespace {class file_bodyforce006_1559584371m788 : public Loci::apply_rule< param<std::vector<flowPsi::real> > ,vectorSumOp<std::vector<flowPsi::real> > >  {
#line 186 "bodyforce.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 186 "bodyforce.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 186 "bodyforce.loci"
    Loci::const_blackbox<flowPsi::actuatorDiskList>  L_ActuatorDiskList_ ; 
#line 186 "bodyforce.loci"
    Loci::const_store<int>  L_ActuatorID_ ; 
#line 186 "bodyforce.loci"
    Loci::const_store<flowPsi::real>  L_ActuatorWeightTorque_ ; 
#line 186 "bodyforce.loci"
    Loci::param<std::vector<flowPsi::real> >  L_ActuatorTorqueFactor_ ; 
#line 186 "bodyforce.loci"
public:
#line 186 "bodyforce.loci"
    file_bodyforce006_1559584371m788() {
#line 186 "bodyforce.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 186 "bodyforce.loci"
       name_store("vol",L_vol_) ;
#line 186 "bodyforce.loci"
       name_store("ActuatorDiskList",L_ActuatorDiskList_) ;
#line 186 "bodyforce.loci"
       name_store("ActuatorID",L_ActuatorID_) ;
#line 186 "bodyforce.loci"
       name_store("ActuatorWeightTorque",L_ActuatorWeightTorque_) ;
#line 186 "bodyforce.loci"
       name_store("ActuatorTorqueFactor",L_ActuatorTorqueFactor_) ;
#line 186 "bodyforce.loci"
       input("ActuatorWeightTorque,ActuatorID,vol,ActuatorDiskList,cellcenter") ;
#line 186 "bodyforce.loci"
       output("ActuatorTorqueFactor") ;
#line 186 "bodyforce.loci"
    }
#line 186 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 187 "bodyforce.loci"
    int id = L_ActuatorID_[_e_];
    if (id >=0) {
      vect3d pt = L_cellcenter_[_e_]- L_ActuatorDiskList_[_e_][id ].center ;
      real dist = dot (pt ,L_ActuatorDiskList_[_e_][id ].axis ) ;
      real r = norm (pt -dist *L_ActuatorDiskList_[_e_][id ].axis ) ;
      L_ActuatorTorqueFactor_[_e_][id ]+= r *L_vol_[_e_]*L_ActuatorWeightTorque_[_e_];
    }
  }    void compute(const Loci::sequence &seq) { 
#line 194 "bodyforce.loci"
      do_loop(seq,this) ;
#line 194 "bodyforce.loci"
    }
#line 194 "bodyforce.loci"
} ;
#line 194 "bodyforce.loci"
Loci::register_rule<file_bodyforce006_1559584371m788> register_file_bodyforce006_1559584371m788 ;
#line 194 "bodyforce.loci"
}
#line 194 "bodyforce.loci"


  // $type firstOrderCells store<char> 
  namespace {class file_bodyforce007_1559584371m788 : public Loci::apply_rule< store<char> ,Loci::Maximum<char> >  {
#line 198 "bodyforce.loci"
    Loci::const_Map L_cl_ ; 
#line 198 "bodyforce.loci"
    Loci::const_Map L_cr_ ; 
#line 198 "bodyforce.loci"
    Loci::const_store<int>  L_ActuatorID_ ; 
#line 198 "bodyforce.loci"
    Loci::store<char>  L_firstOrderCells_ ; 
#line 198 "bodyforce.loci"
public:
#line 198 "bodyforce.loci"
    file_bodyforce007_1559584371m788() {
#line 198 "bodyforce.loci"
       name_store("cl",L_cl_) ;
#line 198 "bodyforce.loci"
       name_store("cr",L_cr_) ;
#line 198 "bodyforce.loci"
       name_store("ActuatorID",L_ActuatorID_) ;
#line 198 "bodyforce.loci"
       name_store("firstOrderCells",L_firstOrderCells_) ;
#line 198 "bodyforce.loci"
       input("(cl,cr)->ActuatorID") ;
#line 198 "bodyforce.loci"
       output("(cl,cr)->firstOrderCells") ;
#line 198 "bodyforce.loci"
       constraint("(cl,cr)->geom_cells") ;
#line 198 "bodyforce.loci"
    }
#line 198 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 199 "bodyforce.loci"
    if (min (L_ActuatorID_[L_cl_[_e_]],L_ActuatorID_[L_cr_[_e_]]) < 0) {
      char tmp = 1 ;
      join (L_firstOrderCells_[L_cl_[_e_]],tmp ) ;
      join (L_firstOrderCells_[L_cr_[_e_]],tmp ) ;
    }
  }    void compute(const Loci::sequence &seq) { 
#line 204 "bodyforce.loci"
      do_loop(seq,this) ;
#line 204 "bodyforce.loci"
    }
#line 204 "bodyforce.loci"
} ;
#line 204 "bodyforce.loci"
Loci::register_rule<file_bodyforce007_1559584371m788> register_file_bodyforce007_1559584371m788 ;
#line 204 "bodyforce.loci"
}
#line 204 "bodyforce.loci"

  
  namespace {class file_bodyforce008_1559584371m788 : public Loci::pointwise_rule {
#line 208 "bodyforce.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 208 "bodyforce.loci"
    Loci::const_blackbox<flowPsi::actuatorDiskList>  L_ActuatorDiskList_ ; 
#line 208 "bodyforce.loci"
    Loci::const_store<int>  L_ActuatorID_ ; 
#line 208 "bodyforce.loci"
    Loci::const_store<flowPsi::real>  L_ActuatorWeightThrust_ ; 
#line 208 "bodyforce.loci"
    Loci::const_store<flowPsi::real>  L_ActuatorWeightTorque_ ; 
#line 208 "bodyforce.loci"
    Loci::const_param<std::vector<flowPsi::real> >  L_ActuatorThrustFactor_ ; 
#line 208 "bodyforce.loci"
    Loci::const_param<std::vector<flowPsi::real> >  L_ActuatorTorqueFactor_ ; 
#line 208 "bodyforce.loci"
    Loci::store<flowPsi::vect3d>  L_bodyForce_ ; 
#line 208 "bodyforce.loci"
public:
#line 208 "bodyforce.loci"
    file_bodyforce008_1559584371m788() {
#line 208 "bodyforce.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 208 "bodyforce.loci"
       name_store("ActuatorDiskList",L_ActuatorDiskList_) ;
#line 208 "bodyforce.loci"
       name_store("ActuatorID",L_ActuatorID_) ;
#line 208 "bodyforce.loci"
       name_store("ActuatorWeightThrust",L_ActuatorWeightThrust_) ;
#line 208 "bodyforce.loci"
       name_store("ActuatorWeightTorque",L_ActuatorWeightTorque_) ;
#line 208 "bodyforce.loci"
       name_store("ActuatorThrustFactor",L_ActuatorThrustFactor_) ;
#line 208 "bodyforce.loci"
       name_store("ActuatorTorqueFactor",L_ActuatorTorqueFactor_) ;
#line 208 "bodyforce.loci"
       name_store("bodyForce",L_bodyForce_) ;
#line 208 "bodyforce.loci"
       input("ActuatorID,cellcenter,ActuatorDiskList,ActuatorWeightThrust,ActuatorThrustFactor,ActuatorWeightTorque,ActuatorTorqueFactor") ;
#line 208 "bodyforce.loci"
       output("bodyForce") ;
#line 208 "bodyforce.loci"
    }
#line 208 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 209 "bodyforce.loci"
    vect3d bf (0,0,0) ;
    int id = L_ActuatorID_[_e_];
    if (id >=0) {
      const vect3d pt = L_cellcenter_[_e_]- L_ActuatorDiskList_[_e_][id ].center ;
      const vect3d n = L_ActuatorDiskList_[_e_][id ].axis ;
      const real dist = dot (pt ,n ) ;
      const vect3d radial = pt -dist *n ;
      vect3d nt = cross (n ,radial ) ;
      nt *= 1./max <real >(norm (nt ),1e-30) ;
      const real tr = L_ActuatorDiskList_[_e_][id ].thrust ;
      const real tq = L_ActuatorDiskList_[_e_][id ].torque ;
      bf = n *tr *L_ActuatorWeightThrust_[_e_]/L_ActuatorThrustFactor_[_e_][id ];
      bf += nt *tq *L_ActuatorWeightTorque_[_e_]/L_ActuatorTorqueFactor_[_e_][id ];
    }
    L_bodyForce_[_e_]= bf ;
  }    void compute(const Loci::sequence &seq) { 
#line 224 "bodyforce.loci"
      do_loop(seq,this) ;
#line 224 "bodyforce.loci"
    }
#line 224 "bodyforce.loci"
} ;
#line 224 "bodyforce.loci"
Loci::register_rule<file_bodyforce008_1559584371m788> register_file_bodyforce008_1559584371m788 ;
#line 224 "bodyforce.loci"
}
#line 224 "bodyforce.loci"


  
  namespace {class file_bodyforce009_1559584371m789 : public Loci::apply_rule< store<Loci::Array<flowPsi::real,5> > ,Loci::Summation<Loci::Array<flowPsi::real,5> > >  {
#line 228 "bodyforce.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 228 "bodyforce.loci"
    Loci::const_store<flowPsi::vect3d>  L_u_ ; 
#line 228 "bodyforce.loci"
    Loci::const_store<flowPsi::vect3d>  L_bodyForce_ ; 
#line 228 "bodyforce.loci"
    Loci::store<Loci::Array<flowPsi::real,5> >  L_src_ ; 
#line 228 "bodyforce.loci"
public:
#line 228 "bodyforce.loci"
    file_bodyforce009_1559584371m789() {
#line 228 "bodyforce.loci"
       name_store("vol",L_vol_) ;
#line 228 "bodyforce.loci"
       name_store("u",L_u_) ;
#line 228 "bodyforce.loci"
       name_store("src",L_src_) ;
#line 228 "bodyforce.loci"
       name_store("bodyForce",L_bodyForce_) ;
#line 228 "bodyforce.loci"
       input("bodyForce,vol,u") ;
#line 228 "bodyforce.loci"
       output("src") ;
#line 228 "bodyforce.loci"
       constraint("bodyForce,geom_cells") ;
#line 228 "bodyforce.loci"
    }
#line 228 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 229 "bodyforce.loci"
    const int mi = 1 ;
    const int ei = 4 ;
    L_src_[_e_][mi +0]+= L_bodyForce_[_e_].x *L_vol_[_e_];
    L_src_[_e_][mi +1]+= L_bodyForce_[_e_].y *L_vol_[_e_];
    L_src_[_e_][mi +2]+= L_bodyForce_[_e_].z *L_vol_[_e_];
    L_src_[_e_][ei ]+= dot (L_bodyForce_[_e_],L_u_[_e_])*L_vol_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 235 "bodyforce.loci"
      do_loop(seq,this) ;
#line 235 "bodyforce.loci"
    }
#line 235 "bodyforce.loci"
} ;
#line 235 "bodyforce.loci"
Loci::register_rule<file_bodyforce009_1559584371m789> register_file_bodyforce009_1559584371m789 ;
#line 235 "bodyforce.loci"
}
#line 235 "bodyforce.loci"


  // Add source term jacobian (energy dependency on u)
  // Since both primitive variable formulations are in terms of u the
  // jacobian is the same.
  namespace {class file_bodyforce010_1559584371m789 : public Loci::apply_rule< storeMat<flowPsi::real_fj> ,Loci::Summation<Mat<flowPsi::real_fj>  > >  {
#line 240 "bodyforce.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 240 "bodyforce.loci"
    Loci::const_store<flowPsi::vect3d>  L_bodyForce_ ; 
#line 240 "bodyforce.loci"
    Loci::storeMat<flowPsi::real_fj>  L_srcJ_ ; 
#line 240 "bodyforce.loci"
public:
#line 240 "bodyforce.loci"
    file_bodyforce010_1559584371m789() {
#line 240 "bodyforce.loci"
       name_store("vol",L_vol_) ;
#line 240 "bodyforce.loci"
       name_store("srcJ",L_srcJ_) ;
#line 240 "bodyforce.loci"
       name_store("bodyForce",L_bodyForce_) ;
#line 240 "bodyforce.loci"
       input("bodyForce,vol") ;
#line 240 "bodyforce.loci"
       output("srcJ") ;
#line 240 "bodyforce.loci"
    }
#line 240 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 241 "bodyforce.loci"
    const int mi = 1 ;
    const int ei = 4 ;
    L_srcJ_[_e_][ei ][mi +0] += L_vol_[_e_]*L_bodyForce_[_e_].x ;
    L_srcJ_[_e_][ei ][mi +1] += L_vol_[_e_]*L_bodyForce_[_e_].y ;
    L_srcJ_[_e_][ei ][mi +2] += L_vol_[_e_]*L_bodyForce_[_e_].z ;
  }    void compute(const Loci::sequence &seq) { 
#line 246 "bodyforce.loci"
      do_loop(seq,this) ;
#line 246 "bodyforce.loci"
    }
#line 246 "bodyforce.loci"
} ;
#line 246 "bodyforce.loci"
Loci::register_rule<file_bodyforce010_1559584371m789> register_file_bodyforce010_1559584371m789 ;
#line 246 "bodyforce.loci"
}
#line 246 "bodyforce.loci"


  namespace {class file_bodyforce011_1559584371m789 : public Loci::pointwise_rule {
#line 248 "bodyforce.loci"
    Loci::store<flowPsi::vect3d>  L_bodyForce_f_ ; 
#line 248 "bodyforce.loci"
public:
#line 248 "bodyforce.loci"
    file_bodyforce011_1559584371m789() {
#line 248 "bodyforce.loci"
       name_store("bodyForce_f",L_bodyForce_f_) ;
#line 248 "bodyforce.loci"
       output("bodyForce_f") ;
#line 248 "bodyforce.loci"
       constraint("ci->vol") ;
#line 248 "bodyforce.loci"
    }
#line 248 "bodyforce.loci"
    void calculate(Entity _e_) { 
#line 249 "bodyforce.loci"
    L_bodyForce_f_[_e_]= vect3d (0,0,0) ;
  }    void compute(const Loci::sequence &seq) { 
#line 250 "bodyforce.loci"
      do_loop(seq,this) ;
#line 250 "bodyforce.loci"
    }
#line 250 "bodyforce.loci"
} ;
#line 250 "bodyforce.loci"
Loci::register_rule<file_bodyforce011_1559584371m789> register_file_bodyforce011_1559584371m789 ;
#line 250 "bodyforce.loci"
}
#line 250 "bodyforce.loci"


  
  OUTPUT_VECTOR("cell2node_v3d(bodyForce)",bodyForce) ;
}

