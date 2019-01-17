%module buckinghamplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "BuckinghamForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.
*/
%pythonappend BuckinghamPlugin::BuckinghamForce::getBondParameters(int index, double& c6, double& c8,
                                                 double& c10, double& a, double& b, double& gamma) const %{
    val[0]=unit.Quantity(val[0], unit.kilojoule_per_mole)
    val[1]=unit.Quantity(val[1], unit.kilojoule_per_mole*unit.nanometer**-6)
    val[2]=unit.Quantity(val[2], unit.kilojoule_per_mole*unit.nanometer**-8)
    val[3]=unit.Quantity(val[3], unit.kilojoule_per_mole*unit.nanometer**-10)
    val[5]=unit.Quantity(val[5], unit.kilojoule_per_mole)
    val[6]=unit.Quantity(val[6], unit.nanometer)
    val[7]=unit.Quantity(val[7], unit.nanometer**-1)
%}


namespace BuckinghamPlugin {

class BuckinghamForce : public OpenMM::Force {
public:
    BuckinghamForce();

    int getNumParticles() const;

    int addParticle(double a, double b, double c6, double c8, double c10, double gamma);

    void setParticleParameters(int index, double a, double b, double c6, double c8, double c10, double gamma);

    void updateParametersInContext(OpenMM::Context& context);

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply double& OUTPUT {double& c6};
    %apply double& OUTPUT {double& c8};
    %apply double& OUTPUT {double& C10};
    %apply double& OUTPUT {double& a};
    %apply double& OUTPUT {double& b};
    %apply double& OUTPUT {double& gamma};
    void getParticleParameters(int index, double& c6, double& c8, double& c10, double& a, double& b, double& gamma) const;
    %clear double& c6;
    %clear double& c8;
    %clear double& c10;
    %clear double& a;
    %clear double& b;
    %clear double& gamma;
};

}
