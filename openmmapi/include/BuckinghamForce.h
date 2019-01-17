#ifndef OPENMM_BUCKINGHAMFORCE_H_
#define OPENMM_BUCKINGHAMFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Memorial University and the Authors.           *
 * Authors: Mohamad Mohebifar                                                 *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportBuckingham.h"

namespace BuckinghamPlugin {

class OPENMM_EXPORT_BUCKINGHAM BuckinghamForce : public OpenMM::Force {
public:
    /**
     * Create an BuckinghamForce.
     */
    BuckinghamForce();
    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Add a particle to the force.
     *
     * @param a         the A coefficient in the buckingham potential
     * @param b         the b coefficient in the buckingham potential
     * @param c6        the C6 coefficient
     * @param c8        the C8 coefficient
     * @param c10       the C10 coefficient
     * @param gamma     the Tang-Toennis damping coefficient
     * @return the index of the particle that was added
     */
    int addParticle(double a, double b, double c6, double c8, double c10, double gamma);
    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the particle for which to get parameters
     * @param a         the A coefficient in the buckingham potential
     * @param b         the b coefficient in the buckingham potential
     * @param c6        the C6 coefficient
     * @param c8        the C8 coefficient
     * @param c10       the C10 coefficient
     */
    void getParticleParameters(int index, double& a, double& b, double& c6, double& c8, double& c10, double& gamma) const;
    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the particle for which to set parameters
     * @param a         the A coefficient in the buckingham potential
     * @param b         the b coefficient in the buckingham potential
     * @param c6        the C6 coefficient
     * @param c8        the C8 coefficient
     * @param c10       the C10 coefficient
     * @param gamma     the Tang-Toennis damping coefficient
     */
    void setParticleParameters(int index, double a, double b, double c6, double c8, double c10, double gamma);
    /**
     * Get the cutoff distance (in nm) being used for buckingham interactions.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for buckingham interactions.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
        return true;
    }
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    double cutoffDistance;
    std::vector<ParticleInfo> particles;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class BuckinghamForce::ParticleInfo {
public:
    double a, b, c6, c8, c10, gamma;
    ParticleInfo() {
        a = b = c6 = c8 = c10 = gamma = 0.0;
    }
    ParticleInfo(double a, double b, double c6, double c8, double c10, double gamma) :
        a(a), b(b), c6(c6), c8(c8), c10(c10), gamma(gamma) {
    }
};

} // namespace BuckinghamPlugin

#endif /*OPENMM_BUCKINGHAMFORCE_H_*/
