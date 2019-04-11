/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

/**
 * This tests the Reference implementation of BuckinghamForce.
 */

#include "BuckinghamForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace BuckinghamPlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerBuckinghamReferenceKernelFactories();

void testForce() {
    // Create a chain of particles connected by bonds.
    
    const double boxSize = 2.0;
    const int numParticles = 2;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(8.0);
        positions[i] = Vec3(0.5*i, 0.1*i, 0.1*i);
    }
    BuckinghamForce* force = new BuckinghamForce();
    force->setCutoffDistance(0.9);
    
    system.addForce(force);
    for (int i = 0; i < numParticles; i++)
        force->addParticle(8.0e+05, 3.8e+01, 0.0026399910, 0.0001937841, 0.0000166944, 35.8967);
    
    // Compute the forces and energy.

    VerletIntegrator integ(1.0);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);

    context.setPeriodicBoxVectors(Vec3(2*boxSize, 0, 0), Vec3(0, 2*boxSize, 0), Vec3(0, 0, 2*boxSize));

    State state = context.getState(State::Energy | State::Forces);
    
    // See if the energy is correct.
    
    double expectedEnergy = -0.180097;
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);
}

int main() {
    try {
        registerBuckinghamReferenceKernelFactories();
        testForce();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
