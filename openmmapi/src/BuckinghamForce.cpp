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

#include "BuckinghamForce.h"
#include "internal/BuckinghamForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"

using namespace BuckinghamPlugin;
using namespace OpenMM;
using namespace std;

BuckinghamForce::BuckinghamForce() {
}

int BuckinghamForce::addParticle(double a, double b, double c6, double c8, double c10, double gamma) {
    particles.push_back(ParticleInfo(a, b, c6, c8, c10, gamma));
    return particles.size()-1;
}

void BuckinghamForce::getParticleParameters(int index, double& a, double& b, double& c6, double& c8, double& c10, double& gamma) const {
    ASSERT_VALID_INDEX(index, particles);
    a = particles[index].a;
    b = particles[index].b;
    c6 = particles[index].c6;
    c8 = particles[index].c8;
    c10 = particles[index].c10;
    gamma = particles[index].gamma;
}

void BuckinghamForce::setParticleParameters(int index, double a, double b, double c6, double c8, double c10, double gamma) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].a = a;
    particles[index].b = b;
    particles[index].c6 = c6;
    particles[index].c8 = c8;
    particles[index].c10 = c10;
    particles[index].gamma = gamma;
}

double BuckinghamForce::getCutoffDistance() const {
    return cutoffDistance;
}

void BuckinghamForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

ForceImpl* BuckinghamForce::createImpl() const {
    return new BuckinghamForceImpl(*this);
}

void BuckinghamForce::updateParametersInContext(Context& context) {
    dynamic_cast<BuckinghamForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
