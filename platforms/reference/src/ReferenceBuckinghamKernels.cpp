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

#include "ReferenceBuckinghamKernels.h"
#include "BuckinghamForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include <iostream>

using namespace BuckinghamPlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (Vec3*) data->periodicBoxVectors;
}

ReferenceCalcBuckinghamForceKernel::~ReferenceCalcBuckinghamForceKernel() {
    if (neighborList) {
        delete neighborList;
    } 
}

void ReferenceCalcBuckinghamForceKernel::initialize(const System& system, const BuckinghamForce& force) {
    // Initialize particles parameters.
    
    numParticles = system.getNumParticles();
    a.resize(numParticles);
    b.resize(numParticles);
    c6.resize(numParticles);
    c8.resize(numParticles);
    c10.resize(numParticles);
    gamma.resize(numParticles);
    exclusions.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
        force.getParticleParameters(i, a[i], b[i], c6[i], c8[i], c10[i], gamma[i]);
    usePeriodic = force.usesPeriodicBoundaryConditions();
    cutoff = force.getCutoffDistance();
    neighborList = new OpenMM::NeighborList();
}

double ReferenceCalcBuckinghamForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& forces = extractForces(context);
    RealVec* boxVectors = extractBoxVectors(context);

    cout << "asdasd" << pos[1] << endl;
    cout << "asdasd" << boxVectors[1] << endl;
    computeNeighborListVoxelHash(*neighborList, numParticles, pos, exclusions, boxVectors, usePeriodic, cutoff, 0.0);

    double energy = 0.0;

    if (usePeriodic) {
        double minAllowedSize = 1.999999*cutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize) {
            throw OpenMMException("The periodic box size has decreased to less than twice the cutoff.");
        }

        OpenMM::NeighborList nl = *neighborList;
        cout << "PERIODIC -> " << endl;
        for (unsigned int ii = 0; ii < nl.size(); ii++) {
            OpenMM::AtomPair pair = nl[ii];
            int siteI = pair.first;
            int siteJ = pair.second;

            // Calculate combined parameters
            double combinedA = sqrt(a[siteI] * a[siteJ]);
            double combinedB = sqrt(b[siteI] * b[siteJ]);
            double combinedC6 = sqrt(c6[siteI] * c6[siteJ]);
            double combinedC8 = sqrt(c8[siteI] * c8[siteJ]);
            double combinedC10 = sqrt(c10[siteI] * c10[siteJ]);
            double d = gamma[siteI] + gamma[siteJ];

            // Calculate distance
            double deltaR[ReferenceForce::LastDeltaRIndex];
            ReferenceForce::getDeltaRPeriodic(pos[siteI], pos[siteJ], boxVectors, deltaR);

            double r2 = deltaR[ReferenceForce::R2Index];
            double r = deltaR[ReferenceForce::RIndex];
            double invR = 1.0f / r;
            double invR2 = invR * invR;
            double invR3 = invR2 * invR;
            double invR4 = invR3 * invR;
            double invR5 = invR4 * invR;
            double invR6 = invR5 * invR;
            double invR7 = invR6 * invR;
            double invR8 = invR7 * invR;
            double invR9 = invR8 * invR;
            double invR10 = invR9 * invR;

            // Calculate energy and forces
            double tempEnergy, tempForce, dEdR;

            double d2 = d * d * 0.5f;
            double d3 = d2 * d * 0.3333333333f;
            double d4 = d3 * d * 0.25f;
            double d5 = d4 * d * 0.2f;
            double d6 = d5 * d * 0.1666666667f;
            double d7 = d6 * d * 0.1428571429f;
            double d8 = d7 * d * 0.125f;
            double d9 = d8 * d * 0.1111111111f;
            double d10 = d9 * d * 0.1f;

            double mdr = -d * r;
            double expTerm = exp(mdr);

            double c6Deriv = 6.0f * invR6 + expTerm * (
                invR6 * (mdr - 6.0f) +
                d * invR5  * (mdr - 5.0f) +
                d2 * invR4 * (mdr - 4.0f) +
                d3 * invR3 * (mdr - 3.0f) +
                d4 * invR2 * (mdr - 2.0f) +
                d5 * invR * (mdr - 1.0f) +
                d6 * mdr
            );

            double c8Deriv = 8.0f * invR8 + expTerm * (
                invR8 * (mdr - 8.0f) +
                d * invR7 * (mdr - 7.0f) +
                d2 * invR6 * (mdr - 6.0f) +
                d3 * invR5 * (mdr - 5.0f) +
                d4 * invR4 * (mdr - 4.0f) +
                d5 * invR3 * (mdr - 3.0f) +
                d6 * invR2 * (mdr - 2.0f) +
                d7 * invR * (mdr - 1.0f) +
                d8 * mdr
            );

            double c10Deriv = 10.0f * invR10 + expTerm * (
                invR10 * (mdr - 10.0f) +
                d * invR9 * (mdr - 9.0f) +
                d2 * invR8 * (mdr - 8.0f) +
                d3 * invR7 * (mdr - 7.0f) +
                d4 * invR6 * (mdr - 6.0f) +
                d5 * invR5 * (mdr - 5.0f) +
                d6 * invR4 * (mdr - 4.0f) +
                d7 * invR3 * (mdr - 3.0f) +
                d8 * invR2 * (mdr - 2.0f) +
                d9 * invR * (mdr - 1.0f) +
                d10 * mdr
            );

            double c6E = invR6 - expTerm * (
                invR6 +
                d * invR5 +
                d2 * invR4 +
                d3 * invR3 +
                d4 * invR2 +
                d5 * invR +
                d6
            );

            double c8E = invR8 - expTerm * (
                invR8 +
                d * invR7 +
                d2 * invR6 +
                d3 * invR5 +
                d4 * invR4 +
                d5 * invR3 +
                d6 * invR2 +
                d7 * invR +
                d8
            );

            double c10E = invR10 - expTerm * (
                invR10 +
                d * invR9 +
                d2 * invR8 +
                d3 * invR7 +
                d4 * invR6 +
                d5 * invR5 +
                d6 * invR4 +
                d7 * invR3 +
                d8 * invR2 +
                d9 * invR +
                d10
            );

            double buckinghamExp = -1.0f * combinedB * r;
            double buckinghamRepulsion = combinedA * exp(buckinghamExp);

            tempForce = -buckinghamExp * buckinghamRepulsion - c6Deriv * combinedC6 - c8Deriv * combinedC8 - c10Deriv * combinedC10;
            tempEnergy = buckinghamRepulsion - c6E * combinedC6 - c8E * combinedC8 -c10E * combinedC10;
            dEdR = tempForce * invR * invR;

            energy += tempEnergy;

            for (int kk = 0; kk < 3; kk++) { // x, y, z
                double _force = dEdR * deltaR[kk];
                forces[siteI][kk] += _force;
                forces[siteJ][kk] -= _force;
            }
        }
    }

    return energy;
}

void ReferenceCalcBuckinghamForceKernel::copyParametersToContext(ContextImpl& context, const BuckinghamForce& force) {
    if (force.getNumParticles() != numParticles)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    for (int i = 0; i < force.getNumParticles(); i++) {
        force.getParticleParameters(i, a[i], b[i], c6[i], c8[i], c10[i], gamma[i]);
    }
}
