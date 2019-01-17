/* -------------------------------------------------------------------------- *
 *                                OpenMMExample                                 *
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

#include "BuckinghamForceProxy.h"
#include "BuckinghamForce.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace BuckinghamPlugin;
using namespace OpenMM;
using namespace std;

BuckinghamForceProxy::BuckinghamForceProxy() : SerializationProxy("BuckinghamForce") {
}

void BuckinghamForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const BuckinghamForce& force = *reinterpret_cast<const BuckinghamForce*>(object);
    SerializationNode& particles = node.createChildNode("Particles");

    for (int i = 0; i < force.getNumParticles(); i++) {
        double a, b, c6, c8, c10, gamma;
        force.getParticleParameters(i, a, b, c6, c8, c10, gamma);
        particles.createChildNode("Atom").setDoubleProperty("a", a).setDoubleProperty("b", b).setDoubleProperty("c6", c6)
            .setDoubleProperty("c8", c8).setDoubleProperty("c10", c10).setDoubleProperty("gamma", gamma);
    }
}

void* BuckinghamForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    BuckinghamForce* force = new BuckinghamForce();
    try {
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            force->addParticle(particle.getDoubleProperty("a"), particle.getDoubleProperty("b"), particle.getDoubleProperty("c6"), particle.getDoubleProperty("c8"), particle.getDoubleProperty("c10"), particle.getDoubleProperty("gamma"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
