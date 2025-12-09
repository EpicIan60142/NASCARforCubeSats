/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Justin Kottinger */

#pragma once

#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/spaces/SO2StateSpace.h>
#include <cmath>

namespace ob = ompl::base;

ob::StateSpacePtr createBoundedCubesatStateSpace(const std::pair<double, double> x_limits, const std::pair<double, double> y_limits, const std::pair<double, double> z_limits)
{
    // Compound state space
    ob::StateSpacePtr space = std::make_shared<ob::CompoundStateSpace>();
    //space->as<ob::CompoundStateSpace>()->addSubspace(ob::StateSpacePtr(new ob::RealVectorStateSpace(6)), 1.0);

    // Position space
    auto posSpace = std::make_shared<ob::RealVectorStateSpace>(3);
    ob::RealVectorBounds posBounds(3);
    posBounds.setLow(0, x_limits.first); //  x lower bound
    posBounds.setHigh(0, x_limits.second); // x upper bound
    posBounds.setLow(1, y_limits.first);  // y lower bound
    posBounds.setHigh(1, y_limits.second); // y upper bound
    posBounds.setLow(2, z_limits.first);  // z lower bound
    posBounds.setHigh(2, z_limits.second); // z upper bound
    posSpace->setBounds(posBounds);

    // Velocity space
    auto velSpace = std::make_shared<ob::RealVectorStateSpace>(3);
    ob::RealVectorBounds velBounds(3);
    velBounds.setLow(-10);
    velBounds.setHigh(10);
    velSpace->setBounds(velBounds);

    // Calculate weights
    double maxPosRange = std::max({
        x_limits.second - x_limits.first,
        y_limits.second - y_limits.first,
        z_limits.second - z_limits.first
    });
    double maxVelRange = 20;

    double posWeight = 1;
    double velWeight = maxPosRange/maxVelRange;

    // Add spaces to compound space with weights
    space->as<ob::CompoundStateSpace>()->addSubspace(posSpace, posWeight);
    space->as<ob::CompoundStateSpace>()->addSubspace(velSpace, velWeight);

    //space->as<ob::CompoundStateSpace>()->as<ob::RealVectorStateSpace>(0)->setBounds(bounds);

    return space;
}
