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

#include <ompl/control/ODESolver.h>
#include <ompl/base/spaces/SO2StateSpace.h>
#include <ompl/control/StatePropagator.h>
#include <ompl/base/Goal.h>
#include <ompl/control/spaces/RealVectorControlSpace.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


namespace ob = ompl::base;
namespace oc = ompl::control;

void CubesatCWHODE (const oc::ODESolver::StateType& q, const oc::Control* control, oc::ODESolver::StateType& qdot, double n)
{
    // q = x, y, z, v_x, v_y, v_z
    // u = [u_x, u_y, u_z]
    const double *u = control->as<oc::RealVectorControlSpace::ControlType>()->values;

    // state params
    const double x = q[0];
    const double y = q[1];
    const double z = q[2];
    const double v_x = q[3];
    const double v_y = q[4];
    const double v_z = q[5];

    // Zero out qdot
    qdot.resize (q.size(), 0);
 
    // vehicle model
    qdot[0] = v_x;
    qdot[1] = v_y;
    qdot[2] = v_z;
    qdot[3] = 2*n*v_y + 3*n*n*x + u[0];
    qdot[4] = -2*n*v_x + u[1];
    qdot[5] = -n*n*z + u[2];
}
