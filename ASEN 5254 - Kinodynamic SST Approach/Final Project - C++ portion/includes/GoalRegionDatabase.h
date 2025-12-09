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

#include <ompl/base/goals/GoalRegion.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <utility>

namespace ob = ompl::base;

class GoalRegionCubesat: public ob::GoalRegion
{
public:
    GoalRegionCubesat(const ob::SpaceInformationPtr &si, Eigen::VectorXd goal, Ring ring)
                        : ob::GoalRegion(si), goal_(goal), ring_(ring)
    {
        threshold_ = std::min(ring.a_, ring.b_);
    }
    
    double distanceGoal(const ob::State *st) const override
    {
        const double *sat_pos_OMPL = st->as<ob::CompoundStateSpace::StateType>()->as<ob::RealVectorStateSpace::StateType>(0)->values;
        Eigen::VectorXd sat_pos(3);
        sat_pos << sat_pos_OMPL[0], sat_pos_OMPL[1], sat_pos_OMPL[2];

        Eigen::VectorXd goal_pos(3);
        goal_pos << goal_[0], goal_[1], goal_[2];

        Eigen::VectorXd diff = sat_pos - goal_pos;
        return diff.norm();
    }

    bool isSatisfied(const ob::State *st) const override
    {
        const double *sat_vel_OMPL = st->as<ob::CompoundStateSpace::StateType>()->as<ob::RealVectorStateSpace::StateType>(1)->values;

        OMPL_INFORM("isSatisfied called!");

        Eigen::VectorXd heading(3);
        heading << sat_vel_OMPL[0], sat_vel_OMPL[1], sat_vel_OMPL[2];
        heading = heading/heading.norm();

        double angleDiff = std::acos(heading.dot(ring_.normal_)); // angle between heading and ring normal in radians

        return (distanceGoal(st) <= std::min(ring_.a_, ring_.b_)) && (angleDiff <= angleEpsilon_);
    }

private:
    Eigen::VectorXd goal_; // Goal state
    Ring ring_; // Goal ring
    double angleEpsilon_ = 89*M_PI/180; // Angle tolerance on final velocity heading - want to fly through the ring

};
