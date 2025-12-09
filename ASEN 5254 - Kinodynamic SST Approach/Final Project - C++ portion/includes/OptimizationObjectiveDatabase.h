//
// Created by ianmf on 12/2/25.
//

#pragma once

#include <ompl/base/OptimizationObjective.h>
#include <ompl/base/objectives/StateCostIntegralObjective.h>
#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>

#include "World.h"
#include "Eigen/Dense"



class CenterObjective : public ob::StateCostIntegralObjective
{
public:
    CenterObjective(const ob::SpaceInformationPtr &si, const Ring ring) :
                    ob::StateCostIntegralObjective(si, true),
                    ring_(ring)
    {}

    ob::Cost stateCost(const ob::State *state) const override
    {
        const double *q = state->as<ob::CompoundStateSpace::StateType>()->as<ob::RealVectorStateSpace::StateType>(0)->values;

        const double x = q[0];
        const double y = q[1];
        const double z = q[2];

        Eigen::Vector3d sat_pos(x, y, z);
        Eigen::Vector3d ring_center = ring_.center_;

        ob::Cost cost = ob::Cost((sat_pos - ring_.center_).norm());

        return(ob::Cost((sat_pos - ring_.center_).norm()));
    }
private:
    Ring ring_;
};

ob::OptimizationObjectivePtr getTimeObjective(const ob::SpaceInformationPtr& si)
{
    return ob::OptimizationObjectivePtr(new ob::PathLengthOptimizationObjective(si));
}

ob::OptimizationObjectivePtr getCenterObjective(const ob::SpaceInformationPtr& si, Ring ring)
{
    return ob::OptimizationObjectivePtr(new CenterObjective(si, ring));
}

ob::OptimizationObjectivePtr getMinimizeCostObjective(const ob::SpaceInformationPtr& si, const Ring ring)
{
    ob::OptimizationObjectivePtr timeObj = getTimeObjective(si);
    ob::OptimizationObjectivePtr centerObj = getCenterObjective(si, ring);

    auto opt = std::make_shared<ob::MultiOptimizationObjective>(si);
    opt->addObjective(timeObj, 1.0);
    opt->addObjective(centerObj, 1.0);
    return ob::OptimizationObjectivePtr(opt);
}
