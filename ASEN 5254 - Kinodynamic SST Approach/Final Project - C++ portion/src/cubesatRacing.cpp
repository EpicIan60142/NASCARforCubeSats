//
//  Main source file for the ASEN 5254 Final Project
//
// Created by ianmf on 12/2/25.
//

#include <ompl/control/SimpleSetup.h>
#include <ompl/control/planners/rrt/RRT.h>
#include <ompl/control/planners/sst/SST.h>
#include "World.h"
#include "StateSpaceDatabase.h"
#include "ControlSpaceDatabase.h"
#include "StateValidityCheckerDatabase.h"
#include "StatePropagatorDatabase.h"
#include "GoalRegionDatabase.h"
#include "OptimizationObjectiveDatabase.h"
#include "PostProcessing.h"

namespace fs = std::filesystem;
namespace ob = ompl::base;
namespace oc = ompl::control;

// this function sets-up an ompl planning problem for an arbtrary number of agents
oc::SimpleSetupPtr controlSimpleSetUp(const World *w, const size_t ring_idx, Cubesat *c)
{
    // Create state and control spaces
    ob::StateSpacePtr space = createBoundedCubesatStateSpace(w->getWorldDimensions()[0], w->getWorldDimensions()[1], w->getWorldDimensions()[2]);
    oc::ControlSpacePtr cspace = createUniform3DRealVectorControlSpace(space, c);

    // Define a simple setup class from the state and control spaces
    auto ss = std::make_shared<oc::SimpleSetup>(cspace);

    // Set state validity checker
    ss->setStateValidityChecker(std::make_shared<isStateValid_3D>(ss->getSpaceInformation(), w, c));

    // Use the ODESolver to propagate the system
    auto odeFunction = [w](const oc::ODESolver::StateType &q, const oc::Control *control, oc::ODESolver::StateType &qdot)
    {
        CubesatCWHODE(q, control, qdot, w->getMeanMotion());
    };
    auto odeSolver(std::make_shared<oc::ODEBasicSolver<>>(ss->getSpaceInformation(), odeFunction));
    ss->setStatePropagator(oc::ODESolver::getStatePropagator(odeSolver));

    // Set up optimization objectives
    ss->setOptimizationObjective(getMinimizeCostObjective(ss->getSpaceInformation(), w->getRings()[ring_idx]));

    // Create start state
    ob::ScopedState<> start(space);
    start[0] = c->getStartLocation()[0];
    start[1] = c->getStartLocation()[1];
    start[2] = c->getStartLocation()[2];
    start[3] = c->getStartLocation()[3];
    start[4] = c->getStartLocation()[4];
    start[5] = c->getStartLocation()[5];

    // Create goal region
    ob::GoalPtr goal(new GoalRegionCubesat(ss->getSpaceInformation(), c->getGoalLocation(), w->getRings()[ring_idx]));

    // Set propagation step size and min/max number of steps
    ss->getSpaceInformation()->setPropagationStepSize(0.05);
    ss->getSpaceInformation()->setMinMaxControlDuration(5, 20);

    // Add start and goal to problem setup
    ss->setStartState(start);
    ss->setGoal(goal);

    OMPL_INFORM("    Successfully Setup the problem instance");
    return ss;
}

// Sequential planning function between rings
std::vector<oc::PathControl> planSequential(const std::string plannerString, const std::string problemFile, const int cubesatIdx)
{
    // Extract the overall world
    World *w = yaml2world("../problems/" + problemFile + ".yml", true);
    auto worldDims = w->getWorldDimensions();

    // Extract the rings and cubesat of interest
    auto rings = w->getRings();
    auto cubesat = w->getCubesats()[cubesatIdx];

    std::vector<oc::PathControl> segments;

    OMPL_INFORM("  Starting sequential planning through %zu rings", rings.size());

    // Current state starts at cubesat start location
    Eigen::VectorXd current_state = cubesat->getStartLocation();

    for (size_t ring_idx = 1; ring_idx < rings.size(); ring_idx++)
    {
        OMPL_INFORM("    Planning segment %zu: Ring %zu -> Ring %zu",
                    ring_idx, ring_idx-1, ring_idx);

        // Create a simplified world for this segment
        World segmentWorld = *w;  // Copy world

        double buffer = 100; // m, buffer in each direction in addition to existing coordinate differences

        std::vector<double> x_coords = {rings[ring_idx].center_[0], rings[ring_idx-1].center_[0], current_state[0]};
        std::vector<double> y_coords = {rings[ring_idx].center_[1], rings[ring_idx-1].center_[1], current_state[1]};
        std::vector<double> z_coords = {rings[ring_idx].center_[2], rings[ring_idx-1].center_[2], current_state[2]};

        auto minX = *std::min_element(x_coords.begin(), x_coords.end()) - buffer;
        auto maxX = *std::max_element(x_coords.begin(), x_coords.end()) + buffer;
        auto minY = *std::min_element(y_coords.begin(), y_coords.end()) - buffer;
        auto maxY = *std::max_element(y_coords.begin(), y_coords.end()) + buffer;
        auto minZ = *std::min_element(z_coords.begin(), z_coords.end()) - buffer;
        auto maxZ = *std::max_element(z_coords.begin(), z_coords.end()) + buffer;

        segmentWorld.setWorldDimensions({minX, maxX}, {minY, maxY}, {minZ, maxZ});

        // Set goal to current ring
        Eigen::VectorXd ringGoal(6);
        ringGoal << rings[ring_idx].center_[0], rings[ring_idx].center_[1], rings[ring_idx].center_[2],
                     rings[ring_idx].normal_[0], rings[ring_idx].normal_[1], rings[ring_idx].normal_[2];
        cubesat->setGoalLocation(ringGoal);

        OMPL_INFORM("    === SEGMENT %zu DEBUGGING ===", ring_idx);
        OMPL_INFORM("    Start state: [%.2f, %.2f, %.2f] vel:[%.2f, %.2f, %.2f]",
                    current_state[0], current_state[1], current_state[2],
                    current_state[3], current_state[4], current_state[5]);
        OMPL_INFORM("    Target Ring %zu center: [%.2f, %.2f, %.2f] normal: [%.2f, %.2f, %.2f]",
                    ring_idx, rings[ring_idx].center_[0], rings[ring_idx].center_[1], rings[ring_idx].center_[2],
                    rings[ring_idx].normal_[0], rings[ring_idx].normal_[1], rings[ring_idx].normal_[2]);
        OMPL_INFORM("    Ring %zu dimensions: a=%.2f, b=%.2f, threshold=%.2f",
                    ring_idx, rings[ring_idx].a_, rings[ring_idx].b_, std::min(rings[ring_idx].a_, rings[ring_idx].b_));

        // Set up planning problem for this segment
        oc::SimpleSetupPtr ss = controlSimpleSetUp(&segmentWorld, ring_idx, cubesat);

        // Configure planner
        if (plannerString == "SST")
        {
            auto planner = std::make_shared<oc::SST>(ss->getSpaceInformation());

            planner->setGoalBias(0.3);
            planner->setSelectionRadius(15.0);
            planner->setPruningRadius(8.0);

            ss->setPlanner(planner);

            OMPL_INFORM("    SST configured for segment %zu:", ring_idx);
            OMPL_INFORM("      Goal bias: %.2f, Selection: %.2f, Pruning: %.2f",
                        planner.get()->getGoalBias(), planner.get()->getSelectionRadius(), planner.get()->getPruningRadius());
        }
        else
        {
            auto planner = std::make_shared<oc::RRT>(ss->getSpaceInformation());

            planner->setGoalBias(0.2);

            ss->setPlanner(planner);

            OMPL_INFORM("    RRT configured for segment %zu:", ring_idx);
            OMPL_INFORM("      Goal bias: %.2f", planner.get()->getGoalBias());
        }

        ss->setup();

        // Solve this segment
        double segment_timeout = 60.0;
        bool solved = ss->solve(segment_timeout);

        if (solved)
        {
            auto segment_path = ss->getSolutionPath();
            segments.push_back(segment_path);

            OMPL_INFORM("    Segment %zu solved! Path time: %.2f sec",
                        ring_idx,
                        segment_path.length());

            // Update current_state to the final state of this segment
            auto final_state = segment_path.getState(segment_path.getStateCount() - 1);
            const double *pos_vals = final_state->as<ob::CompoundStateSpace::StateType>()
                                      ->as<ob::RealVectorStateSpace::StateType>(0)->values;
            const double *vel_vals = final_state->as<ob::CompoundStateSpace::StateType>()
                                      ->as<ob::RealVectorStateSpace::StateType>(1)->values;

            Eigen::Vector3d pos(pos_vals[0], pos_vals[1], pos_vals[2]);
            Eigen::Vector3d vel(vel_vals[0], vel_vals[1], vel_vals[2]);

            OMPL_INFORM("      Final distance to ring: %.2f m, final angle to ring normal: %.3f deg",
                        (pos - rings[ring_idx].center_).norm(), std::acos(vel.dot(rings[ring_idx].normal_)/vel.norm()) * 180.0 / M_PI);

            for (int i = 0; i < 3; i++)
            {
                current_state[i] = pos_vals[i];
            }

            for (int i = 0; i < 3; i++)
            {
                current_state[i + 3] = vel_vals[i];
            }

            OMPL_INFORM("    Updated start state for next segment: [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]",
                        current_state[0], current_state[1], current_state[2],
                        current_state[3], current_state[4], current_state[5]);

            // Update cubesat start to current state
            cubesat->setStartLocation(current_state);

        } else {
            OMPL_ERROR("    Failed to solve segment %zu", ring_idx);
            break;  // Stop if any segment fails
        }
    }

    return segments;

}

// Function to combine segment paths into a single solution
oc::PathControl concatenatePaths(const std::vector<oc::PathControl>& segments,
                                 const ob::SpaceInformationPtr& si) {

    if (segments.empty()) {
        throw std::runtime_error("    No segments to concatenate");
    }

    // Start with first segment
    oc::PathControl combined_path = segments[0];

    // Append remaining segments
    for (size_t i = 1; i < segments.size(); i++) {
        // Get states and controls from current segment
        for (size_t j = 1; j < segments[i].getStateCount(); j++) {  // Skip first state to avoid duplication
            combined_path.append(segments[i].getState(j),
                               segments[i].getControl(j-1),
                               segments[i].getControlDuration(j-1));
        }
    }

    OMPL_INFORM("    Combined path has %u states, total length: %.2f",
                combined_path.getStateCount(), combined_path.length());

    return combined_path;
}


// main planning function -- uses simple setup
void planControl(std::string planner_string, std::string problem_file)
{
    OMPL_INFORM("Starting sequential ring-by-ring planning");

    for (int i = 0; i < 4; i++)
    {
        // Create a dummy world for output formatting and debug logging
        World *w = yaml2world("../problems/" + problem_file + ".yml", false);

        // Select cubesat
        int cubesat_idx = i;

        OMPL_INFORM("  Planning for cubesat %s", w->getCubesats()[cubesat_idx]->getName().c_str());

        // Plan segments
        std::vector<oc::PathControl> segments = planSequential(planner_string, problem_file, cubesat_idx);

        if (segments.empty()) {
            OMPL_ERROR("    Sequential planning failed - no segments completed");
            return;
        }

        if (segments.size() < 2) {  // Adjust based on expected number of rings
            OMPL_WARN("      Only completed %zu segments - incomplete solution", segments.size());
        }

        oc::SimpleSetupPtr ss = controlSimpleSetUp(w, 1, w->getCubesats()[cubesat_idx]);

        // Combine segments into single path
        try {
            oc::PathControl combined_path = concatenatePaths(segments, ss->getSpaceInformation());

            // Set the combined path as the solution
            ss->getProblemDefinition()->clearSolutionPaths();
            ss->getProblemDefinition()->addSolutionPath(std::make_shared<oc::PathControl>(combined_path));

            OMPL_INFORM("    Sequential planning SUCCESS!");
            OMPL_INFORM("    Total segments: %zu", segments.size());
            OMPL_INFORM("    Combined path length: %.2f", combined_path.length());

            // Write solution
            write2sys(ss, w->getCubesats()[cubesat_idx], problem_file);

        } catch (const std::exception& e) {
            OMPL_ERROR("    Failed to concatenate paths: %s", e.what());
        }
    }
}

int main(int argc, char ** argv) {
    std::string plannerName = "SST";
    std::string problem = "RaceCourse";
    OMPL_INFORM("Running Race Course Problem with %s", plannerName.c_str());
    planControl(plannerName, problem);
}



