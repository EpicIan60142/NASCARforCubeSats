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

#include <boost/geometry/io/io.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <ompl/tools/config/SelfConfig.h>
#include <yaml-cpp/yaml.h>
#include "Eigen/Dense"
#include <memory>

typedef boost::geometry::model::d2::point_xy<double> point;
typedef boost::geometry::model::polygon<point> polygon;

// Ring defined by a center, normal vector, semi-major and -minor axes, and a DCM
struct Ring {
    // Members
    Eigen::Vector3d center_;
    Eigen::Vector3d normal_;
    double a_;
    double b_;
    Eigen::Matrix3d NR_;

    // Methods
    void printCenter() const
    {
        for (int i = 0; i < center_.size(); i++)
        {
            std::cout << "centerX: " << center_[0] << "centerY: " << center_[1] << "centerZ" << center_[2] << std::endl;
        }
    }
};


// A cubesat has a name, dynamics, max control, and start and goal regions
// Created as class to keep important variables safe
class Cubesat
{
public:
    // Constructor
    Cubesat(std::string name, std::string dyn, double uMax, Eigen::VectorXd s, Eigen::VectorXd g) {
        name_ = name;
        dynamics_ = dyn;
        uMax_ = uMax;
        start_ = s;
        goal_ = g;
    }

    // Getters
    std::string getName() const {return name_;}
    std::string getDynamics() const {return dynamics_;}
    double getUMax() const {return uMax_;}
    Eigen::VectorXd getStartLocation() const {return start_;}
    Eigen::VectorXd getGoalLocation() const {return goal_;}

    // Setters
    void setStartLocation(Eigen::VectorXd s) {start_ = s;}
    void setGoalLocation(Eigen::VectorXd g) {goal_ = g;}

private:
    std::string name_;
    std::string dynamics_;
    double uMax_;
    Eigen::VectorXd start_;
    Eigen::VectorXd goal_;
};

// World class holds all relevent data in the world that is used by OMPL
class World
{
public:
    // Constructors
        // Defgault
    World(){}

        // Copy
    World(const World& other) {
        // Copy basic members
        xBounds_ = other.xBounds_;
        yBounds_ = other.yBounds_;
        zBounds_ = other.zBounds_;
        n_ = other.n_;
        Rings_ = other.Rings_;  // Rings can be copied normally

        // Deep copy cubesats
        for (const auto& cubesat_ptr : other.Cubesats_) {
            auto new_cubesat = std::make_unique<Cubesat>(
                cubesat_ptr->getName(),
                cubesat_ptr->getDynamics(),
                cubesat_ptr->getUMax(),
                cubesat_ptr->getStartLocation(),
                cubesat_ptr->getGoalLocation()
            );
            Cubesats_.push_back(std::move(new_cubesat));
        }
    }

    // Copy assignment operator
    World& operator=(const World& other) {
        if (this != &other) {
            // Copy basic members
            xBounds_ = other.xBounds_;
            yBounds_ = other.yBounds_;
            zBounds_ = other.zBounds_;
            n_ = other.n_;
            Rings_ = other.Rings_;

            // Clear existing cubesats and deep copy new ones
            Cubesats_.clear();
            for (const auto& cubesat_ptr : other.Cubesats_) {
                auto new_cubesat = std::make_unique<Cubesat>(
                    cubesat_ptr->getName(),
                    cubesat_ptr->getDynamics(),
                    cubesat_ptr->getUMax(),
                    cubesat_ptr->getStartLocation(),
                    cubesat_ptr->getGoalLocation()
                );
                Cubesats_.push_back(std::move(new_cubesat));
            }
        }
        return *this;
    }


    // methods for dimensions
    void setWorldDimensions(std::pair<double, double> x, std::pair<double, double> y, std::pair<double, double> z)
    {
        xBounds_ = {x.first, x.second};
        yBounds_ = {y.first, y.second};
        zBounds_ = {z.first, z.second};
    }
    void setMeanMotion(double n){n_ = n;}
    std::vector<std::pair<double, double>> getWorldDimensions() const {return {xBounds_, yBounds_, zBounds_};}
    double getMeanMotion() const {return n_;}
    void printWorldDimensions(){OMPL_INFORM("Space Dimensions: [%0.2f, %0.2f; %0.2f, %0.2f; %0.2f, %0.2f]", xBounds_.first, xBounds_.second, yBounds_.first, yBounds_.second, zBounds_.first, zBounds_.second);}

    // methods for rings
    void addRing(Ring ring){Rings_.push_back(ring);}
    std::vector<Ring> getRings() const {return Rings_;}

    // methods for cubesats
    void addCubesat(std::unique_ptr<Cubesat> c){Cubesats_.push_back(std::move(c)); }
    std::vector<Cubesat*> getCubesats() const
    {
        std::vector<Cubesat*> result;
        for (const auto &c : Cubesats_)
            result.push_back(c.get());
        return result;
    }

    // printing methods for usability
    void printRings()
    {
        OMPL_INFORM("%d Rings (x, y; z): ", Rings_.size());
        for (Ring r: Rings_)
        {
            OMPL_INFORM("   - Ring: [%0.2f, %0.2f, %0.2f]", r.center_[0], r.center_[1], r.center_[2]);
        }
    }
    void printCubesats()
    {
        OMPL_INFORM("%d Cubesats: ", Cubesats_.size());
        for (const auto &cPtr : Cubesats_)
        {
            auto c = cPtr.get();
            OMPL_INFORM("   - Name: %s", c->getName().c_str());
            OMPL_INFORM("     Dynamics: %s", c->getDynamics().c_str());
            OMPL_INFORM("     uMax: %0.3f", c->getUMax());
            OMPL_INFORM("     Start: [%0.2f, %0.2f, %0.2f]", c->getStartLocation()[0], c->getStartLocation()[1], c->getStartLocation()[2]);
            OMPL_INFORM("     Goal: [%0.2f, %0.2f, %0.2f]", c->getGoalLocation()[0], c->getGoalLocation()[1], c->getGoalLocation()[2]);
        }
    }
    void printWorld()
    {
        printWorldDimensions();
        printRings();
        printCubesats();
    }
private:
    std::pair<double, double> xBounds_;
    std::pair<double, double> yBounds_;
    std::pair<double, double> zBounds_;
    double n_; // Course origin mean motion
    std::vector<Ring> Rings_;
    std::vector<std::unique_ptr<Cubesat>> Cubesats_;
};

// function that parses YAML file to world object
inline World* yaml2world(std::string file, bool verbose)
{
    YAML::Node config;
    World *w = new World();
    try
    {
        OMPL_INFORM("    Path to Problem File: %s", file.c_str());
        config = YAML::LoadFile(file);
        std::cout << "" << std::endl;
        OMPL_INFORM("    File loaded successfully. Parsing...");
    }
    catch (const std::exception& e) 
    {
        OMPL_ERROR("    Invalid file path. Aborting prematurely to avoid critical error.");
        exit(1);
    }
    try
    {        
        // grab dimensions from problem definition
        const auto& dims = config["Course"]["Dimensions"];
        const double x_min = dims[0].as<double>();
        const double x_max = dims[1].as<double>();
        const double y_min = dims[2].as<double>();
        const double y_max = dims[3].as<double>();
        const double z_min = dims[4].as<double>();
        const double z_max = dims[5].as<double>();
        w->setWorldDimensions({x_min, x_max}, {y_min, y_max}, {z_min, z_max});

        const double n = config["Course"]["MeanMotion"].as<double>();
        w->setMeanMotion(n);

        // set rings
        const auto& ring = config["Course"]["Rings"];
        for (int i=0; i < ring.size(); i++)
        {
            // Look at current ring
            std::string name = "ring" + std::to_string(i);

            // Pull out ring center
            Eigen::Vector3d center;
            center << ring[name][0].as<double>(), ring[name][1].as<double>(), ring[name][2].as<double>();

            // Pull out ring normal
            Eigen::Vector3d normal;
            normal << ring[name][3].as<double>(), ring[name][4].as<double>(), ring[name][5].as<double>();

            // Pull out semi-major and -minor axes
            const double a = ring[name][6].as<double>();
            const double b = ring[name][7].as<double>();

            // Pull out DCM to Inertial from Ring frame
            Eigen::Matrix3d NR;
            NR << ring[name][8].as<double>(), ring[name][9].as<double>(), ring[name][10].as<double>(),
                  ring[name][11].as<double>(), ring[name][12].as<double>(), ring[name][13].as<double>(),
                  ring[name][14].as<double>(), ring[name][15].as<double>(), ring[name][16].as<double>();

            // Add ring to the world
            Ring r = {center, normal, a, b, NR};
            w->addRing(r);
        }

        // Define cubesats
        const auto& cubesats = config["Cubesats"];
        for (int i = 0; i < cubesats.size(); i++)
        {
            std::string name = "sat" + std::to_string(i);
            Eigen::VectorXd start(6);
            start << cubesats[name]["Start"][0].as<double>(), cubesats[name]["Start"][1].as<double>(), cubesats[name]["Start"][2].as<double>(),
                     cubesats[name]["Start"][3].as<double>(), cubesats[name]["Start"][4].as<double>(), cubesats[name]["Start"][5].as<double>();

            Eigen::VectorXd goal(6);
            goal << w->getRings()[1].center_[0], w->getRings()[1].center_[1], w->getRings()[1].center_[2],
                    w->getRings()[1].normal_[0], w->getRings()[1].normal_[1], w->getRings()[1].normal_[2];
            auto c = std::make_unique<Cubesat>(cubesats[name]["Name"].as<std::string>(), cubesats[name]["Model"].as<std::string>(), cubesats[name]["uMax"].as<double>(),
                                                  start, goal);
            w->addCubesat(std::move(c));
        }
        OMPL_INFORM("    Parsing Complete.");
        std::cout << "" << std::endl;
        if (verbose)
            w->printWorld();

    }
    catch (const std::exception& e) 
    {
        OMPL_ERROR("    Error During Parsing. Aborting prematurely to avoid critical error.");
        exit(1);
    }
    return w;
}
