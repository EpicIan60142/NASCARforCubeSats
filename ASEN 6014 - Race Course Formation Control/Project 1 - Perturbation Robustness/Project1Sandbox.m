%% ASEN 6014 Project 1 Sandbox
% By: Ian Faber, 09/27/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

addpath("../../ASEN6080/Utilities/Dynamics/")

%% Setup
pConst.mu = 398600.4415;
pConst.Ri = 6378;
pConst.J2 = 1.08264e-3;

scConst = 0;

oe.mu = pConst.mu;
oe.a = 7500;
oe.e = 0.1;
oe.i = deg2rad(35);
oe.RAAN = deg2rad(47);
oe.argPeri = deg2rad(13);
oe.f = deg2rad(45);

cart = convClassicOE2Cart(oe);

X0 = [cart.rVec; cart.vVec];

tspan = 0:10:86400;

opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

%% Compare 6080 J2 acceleration to 6014 J2 acceleration
    % 6014 approach
[t_6014, X_6014] = ode45(@(t,X)orbitEOM(t,X,pConst,scConst,true), tspan, X0, opt);

    % 6080 approach
[t_6080, X_6080] = ode45(@(t,X)orbitEOM_MuJ2(t,X,pConst.mu,pConst.J2,pConst.Ri), tspan, X0, opt);

    % plot differences
figure;
hold on; grid on;
title("J2 perturbed orbits")
plot3(X_6014(:,1), X_6014(:,2), X_6014(:,3), 'b-', X_6080(:,1), X_6080(:,2), X_6080(:,3), 'r-')
legend("6014", "6080")
view([30 35])

figure;
hold on; grid on;
title("Difference between J2 perturbed orbits")
plot3(X_6014(:,1) - X_6080(:,1), X_6014(:,2) - X_6080(:,2), X_6014(:,3) - X_6080(:,3))
view([30 35])

figure; tl = tiledlayout(3,1);
title(tl, "Difference between ASEN 6014 and 6080 J2 methods")
nexttile;
    hold on; grid on
    plot(tspan, X_6014(:,1) - X_6080(:,1), 'b-')
    title("\Deltax vs time")
nexttile;
    hold on; grid on
    plot(tspan, X_6014(:,2) - X_6080(:,2), 'r-')
    title("\Deltay vs time")
nexttile;
    hold on; grid on
    plot(tspan, X_6014(:,3) - X_6080(:,3), 'k-')
    title("\Deltaz vs time")




