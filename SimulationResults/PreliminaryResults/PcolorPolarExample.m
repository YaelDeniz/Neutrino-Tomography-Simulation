addpath VisToolBox
clear all
close all

R = linspace(0.6, 1,100);
theta = linspace(0,100,100);
Z = linspace(0,10,100)'*linspace(0,10,100);
figure
polarPcolor(R,theta,Z)