clearvars
close all
clc

vertex0 = [-4.704, -0.96, 0];
vertex1 = [4.704, -0.98613, 0];
vertex2 = [4.704, -0.98613, 1.32];

A = vertex1 - vertex0;
B = vertex2 - vertex0;


facetNormal = cross(A, B);
facetNormalNorm = normalize(facetNormal, 'norm')