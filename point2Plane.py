#---------------------------------------------------------------------------------------------------
# Name: point2Plane.py
# Author: Yolanda
# Instruction: To calculate the distance of a point to a plane, which is defined by 3 other points,
#    user should input the coordinates of 3 points in the plane into (x1,y1,z1)(x2,y2,z2)(x3,y3,z3),
#    also input the coordinates of the point out of the plane into (x0,y0,z0). This program will
#    print the distance.
# Application: Measure the distance of a atom to a benzene ring.
#---------------------------------------------------------------------------------------------------



import numpy as np
from scipy import linalg

# 3 points to define a plane. For example, 3 atoms in a benzene ring.
x1, y1, z1 = 0.421, 9.340, 10.017
x2, y2, z2 = -0.042, 8.673, 8.866
x3, y3, z3 = 0.785, 8.316, 7.853

# Train the equation of the plane. Equation: Ax + By + Cz + 1 = 0
A = np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])
b = np.array([[-1],[-1],[-1]])
pr = list(linalg.inv(A).dot(b).flat) + [1] # pr == [A, B, C, 1]

# The point out of the plane.
x0, y0, z0 = 2.691, 11.980, 9.187

# Calculte the distance of the point to the plane.
d = np.abs(sum([p*x for p,x in zip(pr, [x0,y0,z0,1])])) / np.sqrt(sum([a**2 for a in pr[:3]]))

print "Distance: ", d
