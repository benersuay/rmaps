# Bener Suay, Feb 2014
#
# benersuay@wpi.edu
#

## OPENRAVE ##
import openravepy
if not __openravepy_build_doc__:
    from openravepy import *
    from numpy import *

from openravepy.misc import OpenRAVEGlobalArguments

## ROBOT PLACEMENET ##
from Reachability import *

## MATH ##
from random import *

## SYSTEM - FILE OPS ##
import sys
import os
from datetime import datetime
import time
import commands

## Constraint Based Manipulation ##
from rodrigues import *
from TransformMatrix import *
from str2num import *
from TSR import *

env = Environment()
env.SetViewer('qtcoin')

print "Press Enter to generate the reachability map..."
sys.stdin.readline()

robot = env.ReadRobotURI('robots/barrettwam.robot.xml')
env.Add(robot)

rm = ReachabilityMap(robot,robot.GetManipulators()[0].GetName())

rm.r = 0
rm.g = 0
rm.b = 1

rm.xmin = 0.6
rm.xmax = 0.7

rm.ymin = 0.1
rm.ymax = 0.2

rm.zmin = 0.1
rm.zmax = 0.2

rm.generate(env)

print "Done. Press Enter to save the reachability map."
sys.stdin.readline()

# Save the reachability map
print "Saving the reachability map"
rm.name = "barrettwam_arm"
rm.save()

# Clean the object, just in case
del rm

print "Done. Press Enter to exit..."
sys.stdin.readline()

env.Destroy()
RaveDestroy()

