# Bener Suay, April 2013
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

robot = env.ReadRobotURI('robots/barrettwam.robot.xml')
env.Add(robot)

rm = ReachabilityMap(robot,robot.GetManipulators()[0].GetName())


print "Press Enter to load barrettwam's reachability map."
sys.stdin.readline()

# Note to self: Add attribute "name" in params class
print "Loading..."
rm.load("barrettwam_arm")
rm.name = "barrettwam_arm"

print "Loaded the reachability map for barrettwam_arm."

# Implement a rm.info() method
# Until then...
print "number of reachability spheres:"
print len(rm.map)

print "Press Enter to show the map..."
sys.stdin.readline()
rm.show(env)

print "Done. Press Enter to exit..."
sys.stdin.readline()

# Clean the object
del rm

# Clean all
env.Destroy()
RaveDestroy()

