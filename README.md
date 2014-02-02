Reachability Maps
=================

This is a small OpenRAVE library that tests different kinematic configuration for a manipulator and saves successful (collision free) end effector positions wrt base of the manipulator for future reference. This might be handy when trying to solve the robot base placement problem for instance. The concept is based on Franziska Zacharias' capability maps paper (Zacharias et al., The Capability Map: A Tool for Reasoning in Mobile Manipulation).

Dependency
----------

The scripts in this repo depends on OpenRAVE and Constrained Manipulation Planning Suite scripts for matrix operations (or CoMPS: http://sourceforge.net/p/comps/wiki/Installation/).

Run
---

To run the sample rmap generator for the barrett arm in OpenRAVE

`python barrettwam_generate_and_save_reachability_map.py`

The script above will generate a very tiny portion of the full reachability map defined by the limits written in the script, and it will save the map and the map parameters in pickled format. 

To load a previously generated map run

`python barrettwam_load_reachability_map.py`



Todo
----

- Change orientation sampling to uniformly distributed points on a sphere.
- Make Reachability.py an independent OpenRAVE Plug-in.

See it in action on YouTube
---------------------------

<a href="http://www.youtube.com/watch?feature=player_embedded&v=3cpZh7JDMf0
" target="_blank"><img src="http://img.youtube.com/vi/3cpZh7JDMf0/0.jpg" 
alt="Speed test." width="240" height="180" border="10" /></a>
