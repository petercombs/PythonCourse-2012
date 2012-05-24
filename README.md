PythonCourse-2012
=================

Code for the QB3 Summer Python Course, run in July/August 2012.

Hey everyone!  Here's where I'm storing all the code I've written for the
class. If you go click on the "Commits" link, you can see my progress as I
wrote things, including the inevitable missteps.  As Einstein once said, "if we
knew what we were doing, it wouldn't be called research".


Github has options for making all of this stuff Private, so the lab you're
competing with on a paper won't be able to lift your code and figure out what
analyses you're doing, but on the other hand, your code *is* your method, so if
you're doing anything non-trivial, there's definitely a case to be made for
making it available after publication, at the very least. See, for example,
http://sciencecodemanifesto.org/

CalcFalloff.py
-------------

This is code for finding the average(?) expression across all genes as a
function of position.  We'd like to see the upstream area (which should be
relatively low), the gene itself (relatively high, and nearly constant across
the gene), and then the downstream, which ought to fall off quickly in the
no-drug case, and fall off much more slowly in the +drug case.

PlotUtils.py
------------

Throughout the second week, we're making a bunch of different plots.  Some of
these will be easy, others not so much. For the more complicated graphs, we'd
like to have a dedicated module for making them, so we can call them quickly
and easily from the interactive prompt, as well as making everything more
clear.


