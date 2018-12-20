# FeedbackLoopDetection
Tools in R and Python to find feedback loops (FBLs) in a system of ordinary differential equations (ODEs)

A biological system can be mathematically represented as a system of ordinary differential equations (ODEs). The Jacobian matrix of a system of ODEs is the matrix of the partial derivatives. The Jacobian matrix can be represented as a graph where each species is a node and an entry unequal to zero in the matrix is an edge in the graph. 
With directed path detection in the graph it is possible to determine if a species is affecting itself either directly (self-loop) or via other species (loop). The tools can calculate self-loops and loops, determine the length of each loop and show if these are positive or negative.

In this repository are two codes provided, one in R and one in Python, to analyse ODEs for FBLs. The files can be read in in R with the command *source("FBL_detection.r")* and in Python with the command *exec(open("FBL_detection.py").read())*. The usage and functionalities are shown in the documentations.

The tools were developed in the course of a research internship at the Wolf Lab of the Max-Delbrück-Centrum for Molecular Medicine (MDC) in Berlin.
I want to thank Katharina Baum for the good supervision, support and helpful tips and the group leader Jana Wolf for the opportunity of an internship in her group.
