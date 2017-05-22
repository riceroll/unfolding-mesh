============================================================
[Overview]

This software uses the following libraries:

OpenGL&GLUT (Mesa)
CGAL
Boost Graph library (BGL)
GNU Scientific Library (GSL)

Before compilig the program, please make sure that the above
software has been installed in your computational
enviroment.
============================================================
[Compilation]

On linux
% make linux

On Mac 
% make mac
============================================================
[Usage]
Let us try to unfold the "bunny-128.off" model.

1.) Start the program:

% ./unfold bunny-128.off

Press the left button for menu.  In the "Mesh" window, you
can change the viewpoint and distance by dragging the right
button and the middle button (vertically), respectively.


2.) Select the menu

[Split]->[Spanning Trees]
(or press 's')

for decomposing the mesh into a set of patches with a small
number of faces.

3.) Select the menu

[Stitch]->[GA-based stitch]
(or press 'o')

for composing a single unfolded patch (or a few unfolded
patches).

4.) If you still have a few unfolded patches, select the menu

[Stitch]->[Greedy Remerge]

several times until you can merge them into one single patch.

5.) Select the menu

[Render]->[Papercraft]
(or press '6')

for rendering the boundary of the mesh unfolding with arrows
of different colors.

[Render]->[Pattern]
(or press '9')
for rendering edges in black and white.

Note that solid and broken lines correspond to mountain and
valley folds.

6.) You can control the degree of freedom in the initial set
of small patches by selecting the menu

[ST Deviation]->[??%]

where ??% represents the ?? percents of the dual edges where
split in this first decomposition stage.

7.) For finding multiple patches using MST, select the menu

[MST] 
(or press 'm') 

after Step 1.)




============================================================
[Reference]

Shigeo Takahashi Hsiang-Yun Wu, Seow Hui Saw, Chun-Cheng
Lin, and Hsu-Chun Yen: "Optimized Topological Surgery for
Unfolding 3D Meshes," Computer Graphics Forum (Proceedings
of Pacific Graphics 2011), Vol. 30, No. 7, pp. 2077-2086,
2011.

http://www.tak-lab.org/research/unfolding/
