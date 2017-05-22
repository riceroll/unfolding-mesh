# meshed_unfolder
This code is modified from [Optimized Topological Surgery for Unfolding 3D Meshes](http://xueshu.baidu.com/s?wd=paperuri%3A%289a2c9ddde2f4b48f4cedfb1fa13793a7%29&filter=sc_long_sign&tn=SE_xueshusource_2kduw22v&sc_vurl=http%3A%2F%2Fonlinelibrary.wiley.com%2Fdoi%2F10.1111%2Fj.1467-8659.2011.02053.x%2Fpdf&ie=utf-8&sc_us=15904426517141857991)

## Requirement
- OpenGL&GLUT (Mesa)
- CGAL
- Boost Graph library (BGL)
- GNU Scientific Library (GSL)

## Compilation

On linux
```bash
make linux
```
On Mac 
```bash
make mac
```
## Quick Start

To unfold the "bunny-128.off" model.

### I. Start the program:
```bash
./unfold bunny-128.off
```
Press the left button for menu.  In the "Mesh" window, you
can change the viewpoint and distance by dragging the right
button and the middle button (vertically), respectively.

### II. Select the menu
[Split]->[Spanning Trees]
(or press 's')

for decomposing the mesh into a set of patches with a small
number of faces.

### III. Select the menu

[Stitch]->[GA-based stitch]
(or press 'o')

for composing a single unfolded patch (or a few unfolded
patches).

### IV. If you still have a few unfolded patches, select the menu

[Stitch]->[Greedy Remerge]

several times until you can merge them into one single patch.

## Reference

Shigeo Takahashi Hsiang-Yun Wu, Seow Hui Saw, Chun-Cheng
Lin, and Hsu-Chun Yen: ["Optimized Topological Surgery for
Unfolding 3D Meshes,"](http://xueshu.baidu.com/s?wd=paperuri%3A%289a2c9ddde2f4b48f4cedfb1fa13793a7%29&filter=sc_long_sign&tn=SE_xueshusource_2kduw22v&sc_vurl=http%3A%2F%2Fonlinelibrary.wiley.com%2Fdoi%2F10.1111%2Fj.1467-8659.2011.02053.x%2Fpdf&ie=utf-8&sc_us=15904426517141857991) Computer Graphics Forum (Proceedings
of Pacific Graphics 2011), Vol. 30, No. 7, pp. 2077-2086,
2011.

