

from ete3 import Tree, TreeStyle

# When you see this error: This is due to ete3 is dependent on X11, 

# qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
# This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
# Available platform plugins are: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-egl, wayland-xcomposite-glx, webgl, xcb.

# Set these environmental variables to solve it!
# export QT_QPA_PLATFORM=offscreen
# export XDG_RUNTIME_DIR=/tmp/runtime-runner


t = Tree("/medstore/projects/P23-044/Intermediate/SCOPE-Batch1_HBVintegrationPipeline/VirusPipelineOut/barcode04/barcode04_Merged_NJ.newick")



ts=TreeStyle()
ts.show_leaf_name=True
ts.mode='c'


t.render("mytree.pdf", w=1600, units="px",tree_style=ts)



#print(t)
