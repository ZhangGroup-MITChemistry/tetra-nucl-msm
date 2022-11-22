# the example provided by Xingcheng

# steps to render:

In command line, run command "vmd -e surface_render_forBin.vmd"

Then in vmd command line, run command "source set_color.tcl", then run command "source set_user.tcl"

Next there are something we need to fix manually in vmd. Go to vmd-Graphics-Representations, in the final two representations, we need to manually fix Trajectory-Color Scale Data Range. For the final two representations, according to surface_render_forBin.vmd, parameters for "scaleminmanx" tell us the last two representations should have color scales 0-2086 and 3903-5989, respectively. Manually fix them. 

To write the plot, render with Tachyon. 

