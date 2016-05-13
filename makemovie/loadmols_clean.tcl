#!/usr/bin/tclsh
#!/usr/local/bin/vmd
# vmd_select.tcl
#Name:  animatepdbs 
#Synopsis:
#  A Tcl script to load a consecutively numbered 
#  sequence of PDB files into VMD for animation purposes.
#Version: 1.0
#Uses VMD Version: 1.1
#Parameters:
# count - should be equal to count2 - not actually useless!!!
#  count2 - first index number for output .tga file - should be 0 for first section - distance from 0 in increment of 1 for others e.g. 11
#  start  - "frame number" of first cube file in sequence
#  end    - "frame number" of last cube file in sequence
#  step   - increment of input file index number
#  fileformat - a Tcl format string which describes the filename/numbering
#           used.
#
#Examples:
#  To load a sequence of PDB files named 0.pdb 1.pdb 2.pdb 3.pdb 
#  one would call this proc with:  animatepdbs 0 3 "%d.pdb"
#
#  To load a sequence of PDB files named foo0000.pdb foo0001.pdb foo0002.pdb
#  one would call this proc with:  animatepdbs 0 2 "foo%04d.pdb"
#
# Author: John Stone $lt;johns@ks.uiuc.edu&gt;

# Here begins my script
puts "Opening files"  


puts "Setting args" 
  
  set count 0
  set count2 0
  set start 0
  set end 10
  set step 1
  set fileformat "cat_tio2_%04d.cube"  
  
# Example values for 12 fs to 200 fs in steps of 20 * 0.1 = 2 fs 
# set count2 101
# set start 120
# set end 2000
# set step 20
# set fileformat "MK2COOH_TiO2ana101H2O_bridge_LUMO+0+1_2_xfs_%d.cube"

  for {set i $start} {$i <= $end} {incr i $step} {
    set filename [format $fileformat [expr $i]]
    set fileout [format "%04d" $count2].tga
    incr count     
    incr count2

puts $filename 

mol load cube $filename  

#Rotate an axis according to convenience: what viewpoint do you want to have?
# rotate x by = rotate around red
# rotate y by = rotate around green
# rotate z by = rotate around blue

 #rotate x by 90
 rotate y by 90
 rotate z by 270

#scale to 0.1

#Clone representation 0
::CloneRep::clone_reps 0 all   

#It turns all the molecules off
mol off all   

#Runs only molecule that we are seeing
mol on $count 

puts "NEW MOL"
#Generates a .tga file
render snapshot $fileout
mol delete $count 
#}   

