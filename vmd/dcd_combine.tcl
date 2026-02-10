# VMD script to combine all DCD files in a directory

# Load topology file (PDB)
mol new filter_4l.pdb  ;

# Get list of all DCD files in current directory
set dcd_files [glob *.dcd]
set dcd_files [lsort $dcd_files]  ;# Sort alphabetically

puts "Found [llength $dcd_files] DCD files:"
foreach file $dcd_files {
    puts "  $file"
}

# Load all DCD files
foreach dcd_file $dcd_files {
    puts "Loading $dcd_file..."
    mol addfile $dcd_file waitfor all
}

# Write combined trajectory
puts "Writing combined trajectory..."
animate write dcd combined_4l72.dcd beg 0 end -1 waitfor all

puts "Combined trajectory saved as combined_4l72.dcd"
quit