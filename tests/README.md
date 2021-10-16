## Simple testing

A PSF file and a DCD file from a short simulation of ~100 water droplet with one fluoride ion is given.

Run the following command:

`Water_order -f testfile.dcd -f testfile.psf -t <task> -c OH2 --rmax 20 --bin-width 0.2`

Replace <task> with OTO, d5, Sk and rhoV.

Then compare the result with the test results I found. (A quick way to do this would be run a checksum on the CSV files. Ensure that the line endings are the same, becuase it might have a different checksum due to the different line endings.)

If the results do not match, then open an issue.
