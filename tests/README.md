## Simple testing

A PSF file and a DCD file from a short simulation of ~100 water droplet with one fluoride ion is given.

Run the following command:
```
Water_order -f testfile.dcd -f testfile.psf -t <task> -c OH2 --rmax 20 --bin-width 0.2
```

Replace `<task>` with `OTO`, `d5`, `Sk` and `rhoV`.

Then compare the result with the test results from my benchmark. Keep in mind that I used the Visual C/C++ compiled version on Windows to make these test results. There might be minor differences in the floating point precision if you use a different compiler and different system, but it is unlikely because only 6 decimal points are printed in the output. I have tested with Debian 10 with gcc/g++, and the numbers are the exact same there. (So a quick way to compare would be to run checksum on the CSV files, although keep in mind the differences in line endings: git usually handles that by default)

If the results do not match, then please open an issue.
