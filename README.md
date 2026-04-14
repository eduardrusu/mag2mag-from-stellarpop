Fork of [LensPop](https://github.com/tcollett/LensPop), but removing all functionality unrelated to `mag2mag.py`.
The original code is now properly linted, and I have added error bar estimation for the function `ab_filter_magnitude()`.

From the original README:
- Our SED code gives correct magnitudes if the SED is in units of ergs/s/cm/cm/Angstrom, and our distance class gives distances in Mpc/h, so we will need to use the conversion Mpc = 3.08568e24cm and assume a value for h.
- The code is open access, but please email thomas.collett@port.ac.uk to say that you are using it. Please cite Collett 2015, if you make use of these codes in your work.

Example usage:

	python mag2mag.py -T Sbc_cww -m1 24 -m1_err 0.01 -f1 V_Johnson -z1 1.00 -f2 U_Johnson -z2 0.00 -errbar -plot
