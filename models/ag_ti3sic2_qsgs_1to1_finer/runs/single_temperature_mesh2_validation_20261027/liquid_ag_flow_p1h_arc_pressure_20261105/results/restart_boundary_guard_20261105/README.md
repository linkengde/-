# Restart boundary-mode guard

A v3 closed-top checkpoint was restored with the open-top executable and otherwise matching Marangoni-off, uniform 1.8 MPa pressure and 0.02 ps source/pressure ramp settings. The solver exited with status 2 and reported `restart case does not match requested sigma/Marangoni/top-boundary mode` before time integration. This confirms cross-boundary restarts are rejected.
