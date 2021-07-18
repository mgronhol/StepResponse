# StepResponse
Compute control parameters from step response



## Example 

```
#!/usr/bin/env python3

import libFOTD
import sys

## Data must contain one step response measurement
## Must be long enough for the measured value to stabilize 

T = [] # time
Y = [] # measured value
u = [] # control output
    
with open( sys.argv[1], 'r' ) as handle:
  for line in handle:
    parts = line.split()
    T.append( float( parts[0] ) )
    Y.append( float( parts[1] ) )
    u.append( float( parts[2] ) )
        
   
gain, lag, tau = libFOTD.compute_gain_lag_tau( T, Y, u )
print( "Regulator type:", libFOTD.get_regulation_type( lag, tau ) )
   
# Get some PID gains
Kp, Ti, Td = libFOTD.get_PID_params_Chien0_setpoint( gain, lag, tau )
print( "Kp:", Kp )
print( "Ti:", Ti )
print( "Td:", Td )
```
