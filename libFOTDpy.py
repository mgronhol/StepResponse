#!/usr/bin/env python3


import numpy as np
import math

EPS = 1e-3

def compute_stable_gain( Y, u ):
    Y = np.array( Y )
    u = np.array( u )

    start_idx = np.where( abs(u - u[0] ) < EPS )

    u0 = np.mean( u[start_idx] )
    Y0 = np.mean( Y[start_idx] )

    N = len(Y)

    lastdecile = int(0.95*N) - 2

    Ylast = Y[lastdecile:]
    ulast = u[lastdecile:]

    Y1 = np.mean( Ylast )
    u1 = np.mean( ulast )


    return (Y1 - Y0) / (u1 - u0)



def compute_gain_lag_tau( T, Y, u ):

    T = np.array( T )
    Y = np.array( Y )
    u = np.array( u )



    gain = compute_stable_gain( Y, u )

    idx = np.where( abs(u - u[0] ) > EPS )
    
    start_idx = idx[0][0]

    Y0 = np.mean( Y[:start_idx] )
    Tt = T[start_idx:] - T[start_idx]

    dY = Y - Y0
    dY = dY[start_idx:]
    N = len(dY)
    #print( "N", N)
    
    lastdecile = int(0.95*N) - 2
    Ylast = dY[lastdecile:]
    Y1 = np.mean( Ylast )

    T0 = T[start_idx]
    T1 = T[-1]

    #print( "Y0", Y0, "Y1", Y1)

    
    Amax = Y1 * ( T1 - T0 )
    #print( "T1-T0", T1-T0)

    A2 = np.trapz( dY, x = Tt )
    A2 = Amax - A2
    #print( "A2", A2)
    #print( "Amax", Amax)

    LT = A2 / Y1
    #print( "LT", LT)



    lt_idx = np.where( Tt < LT )

    A1 = np.trapz( dY[lt_idx], Tt[lt_idx] )

    #print( "lt_idx", lt_idx)
    

    T = A1/Y1 * math.exp(1)

    L = LT - T

    return gain, L, T


def get_PI_params_AstromHagglund( gain, L, T ):
    Kp = (0.63 * T) / ( gain * L)
    Ti = 3.2 * L
    return Kp, Ti



def get_P_params_Chien20_setpoint( gain, L, T ):
    a = gain * L / T
    Kp = 0.7 / a

    return Kp
    
def get_PI_params_Chien20_setpoint( gain, L, T ):
    Kp = (0.6 * T) / ( gain * L )
    Ti = T
    return Kp, Ti

def get_PID_params_Chien20_setpoint( gain, L, T ):
    a = gain * L / T
    Kp = 0.95 * a
    Ti = 1.4 * T
    Td = 0.47 * L
    return Kp, Ti, Td



def get_P_params_Chien0_setpoint( gain, L, T ):
    a = gain * L / T
    Kp = 0.3 / a

    return Kp

def get_PI_params_Chien0_setpoint( gain, L, T ):
    Kp = (0.35 * T) / ( gain * L )
    Ti = 1.17 * T
    return Kp, Ti    

def get_PID_params_Chien0_setpoint( gain, L, T ):
    a = gain * L / T
    Kp = 0.6 / a
    Ti = T
    Td = 0.5*L
    return Kp, Ti, Td

def get_P_params_Chien0_disturbance( gain, L, T ):
    a = gain * L / T
    Kp = 0.3 / a

    return Kp


def get_PI_params_Chien0_disturbance( gain, L, T ):
    a = gain * L / T
    Kp = 0.6 / a
    Ti = 4 * L
    return Kp, Ti    

def get_PID_params_Chien0_disturbance( gain, L, T ):
    a = gain * L / T
    Kp = 0.95 / a
    Ti = 2.4 * L
    Td = 0.42 * L
    return Kp, Ti, Td    


def get_P_params_Chien20_disturbance( gain, L, T ):
    a = gain * L / T
    Kp = 0.7 / a
    
    return Kp

def get_PI_params_Chien20_disturbance( gain, L, T ):
    a = gain * L / T
    Kp = 0.7 / a
    Ti = 2.3 * L
    return Kp, Ti    

def get_PID_params_Chien20_disturbance( gain, L, T ):
    a = gain * L / T
    Kp = 1.2 / a
    Ti = 2 * L
    Td = 0.42 * L
    return Kp, Ti, Td    


def _Zhuang_Atherton_compute_optimal_params2( gain, L, T, a1, b1, a2, b2 ):
    Kp = a1 / gain * math.pow( (L/T), b1 )
    Ti = T / (  a2 + b2 * (L/T) )
    return Kp, Ti

def _Zhuang_Atherton_compute_optimal_params3( gain, L, T, a1, b1, a2, b2, a3, b3 ):
    Kp = a1 / gain * math.pow( (L/T), b1 )
    Ti = T / (  a2 + b2 * (L/T) )
    Td = a3 * T * math.pow( (L/T), b3 )
    return Kp, Ti, Td


def get_PI_params_ISE_setpoint( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 0.980, -0.892, 0.690, -0.155)
    else:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.072, -0.560, 0.648, -0.144)

def get_PI_params_ISTE_setpoint( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 0.712, -0.921, 0.986, -0.247)
    else:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 0.786, -0.559, 0.883, -0.158)


def get_PI_params_IST2E_setpoint( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 0.569, -0.951, 1.023, -0.179)
    else:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 0.628, -0.583, 1.007, -0.167)


def get_PID_params_ISE_setpoint( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.048, -0.897, 1.195, -0.368, 0.489, 0.888)
    else:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.154, -0.567, 1.047, -0.220, 0.490, 0.708 )

def get_PID_params_ISTE_setpoint( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.042, -0.897, 0.987, -0.238, 0.385, 0.906)
    else:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.142, -0.579, 0.919, -0.172, 0.348, 0.839 )


def get_PID_params_IST2E_setpoint( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 0.968, -0.904, 0.977, -0.253, 0.316, 0.892 )
    else:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.061, -0.583, 0.892, -0.165, 0.315, 0.832 )




def get_PI_params_ISE_disturbance( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.279, -0.945, 0.535, 0.586 )
    else:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.346, -0.675, 0.552, 0.438 )


def get_PI_params_ISTE_disturbance( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.015, -0.957, 0.667, 0.552 )
    else:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.065, -0.673, 0.687, 0.427 )


def get_PI_params_IST2E_disturbance( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.021, -0.953, 0.629, 0.546 )
    else:
        return _Zhuang_Atherton_compute_optimal_params2( gain, L, T, 1.076, -0.648, 0.650, 0.442 )




def get_PID_params_ISE_disturbance( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.473, -0.970, 1.115, 0.753, 0.550, 0.948 )
    else:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.524, -0.735, 1.130, 0.641, 0.552, 0.933 )


def get_PID_params_ISTE_disturbance( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.468, -0.970, 0.942, 0.725, 0.443, 0.939 )
    else:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.515, -0.730, 0.957, 0.598, 0.444, 0.847 )


def get_PID_params_IST2E_disturbance( gain, L, T ):
    r = L/T
    if r <= 1.1:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.531, -0.960, 0.971, 0.746, 0.413, 0.933 )
    else:
        return _Zhuang_Atherton_compute_optimal_params3( gain, L, T, 1.592, -0.705, 0.957, 0.597, 0.414, 0.850 )



def get_normalized_dead_time( L, T ):
    return L / ( L + T)



def get_regulation_type( L, T ):
    q = T/L

    if q > 20:
        return "Bang-Bang"

    if 10 < q <= 20:
        return "P"
    
    if 5 < q <= 10:
        return "PI"
    
    if 0.5 < q <= 5:
        return "PID"
    
    if q <= 0.5:
        return "NA"



