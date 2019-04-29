#!/home/rsimons/miniconda2/envs/grizli/bin/python

def OH_R23(OH, use = 'M08', prt = False):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('R23 outside of Maiolino+ 08 calibration range...')
        c = [0.7462, -0.7149, -0.9401, -0.6154, -0.2524]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('R23 outside of Jones+ 15 calibration range...')
        c = [-54.1003, 13.9083, -0.8782, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 8.4):
            if prt: print ('R23 outside of Curti+ 17 calibration range...')
        c = [0.527, -1.569, -1.652, -0.421, 0.0]
    if use == 'B18': 
        x = OH
        if (OH > 8.4) | (OH < 7.8):
            if prt: print ('R23 outside of Bian+ 18 calibration range...')
        c = [138.0430, -54.8284, 7.2954, -0.32293, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_R2(OH, use = 'C17', prt = False):
    #taken from table in Patricio+
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.3) | (OH < 7.6):
            if prt: print ('R2 outside of Curti+ 17 calibration range...')
        c = [0.418, -0.961, -3.505, -1.949, 0.0]
    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_R3(OH, use = 'C17', prt = False):
    #taken from table in Patricio+
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 8.3):
            if prt: print ('R2 outside of Curti+ 17 calibration range...')
        cf = [-0.277, -3.549, -3.593, -0.981, 0.0]
    result = cf[0] * x**0. + cf[1] * x**1. + cf[2] * x**2. + cf[3] * x**3. + cf[4] * x**4.
    return result

def OH_O3(OH, use = 'M08', prt = False):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('O3 outside of Maiolino+ 08 calibration range...')
        c = [0.1549, -1.5031, -0.9790, -0.0297, 0.0]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('O3 outside of Jones+ 15 calibration range...')
        c = [-88.4378, 22.7529, -1.4501, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 8.3):
            if prt: print ('O3 outside of Curti+ 17 calibration range...')
        c = [-0.277, -3.549, -3.593, -0.981, 0.0]
    if use == 'B18': 
        x = OH
        if (OH > 8.4) | (OH < 7.8):
            if prt: print ('O3 outside of Bian+ 18 calibration range...')
        c = [43.9836, -21.6211, 3.4277, -0.1747, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_O2(OH, use = 'M08', prt = False):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('O2 outside of Maiolino+ 08 calibration range...')
        c = [0.5603, 0.0450, -1.8017, -1.8434, -0.6549]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('O2 outside of Jones+ 15 calibration range...')
        c = [-154.9571, 36.9128, -2.1921, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.3) | (OH < 7.6):
            if prt: print ('O2 outside of Curti+ 17 calibration range...')
        c = [0.418, -0.961, -3.505, -1.949, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_O32(OH, use = 'M08', prt = False):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('O32 outside of Maiolino+ 08 calibration range...')
        c = [-0.2839, -1.3881, -0.3172, 0., 0.]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('O32 outside of Jones+ 15 calibration range...')
        c = [17.9828, -2.1552, 0.0, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 7.6):
            if prt: print ('O32 outside of Curti+ 17 calibration range...')
        c = [-0.691, -2.944, -1.308, 0.0, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_S2(OH, use = 'C19', prt = False):
    #Curti+ 19 in prep
    return nan

def OH_O3S2(OH, use = 'C19', prt = False):
    #Curti+ 19 in prep
    return nan
