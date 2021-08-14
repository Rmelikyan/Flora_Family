def read_jpl_ephemeris(filePath):
    import csv
    import julian
    ''' Reads standardized JPL Ephemeris in output style 2 with km/s units
    
        Returns [Time, x, y, z, vx, vy, vz] in SI units
    '''
    Time = []
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    with open(filePath, 'r') as f:
        reader = csv.reader(f, delimiter='=')
        SOE = False
        isoe =  None
        for i, row in enumerate(reader):
            if (row!=[] and row[0] == '$$SOE'):
                SOE = True
                isoe = i+1
                continue
            if (row!=[] and row[0] == '$$EOE'):
                break
            if(SOE):
                if((i - isoe) % 3 == 0):
                    JD = float(row[0])
                    Time.append(julian.from_jd(JD, fmt = 'jd'))
                if((i - isoe) % 3  == 1):
                    x.append(float(row[1][:-3])*1000)  # Convert km to m
                    y.append(float(row[2][:-3])*1000)  # Convert km to m
                    z.append(float(row[3])*1000)       # Convert km to m
                if((i - isoe) % 3  == 2):
                    vx.append(float(row[1][:-3])*1000)  # Convert km/s to m/s
                    vy.append(float(row[2][:-3])*1000)  # Convert km/s to m/s
                    vz.append(float(row[3])*1000)       # Convert km/s to m/s
    return [Time, x, y, z, vx, vy, vz]