# -- ©Rachel Alcraft 2021, PsuGeometry --



def tauCategory(psi, phi, NN1,NO2,NCACO2,O2C):
    '''
    These are what we are checking
    1. Is it hydrogen bonded to an oxygen 2 away?
        If so, is that oxygen planar
    2. Is there a heavy atom nearby?
        If so is that heavy atom planar
    3. Is psi 0 despite not getting through the others?
    4. Is phi 0?
    '''

    if NO2 < 3.6 and abs(NCACO2) < 20:
        if O2C <  4.5: #Then it is right hand image
            if abs(psi)<50:
                return 'A1'
            else:
                return 'A2'
        else:
            if abs(psi)<50:
                return 'A3'
            else:
                return 'A4'
    else:
        return 'X'

