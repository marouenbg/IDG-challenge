import pybel

def SMILE2fp(smile):
    return pybel.readstring("smi", smile).calcfp().bits



