# __author__ = 'NataliaDelCoco'


#Projecao 

# tirado da biblioteca skymapper
# https://github.com/pmelchior/skymapper
# usa a projecao Albers

# criado: 21/03/19

from numpy import *

DEG2RAD=pi/180.0
#==  funcoes  =======================================
def init_params(ra,dec):
  ra_ = array(ra)
  ra_[ra_ > 180] -= 360
  ra_[ra_ < -180] += 360
  # weigh more towards the poles because that decreases distortions
  ra0 = (ra_ * dec).sum() / dec.sum()
  if ra0 < 0:
      ra0 += 360
  dec0 = median(dec)
  # determine standard parallels
  dec1, dec2 = dec.min(), dec.max()
  delta_dec = (dec0 - dec1, dec2 - dec0)
  dec1 += delta_dec[0]/6
  dec2 -= delta_dec[1]/6

  return ra0, dec0, dec1, dec2


def nC(dec_1,dec_2):

  n = (sin(dec_1 * DEG2RAD) + sin(dec_2 * DEG2RAD)) / 2
  C = cos(dec_1 * DEG2RAD)**2 + 2 * n * sin(dec_1 * DEG2RAD)
  return n, C


def rho_f(n,C,dec):
  aux=sqrt(C - 2 * n * sin(dec * DEG2RAD)) / n

  return aux


def _toArray(x):
    """Convert x to array if needed
    Returns:
        array(x), boolean if x was an array before
    """
    if isinstance(x, ndarray):
        return x, True
    if hasattr(x, '__iter__'):
        return array(x), True
    return array([x]), False



def wrapRA(ra_0, ra):
    """Normalize rectascensions to -180 .. 180, with reference `ra_0` at 0"""
    ra_, isArray = _toArray(ra)
    ra_ = ra_0 - ra_ # inverse for RA
    # check that ra_aux is between -180 and 180 deg
    ra_[ra_ < -180 ] += 360
    ra_[ra_ > 180 ] -= 360
    if isArray:
        return ra_
    return ra_[0]

#== corpo ===========================================

#considera que ra e dec sao dois arrays
def proj(ra,dec):

  #pra achar os parametros iniciais
  #ra_0 =RA that maps onto x = 0
  #dec_0 = Dec that maps onto y = 0 
  #dec_1 = lower standard parallel
  #dec_2 = upper standard parallel (must not be -dec_1)
  ra_0, dec_0, dec_1, dec_2 = init_params(ra,dec)



  #passo 1 => constantes
  n, C = nC(dec_1,dec_2)
  rho_0=rho_f(n,C,dec_0)

  #passo 2
  ra_ = wrapRA(ra_0,ra)
  theta = n*ra_
  rho = rho_f(n,C,dec)

  #passo 3
  Xc = rho*sin(theta * DEG2RAD) 
  Yc = rho_0 - rho*cos(theta * DEG2RAD)

  return Xc,Yc