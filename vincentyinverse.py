# Vincenty Inverse geodetic computation function (Andre Verville)
# generates bearing and distance from two given points in lat lon
# Formulas from my Bachelor Degrees thesis, Laval University, 1978
# Credits: Thaddeus Vincenty, Survey Review, Vol.XXIII, No. 176, April 1975
# https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
# ----------------------------------------------------------------------------

def vincentyinverse(lat1deg, lon1deg, lat2deg, lon2deg):
    import math
    a, b = 6378137.0, 6356752.3141 # GRS80 ellipsoid, edit as required

    lat1rad = math.radians(lat1deg)
    lon1rad = math.radians(lon1deg)
    lat2rad = math.radians(lat2deg)
    lon2rad = math.radians(lon2deg)

    e2 = ((a * a) - (b * b)) / (a * a)
    sinlat1 = math.sin(lat1rad)
    coslat1 = math.cos(lat1rad)
    sinlat2 = math.sin(lat2rad)
    coslat2 = math.cos(lat2rad)

    deltalon = lon2rad - lon1rad
    N1 = a / ((1.0 - e2 * (sinlat1 * sinlat1))) ** 0.5
    M1 = (a * (1.0 - e2)) / ((1 - (e2 * sinlat1 * sinlat1)) ** 1.5)
    N2 = a / ((1.0 - e2 * (sinlat2 * sinlat2))) ** 0.5
    M2 = (a * (1.0 - e2)) / ((1 - (e2 * sinlat2 * sinlat2)) ** 1.5)
    X1 = N1 * coslat1
    Z1 = N1 * (1.0 - e2) * sinlat1
    X2 = N2 * coslat2
    Z2 = N2 * (1.0 - e2) * sinlat2
    Y = math.sin(deltalon) * X2
    temp = (N1 * coslat1) - (math.cos(deltalon) * X2)
    X = (temp * sinlat1) - ((Z1 - Z2) * coslat1)
    brg12 = math.atan2(Y, X)
    Y = math.sin(deltalon) * X1
    temp = (N2 * coslat2) - (math.cos(deltalon) * X1)
    X = ((Z2 - Z1) * coslat2) - (temp * sinlat2)
    brg21 = math.atan2(Y, X)
    X2P = math.cos(deltalon) * X2
    Y2 = math.sin(deltalon) * X2
    d12 = (((X2P - X1) * (X2P - X1)) + (Y2 * Y2) + ((Z2 - Z1) * (Z2 - Z1)))
    d12 = d12 ** 0.5
    brgb12 = (brg12 + brg21) / 2.0
    sinsqb = (math.sin(brgb12)) ** 2.0
    cossqb = (math.cos(brgb12)) ** 2.0
    Y = (N1 + N2) * (M1 + M2)
    X = (((M1 + M2) * sinsqb) + ((N1 + N2) * cossqb)) * 2.0
    RB = Y / X
    dist = d12 + ((d12 * d12 * d12) / (24.0 * RB * RB))
    brg12deg = math.degrees(brg12)
    brg21deg = math.degrees(brg21)
    if brg12deg < 0.0: brg12deg = brg12deg + 360.0
    if brg21deg < 0.0: brg21deg = brg21deg + 360.0
    brg21deg = (brg21deg + 180.0) % 360.0
    return brg12deg, brg21deg, dist
