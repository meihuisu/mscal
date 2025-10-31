import sys
import netCDF4 as nc
import numpy as np
from scipy.interpolate import interpn


def load_ncdata(fpath, varnm):
    ds = nc.Dataset(fpath)
    loaddata = ds[varnm][:]
    var=loaddata.data
    return var

def interpolate_mod(loni, lati, depi, vari, lonq, latq, depq):
    points = (depi, lati, loni)
    pointq = np.vstack((depq, latq, lonq)).T
    varq = interpn(points, vari, pointq, method='linear', bounds_error=False, fill_value=np.nan)
    return varq


args=sys.argv
ncfpath=args[1]
lonq, latq, depq = float(args[2]), float(args[3]), float(args[4])

loni = load_ncdata(ncfpath, 'longitude')
lati = load_ncdata(ncfpath, 'latitude')
depi = load_ncdata(ncfpath, 'depth')
vpi = load_ncdata(ncfpath, 'vp')
vsi = load_ncdata(ncfpath, 'vs')
rhoi = load_ncdata(ncfpath, 'rho')
dep_scale_factor=111000.
vpq=interpolate_mod(loni, lati, depi/dep_scale_factor, vpi, lonq, latq, depq/dep_scale_factor)
vsq=interpolate_mod(loni, lati, depi/dep_scale_factor, vsi, lonq, latq, depq/dep_scale_factor)
rhoq=interpolate_mod(loni, lati, depi/dep_scale_factor, rhoi, lonq, latq, depq/dep_scale_factor)
print(f"Vp, Vs, Density = {vpq} m/s, {vsq} m/s, {rhoq} kg/m^3")
