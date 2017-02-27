from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
import copy

def make_wcs(wcs_fn):
    #Set wcs transformation
    w = wcs.WCS(wcs_fn, naxis=2)
    
    return w

def xy_to_radec(x, y, w):
    
    #Transformation
    xy = [[x, y]]
    radec = w.wcs_pix2world(xy, 1)
    
    ra = radec[0,0]
    dec = radec[0,1]
    
    return ra, dec

def radec_to_xy(ra, dec, w):
    
    #Transformation
    radec = [[ra, dec]]
    xy = w.wcs_world2pix(radec, 1)
    
    x = xy[0,0]
    y = xy[0,1]
    
    return x, y

def xcutout_data(big_galfa_data, big_galfa_hdr, xstart = 0, xstop = None):
    
    cutout_galfa_hdr = copy.copy(big_galfa_hdr)
    
    if xstop is None:
        cutout_galfa_data = big_galfa_data[:, xstart:]
    else:
        cutout_galfa_data = big_galfa_data[:, xstart:xstop]
    
    ny, nx = cutout_galfa_data.shape
    
    # Need to change the center pixels
    old_crpix1 = big_galfa_hdr["CRPIX1"]
    old_crpix2 = big_galfa_hdr["CRPIX2"]

    cutout_crpix1 = old_crpix1 - xstart

    # Define cutout header
    cutout_galfa_hdr["NAXIS"] = 2
    cutout_galfa_hdr["NAXIS1"] = nx
    cutout_galfa_hdr["NAXIS2"] = ny
    cutout_galfa_hdr["CRPIX1"] = cutout_crpix1
    
    return cutout_galfa_hdr, cutout_galfa_data
    

def ycutout_data(big_galfa_data, big_galfa_hdr, ystart = 0, ystop = None):
    
    cutout_galfa_hdr = copy.copy(big_galfa_hdr)
    
    if ystop is None:
        cutout_galfa_data = big_galfa_data[ystart:, :]
    else:
        cutout_galfa_data = big_galfa_data[ystart:ystop, :]
    
    ny, nx = cutout_galfa_data.shape
    
    # Need to change the center pixels
    old_crpix2 = big_galfa_hdr["CRPIX2"]

    cutout_crpix2 = old_crpix2 - ystart

    # Define cutout header
    cutout_galfa_hdr["NAXIS"] = 2
    cutout_galfa_hdr["NAXIS1"] = nx
    cutout_galfa_hdr["NAXIS2"] = ny
    cutout_galfa_hdr["CRPIX2"] = cutout_crpix2
    
    return cutout_galfa_hdr, cutout_galfa_data
    
def xycutout_data(big_data, big_hdr, xstart = 0, xstop = None, ystart = 0, ystop = None):
    
    xcut_hdr, xcut_data = xcutout_data(big_data, big_hdr, xstart=x_start, xstop=x_stop)
    xycut_hdr, xycut_data = ycutout_data(xcut_data, xcut_hdr, ystart=y_start, ystop=y_stop)
    
    return xycut_hdr, xycut_data


