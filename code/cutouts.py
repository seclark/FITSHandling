from __future__ import division
import numpy as np
#import cPickle
#import matplotlib
#matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
import astropy.coordinates as coord
import copy

def make_wcs(wcs_fn):
    #Set wcs transformation
    w = wcs.WCS(wcs_fn, naxis=2)
    
    return w

def xy_to_radec(x, y, w, origincoord=1):
    """
    origincoord : coordinate in ULH corner of image: 1 for IDL, 0 for python
    """
    
    #Transformation
    xy = [[x, y]]
    radec = w.wcs_pix2world(xy, origincoord)
    
    ra = radec[0,0]
    dec = radec[0,1]
    
    return ra, dec

def radec_to_xy(ra, dec, w, origincoord=1):
    print("translating radec to xy")
    
    #Transformation
    radec = [[ra, dec]]
    xy = w.wcs_world2pix(radec, origincoord)
    
    x = xy[0,0]
    y = xy[0,1]
    
    return x, y
    
def radecs_to_lb(ras, decs):
    """
    Transformation between lists of ras, decs, to ls, bs. Assumes ra, dec in degrees
    Conforms to astropy 0.4.3     
    """
    obj = coord.SkyCoord(ras, decs, unit = "deg", frame = "icrs")
    obj = obj.galactic
    
    ls = obj.l.degree
    bs = obj.b.degree
    
    return ls, bs

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
    
    xcut_hdr, xcut_data = xcutout_data(big_data, big_hdr, xstart=xstart, xstop=xstop)
    xycut_hdr, xycut_data = ycutout_data(xcut_data, xcut_hdr, ystart=ystart, ystop=ystop)
    
    return xycut_hdr, xycut_data

def get_xlabels_ra(hdr, skip = 1.0, origincoord=1):
    
    w = make_wcs(hdr)
    
    num = hdr["NAXIS1"]/skip
    
    xax = np.linspace(0, hdr["NAXIS1"], num)
    yax = np.zeros(len(xax))
    radec = w.wcs_pix2world(xax, yax, origincoord)
    ra = radec[0]
    dec = radec[1]
    
    return xax, ra
    
def get_ylabels_dec(hdr, skip = 1.0, origincoord=1):
    
    w = make_wcs(hdr)
    
    num = hdr["NAXIS2"]/skip
    
    yax = np.linspace(0, hdr["NAXIS2"], num)
    xax = np.zeros(len(yax))
    radec = w.wcs_pix2world(xax, yax, origincoord)
    ra = radec[0]
    dec = radec[1]
    
    return yax, dec
    
def set_all_labels(ax, xax, ra_label, yax, dec_label, roundnum=1):
    """
    convenience function for setting axis labels
    """
    ax.set_xticks(xax)
    ax.set_xticklabels(np.round(ra_label, roundnum))
    ax.set_yticks(yax)
    ax.set_yticklabels(np.round(dec_label, roundnum))
    ax.set_ylim(yax[0], yax[-1])
    
    return ax
    
def radecs_to_xy(ras, decs, w, origincoord=1):
    
    #Transformation
    xs, ys = w.all_world2pix(ras, decs, origincoord)
    
    return xs, ys
    
    
def xys_to_radec(xs, ys, w, origincoord=1):
    
    #Transformation
    ras, decs = w.all_pix2world(xs, ys, origincoord)

    return ras, decs
    
def lbs_to_radec(ls, bs):

    #Transformation. Conforms to astropy 0.4.3     
    obj = coord.SkyCoord(ls, bs, unit = "deg", frame = "galactic")
    obj = obj.icrs
    
    ras = obj.ra.degree
    decs = obj.dec.degree
    
    return ras, decs 
    
def add_b_lines(b, w):
    
    ls = np.arange(0, 100, 1)
    bs = np.zeros(len(ls)) + b
    
    ras, decs = lbs_to_radec(ls, bs)
    xs, ys = radecs_to_xy(ras, decs, w, origincoord=1)

    return xs, ys
    
def bin_by(bin_quantity, bins):
    """
    make an indexed map identifying bin number of bin_quantiy
    """

    hist_item, bin_edges = np.histogram(bin_quantity[~np.isnan(bin_quantity)], bins=bins)
    
    bin_indx = np.zeros(bin_quantity.shape, np.float_)

    # Step through bins
    for i in range(len(bin_edges) - 1):
        
        # Index each bin by bin quantify value
        bin_indx[np.where((bin_quantity > bin_edges[i]) & (bin_quantity < bin_edges[i + 1]))] = i + 1
        
    return bin_indx, bin_edges

