import pandas as pd
import numpy as np
from astropy.io import fits

# Load the ids of the available galaxies
non_AGN_list = np.loadtxt('../Text_Files/nonAGN.txt')

seyfert_before = np.loadtxt('../Text_Files/seyfert_before.txt')
sf_before = np.loadtxt('../Text_Files/sf_before.txt')
liner_before = np.loadtxt('../Text_Files/liner_before.txt')

AGN_list = np.append(seyfert_before, np.append(sf_before,liner_before,axis=0),axis=0)

full_list = np.append(non_AGN_list, AGN_list,axis=0)

def retrieve_map(galaxies, filename, percentile):
    """ Retrieves distributions of values calculated from fits of spectra and
    returns particular percentile (or mean) of that distribution for each
    galaxy.
    
    galaxies: list of galaxy ids to iterate through
    
    filename: name of file that contains particular calculated values
    
    percentile: which percentile value to select from distribution
    
    @returns map_vals: disctribution of percentile values
    """
    map_vals = []
    for gal in galaxies:
        plateid = str(int(gal[0]))
        ifuid = str(int(gal[1]))
        
        path = "../Fits_Files/MPL_5/"+plateid+"/"+ifuid+"/"+plateid+"-"+ifuid+ "-" + filename
        currentHdulist = fits.open(path)
        vals = currentHdulist[0].data
        
        '''SNRPath = "../Fits_Files/MPL_5/"+plateid+"/"+ifuid+"/"+plateid+"-"+ifuid+ "-SNR_Map_Final.fits"
        currentHdulist = fits.open(SNRPath)
        SNRVals = currentHdulist[0].data'''
        
        nonNan = vals[~np.isnan(vals)]
        '''SNRVals = SNRVals[~np.isnan(vals)]
        vals = nonNan[SNRVals > 10]'''
        if np.size(nonNan) != 0:
            map_vals.append(np.percentile(nonNan,percentile))
            #map_vals.append(np.mean(nonNan))
        else:
            map_vals.append(np.nan)
    return map_vals
    
def retrieve_luminosity(plateid, ifuid, percentile):
    """ Retrieves luminosities from drpall fits file.
    
    plateid: plate id for particular galaxy
    
    ifuid: ifu id for particular galaxy
    
    percentile: percentile value of luminosity to return
    
    @returns percentile of luminosity (on NaN)
    """
    drpall = fits.open('../../drpall-v2_0_1.fits')
    data = drpall[1].data
    plates = np.where(data.plate == int(plateid))
    for plate in plates[0]:
        if data[plate]['ifudsgn'] == str(ifuid):
            luminosities = data[plate]['nsa_sersic_absmag']
            return np.percentile(luminosities, percentile)
    return np.nan

def lum_array(galaxies, percentile):
    """ Produces distribution of luminosity values for list of galaxies.
    
    galaxies: list of galaxy ids to iterate through
    
    percentile: which percentile value to select from distribution 
    
    @returns array of luminosities
    """
    lums = []
    for gal in galaxies:
        plateid = str(int(gal[0]))
        ifuid = str(int(gal[1]))
        gal_lum = retrieve_luminosity(plateid, ifuid, percentile)
        lums.append(gal_lum)
    return np.array(lums)


#full_lums_75 = lum_array(full_list, 75)
full_flux_25 = np.array(retrieve_map(full_list, 'Flux_Map_Final.fits', 25))
full_W80_25 = np.array(retrieve_map(full_list, 'W80_Map_Final.fits', 25))
full_vmed_25 = np.array(retrieve_map(full_list, 'Vmed_Map_Final.fits', 25))
#full_mean_gauss = np.array(retrieve_map(full_list, 'Gaussianity_Map_Final.fits', 75))
#print(non_AGN_flux_75.shape)
#print(non_AGN_W80_75.shape)