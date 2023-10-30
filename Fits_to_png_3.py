"""
Takes optical fits cutouts of sources and XXL X-ray mosaic to produce pdf image sets for visual inspection of sources.
False colour images are generated following method in Lupton et al 2004.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os.path
import scipy.ndimage as ndimage
import copy as cp
import matplotlib.gridspec as gridspec

# Astropy imports
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.visualization import make_lupton_rgb

# paramiters
object_id_col = 'Index'

Ra_col   = 'EXT_RA'
Dec_col  = 'EXT_DEC'
ext_col  = 'EXT'
ext_like_col = 'EXT_LIKE'
conf_col = 'Conf mean'
related_object_col = None # 'Asociated source'
class_col = 'C1C2'
z_col = None #'z'
PNT_DET_col = 'PNT_DET_ML'

parent_folder = 'V4.3_Lotto_run/'
source_file = 'Lotto_results_with_PNT_DET_ML.csv'

mosaic_source_file = "../../../../XXLN_mosaics/mosaics/XXLN_0.5_2.0_clean.fits"

Optical_hists = False

minimum = [-0.05,-0.1]
stretch = [0.4,0.8]
Q=[3,2]

# Crop image to center from current size to new smaller size
def Crop_to_Center(image, new_image_size, current_image_size):
    
    size_to_pix = len(image)/current_image_size
    image_center = len(image)/2
    min_cut = int(image_center - (new_image_size/2)*size_to_pix)
    max_cut = int(image_center + (new_image_size/2)*size_to_pix)
    
    return image[min_cut:max_cut, min_cut:max_cut]

# --- Optical data functions -------------------------------------------------------------------------------------------------

#Load in optical data
def load_hsc_images(source_id, folder = ''):
    
    #source_id = int(source_id)
    
    # load G band
    file_name = '{}Fits_files/{}_G.fits'.format(folder, source_id)
    fits_data = fits.open(file_name)
    image_G = fits_data[1].data
    image_G -= np.median(image_G)

    # load R band
    file_name = '{}Fits_files/{}_R.fits'.format(folder, source_id)
    fits_data = fits.open(file_name)
    image_R = fits_data[1].data
    image_R -= np.median(image_R)

    # load I band
    file_name = '{}Fits_files/{}_I.fits'.format(folder, source_id)
    fits_data = fits.open(file_name)
    image_I = fits_data[1].data
    image_I -= np.median(image_I)

    # load Z band
    file_name = '{}Fits_files/{}_Z.fits'.format(folder, source_id)
    fits_data = fits.open(file_name)
    image_Z = fits_data[1].data
    image_Z -= np.median(image_Z)
    
    return image_G, image_R, image_I, image_Z

#Create colour images from indevidual optical bands
def create_colour_images(G, R, I, Z, minimum = [0,0], stretch = [0.8,1.5], Q=[5,3]):

    imageGRI = make_lupton_rgb(I, R, G, minimum=minimum[0], stretch=stretch[0], Q=Q[0])
    imageRIZ = make_lupton_rgb(Z, I, R, minimum=minimum[1], stretch=stretch[1], Q=Q[1])
    
    return imageGRI, imageRIZ

# --- X-ray  data functions --------------------------------------------------------------------------------------------------

# Cutout image of source from X-ray data
def XXL_Cutout(ra, dec, field, wcs, edge_size = 5.):
    
    # Create cutout properties
    position = SkyCoord(ra, dec, unit="deg")
    size = units.Quantity((edge_size, edge_size), units.arcmin)
    
    # Return cutout as array
    return Cutout2D(field[0].data, position, size, wcs=wcs).data

# Calculate contours for ploting X-ray data
def XXL_Contours(image, bins = 500, contour_fractions = [   0.692, 0.841,   0.933, 0.977, 0.99375, 0.9985, 0.999767]):
    
    # min and max values in image, note min value limited to minimum of 0.001 to avoid divide by zero
    min_val = max(np.amin(image),0.001)
    max_val = np.amax(image)
    
    # bin data over logarithmic bins
    hist, bin_edges = np.histogram(image, bins = 10**np.linspace(np.log10(min_val),np.log10(max_val),bins))
    
    # make hist cumulative
    for i in range(1,len(hist)):
        hist[i] += hist[i-1]
    
    # normalise data
    hist = hist/hist[-1]
    
    # calculate threshold for each contour
    contours = []
    i = 0
    for frac in contour_fractions:
        # iterate till end of current hist value is larger than contour fraction
        while( i<len(hist) and hist[i] < frac ):
            i += 1
        
        # append bin center to contour
        contours.append((bin_edges[i]+bin_edges[i+1])/2)
    
    
    i = 0
    while i+1 < len(contours) and contours[i] < contours[i+1]:
        i += 1
    
    return contours[:i]
        
    
    
# --- Ploting functions ------------------------------------------------------------------------------------------------------

def Plot_Cutout(image, axs, title, size = 5., contour_image = None, contours = None, colours = None, cmap = None, y=1.05):
    
    # Plot image
    axs.imshow(image, origin = 'lower', extent = [-size/2, size/2, -size/2, size/2], cmap = cmap)
    
    # Only plot contours if given
    if(contour_image is not None):
        # get positions of each pixle relative to image
        positions = np.linspace(-size/2, size/2, len(contour_image))
        
        # if no contour colours are given use a grey scale with higher values being whiter
        if(colours is None):
            colours = np.linspace(0.6,1,len(contours)).astype(str)
        
        # plot contours
        axs.contour(positions, positions, contour_image, contours, colors = colours, origin = 'lower', alpha = 1, linewidths = 0.3)
    
    # Set plot properties
    axs.set_aspect('equal')
    axs.set_title(title, y=y)
    axs.tick_params(bottom=True, top=True, left=True, right=True)

# --- Main calculation -------------------------------------------------------------------------------------------------------





# Load sources list
sources_data = pd.read_csv('{}{}'.format(parent_folder,source_file))

print(sources_data)

# Remove any duplicate sources from the list
sources_data = sources_data.drop_duplicates()

# Reset indecies to avoid missing ones produced due to removing duplicates
sources_data = sources_data.reset_index()

# Loop through sources and remove any missing optical data
#iterate backwards to avoid problems due to removing rows
for i in range(len(sources_data[object_id_col])-1,-1,-1):

    # if file doesn't exist drop
    if( not os.path.isfile('{}Fits_files/{}_G.fits'.format(parent_folder, sources_data[object_id_col][i])) ):
        print("Droping source: {}".format(sources_data[object_id_col][i]))
        sources_data.drop(i, inplace=True)

print("Remaining source count: {}".format(len(sources_data)))




# Save an orderd list of source ids, along with Ra and Dec.
# This is done so I can create a table to add results to with the same order as the lotto objects.
sources_data[[object_id_col, Ra_col, Dec_col]].to_csv('{}Lotto_objects_orderd.csv'.format(parent_folder))

# -- Get and setup XXL field image ----------------

# Load the XXL mosaics
print("Loading XXL field")
XXLN_field = fits.open(mosaic_source_file)

# Get the world cordinate system for the XXL field
wcs = WCS(XXLN_field[0].header)

# Loop over sources creating lotto data for each
print("Creating images")
j = 0
for index, current in sources_data.iterrows():
    
    # iterate loop counter
    j += 1
    
    # get source info
    source_id = current[object_id_col]
    ra  = current[Ra_col]
    dec = current[Dec_col]
    
    if(related_object_col is None):
        Related = None
    else:
        Related = current[related_object_col]
        
    if(class_col is None):
        Class = None
    else:
        Class = current[class_col]
    
    if(ext_col is None):
        EXT = None
    else:
        EXT = current[ext_col]
    
    if(ext_like_col is None):
        EXT_LIKE = None
    else:
        EXT_LIKE = current[ext_like_col]
    
    if(conf_col is None):
        Conf = None
    else:
        Conf = current[conf_col]
    
    if(z_col is None):
        z = None
    else:
        z = current[z_col]
        
    if(PNT_DET_col is None):
        PNT_DET = None
    else:
        PNT_DET = current[PNT_DET_col]
    
    # print iteration and current source id
    print("{:d} of {:d}: {}".format(j,len(sources_data[object_id_col]),source_id))
    
    # --- X-ray data analysis ---
    
    # Get 5 by 5 cutout
    imageXXL = XXL_Cutout(ra, dec, XXLN_field, wcs)
    
    # Apply power law to raw image to make low counts easier to see
    imageXXL_pow = 1-np.power(imageXXL+1,-0.7)
    
    # Check if the X-ray image contais data
    is_Xray_data = np.amax(imageXXL)>0
    
    # If there is X-ray data process it
    if(is_Xray_data):
    
        # Smooth X-ray data
        imageXXL_smoothed = ndimage.gaussian_filter(imageXXL, sigma = 2, mode = 'constant')
        
        # Calculate contour values for ploting smoothed X-ray image
        #contours = XXL_Contours(imageXXL_smoothed)
        contours = np.median(imageXXL_smoothed)+np.power(2,np.linspace(-5,5,11))
        
        # Create 2 by 2 and 1 by 1 arcmin cutout from smoothed X-ray data
        imageXXL_smoothed_2b2 = Crop_to_Center(imageXXL_smoothed, 2, 5)
        imageXXL_smoothed_1b1 = Crop_to_Center(imageXXL_smoothed, 1, 5)
        
        # Create histograms for X-ray pixle distributions
        min_val = max(np.amin(imageXXL),0.001)
        max_val = np.amax(imageXXL)
        X_ray_bins = 10**np.linspace(np.log10(min_val),np.log10(max_val),100)
        
        # X_ray_bin_edges is the same for all
        X_ray_hist_full, X_ray_bin_edges = np.histogram(imageXXL_smoothed    , bins = X_ray_bins)
        X_ray_hist_2b2,  X_ray_bin_edges = np.histogram(imageXXL_smoothed_2b2, bins = X_ray_bins)
        X_ray_hist_1b1,  X_ray_bin_edges = np.histogram(imageXXL_smoothed_1b1, bins = X_ray_bins)
        
        # Cacl bin sizes
        X_ray_bin_sizes = X_ray_bin_edges[1:]-X_ray_bin_edges[:-1]
        
        
        
    
    # Set values used later without checking if X-ray exists to None
    else:
        imageXXL_smoothed = None
        imageXXL_smoothed_2b2 = None
        contours = None
         
    
    # --- Optical data analysis ---
    
    # Load optical bands into seperate arrays
    G, R, I, Z = load_hsc_images(source_id, folder = parent_folder)
    
    #Combine optical bands into two images
    GRI, RIZ = create_colour_images(G, R, I, Z, minimum = minimum, stretch = stretch, Q=Q)
    
    # Create 2 by 2 arcmin cutouts of optical images
    GRI_2b2 = Crop_to_Center(GRI, 2, 5)
    RIZ_2b2 = Crop_to_Center(RIZ, 2, 5)
    
    # Only needed if plotting the optical histograms
    if(Optical_hists):
        # Create 2 by 2 and 1 by 1 arcmin cutouts of optical chanels
        G_2b2 = Crop_to_Center(G, 2, 5)
        R_2b2 = Crop_to_Center(R, 2, 5)
        I_2b2 = Crop_to_Center(I, 2, 5)
        Z_2b2 = Crop_to_Center(Z, 2, 5)
        
        
        G_1b1 = Crop_to_Center(G, 1, 5)
        R_1b1 = Crop_to_Center(R, 1, 5)
        I_1b1 = Crop_to_Center(I, 1, 5)
        Z_1b1 = Crop_to_Center(Z, 1, 5)
        
        # Get max value in all optical bans so plots can be compared
        optical_max = max((np.amax(G), np.amax(R), np.amax(I), np.amax(Z)))
        
        # Create bins for optical data
        optical_bins = 10**np.linspace(0.1,np.log10(optical_max),50)
        
        # Bin optical channels for all 3 cutout sizes
        # optical_bin_edges is the same for all
        G_hist_full, optical_bin_edges = np.histogram(G,     bins = optical_bins)
        G_hist_2b2,  optical_bin_edges = np.histogram(G_2b2, bins = optical_bins)
        G_hist_1b1,  optical_bin_edges = np.histogram(G_1b1, bins = optical_bins)
        
        R_hist_full, optical_bin_edges = np.histogram(R,     bins = optical_bins)
        R_hist_2b2,  optical_bin_edges = np.histogram(R_2b2, bins = optical_bins)
        R_hist_1b1,  optical_bin_edges = np.histogram(R_1b1, bins = optical_bins)
        
        I_hist_full, optical_bin_edges = np.histogram(I,     bins = optical_bins)
        I_hist_2b2,  optical_bin_edges = np.histogram(I_2b2, bins = optical_bins)
        I_hist_1b1,  optical_bin_edges = np.histogram(I_1b1, bins = optical_bins)
        
        Z_hist_full, optical_bin_edges = np.histogram(Z,     bins = optical_bins)
        Z_hist_2b2,  optical_bin_edges = np.histogram(Z_2b2, bins = optical_bins)
        Z_hist_1b1,  optical_bin_edges = np.histogram(Z_1b1, bins = optical_bins)
        
        #Calc bin sizes
        optical_bin_sizes = optical_bins[1:] - optical_bins[:-1]
        
        #Calc max optical bin count
        optical_bin_max = max((np.amax(G_hist_full), np.amax(R_hist_full), np.amax(I_hist_full), np.amax(Z_hist_full)))
    
    # --- Ploting ---
    fig = plt.figure()
    
    # Create grid for ploting on
    # Number of boxes depends on if we plot the optical histograms
    if(Optical_hists):
        gs = gridspec.GridSpec(4,10, left=0.03, bottom=0.1, right=0.97, top=0.9, wspace=0.5, hspace=0.2, width_ratios=np.ones(10), height_ratios=np.ones(4))
    else:
        gs = gridspec.GridSpec(4,8, left=0.03, bottom=0.1, right=0.97, top=0.9, wspace=0.5, hspace=0.2, width_ratios=np.ones(8), height_ratios=np.ones(4))
    
    #GRI
    axsGRI = plt.subplot(gs[0:2,0:2])
    Plot_Cutout(GRI, axsGRI, "GRI")
    
    #GRI + X-ray contours
    axsGRI_X = plt.subplot(gs[0:2,2:4])
    Plot_Cutout(GRI, axsGRI_X, "GRI + X-ray", contour_image = imageXXL_smoothed, contours = contours)
    
    #2 by 2 arcmin GRI + X-ray contours
    axsGRI_X_2b2 = plt.subplot(gs[0:2,4:6])
    Plot_Cutout(GRI_2b2, axsGRI_X_2b2, "GRI + X-ray", size = 2., contour_image = imageXXL_smoothed_2b2, contours = contours)
    
    #RIZ
    axsRIZ = plt.subplot(gs[2:4,0:2])
    Plot_Cutout(RIZ, axsRIZ, "RIZ", y = -0.2)
    
    #RIZ + X-ray contours
    axsRIZ_X = plt.subplot(gs[2:4,2:4])
    Plot_Cutout(RIZ, axsRIZ_X, "RIZ + X-ray", contour_image = imageXXL_smoothed, contours = contours, y = -0.2)
    
    #2 by 2 arcmin RIZ + X-ray contours
    axsRIZ_X_2b2 = plt.subplot(gs[2:4,4:6])
    Plot_Cutout(RIZ_2b2, axsRIZ_X_2b2, "RIZ + X-ray", size = 2., contour_image = imageXXL_smoothed_2b2, contours = contours, y = -0.2)
    
    #Raw X-ray image
    axsX = plt.subplot(gs[0:2,6:8])
    #Plot_Cutout(imageXXL, axsX, "Raw X-ray", contour_image = imageXXL_smoothed, contours = contours, colours = 'red', cmap = 'Greys')
    Plot_Cutout(imageXXL_pow, axsX, "Raw X-ray", contour_image = imageXXL_smoothed, contours = contours, colours = 'red', cmap = 'Greys')
    
    #X-ray distribution and contour values
    axsX_dist = plt.subplot(gs[2:4,6:8])
    if(is_Xray_data):
        # Plot histograms for different regions
        axsX_dist.bar(X_ray_bin_edges[:-1], X_ray_hist_full, width=X_ray_bin_sizes, align = 'edge', zorder = 1, label = '5 by 5 arcmin')
        axsX_dist.bar(X_ray_bin_edges[:-1], X_ray_hist_2b2,  width=X_ray_bin_sizes, align = 'edge', zorder = 2, label = '2 by 2 arcmin')
        axsX_dist.bar(X_ray_bin_edges[:-1], X_ray_hist_1b1,  width=X_ray_bin_sizes, align = 'edge', zorder = 3, label = '1 by 1 arcmin')
        
        # add contour markers
        for current in contours:
            axsX_dist.plot([current,current],[0, 1.01*np.amax(X_ray_hist_full)], linestyle = ':', c = 'red', zorder = 4)
    
    # Set plot properties
    axsX_dist.set_xlim(X_ray_bin_edges[0],X_ray_bin_edges[-1])
    axsX_dist.set_xscale('log')
    axsX_dist.set_ylim(0, 1.01*np.amax(X_ray_hist_full))
    axsX_dist.set_title('Smoothed X-ray counts\nplus contour values', y = -0.25)
    
    # Only plot optical hist if requested
    if(Optical_hists):
    
        #G distribution
        axsG_dist = plt.subplot(gs[0:1,8:10])
        axsG_dist.bar(optical_bin_edges[:-1], G_hist_full, width=optical_bin_sizes, align = 'edge', zorder = 1, label = '5 by 5 arcmin')
        axsG_dist.bar(optical_bin_edges[:-1], G_hist_2b2,  width=optical_bin_sizes, align = 'edge', zorder = 2, label = '2 by 2 arcmin')
        axsG_dist.bar(optical_bin_edges[:-1], G_hist_1b1,  width=optical_bin_sizes, align = 'edge', zorder = 3, label = '1 by 1 arcmin')

        # Set plot properties
        axsG_dist.set_xlim(optical_bin_edges[0],optical_bin_edges[-1])
        axsG_dist.set_xscale('log')
        axsG_dist.set_ylim(1, 1.01*optical_bin_max)
        axsG_dist.set_yscale('log')
        axsG_dist.set_title('G', y = 0.45, x = 1.05)

        #R distribution
        axsR_dist = plt.subplot(gs[1:2,8:10])
        axsR_dist.bar(optical_bin_edges[:-1], R_hist_full, width=optical_bin_sizes, align = 'edge', zorder = 1, label = '5 by 5 arcmin')
        axsR_dist.bar(optical_bin_edges[:-1], R_hist_2b2,  width=optical_bin_sizes, align = 'edge', zorder = 2, label = '2 by 2 arcmin')
        axsR_dist.bar(optical_bin_edges[:-1], R_hist_1b1,  width=optical_bin_sizes, align = 'edge', zorder = 3, label = '1 by 1 arcmin')

        # Set plot properties
        axsR_dist.set_xlim(optical_bin_edges[0],optical_bin_edges[-1])
        axsR_dist.set_xscale('log')
        axsR_dist.set_ylim(1, 1.01*optical_bin_max)
        axsR_dist.set_yscale('log')
        axsR_dist.set_title('R', y = 0.45, x = 1.05)

        #I distribution
        axsI_dist = plt.subplot(gs[2:3,8:10])
        axsI_dist.bar(optical_bin_edges[:-1], I_hist_full, width=optical_bin_sizes, align = 'edge', zorder = 1, label = '5 by 5 arcmin')
        axsI_dist.bar(optical_bin_edges[:-1], I_hist_2b2,  width=optical_bin_sizes, align = 'edge', zorder = 2, label = '2 by 2 arcmin')
        axsI_dist.bar(optical_bin_edges[:-1], I_hist_1b1,  width=optical_bin_sizes, align = 'edge', zorder = 3, label = '1 by 1 arcmin')

        # Set plot properties
        axsI_dist.set_xlim(optical_bin_edges[0],optical_bin_edges[-1])
        axsI_dist.set_xscale('log')
        axsI_dist.set_ylim(1, 1.01*optical_bin_max)
        axsI_dist.set_yscale('log')
        axsI_dist.set_title('I', y = 0.45, x = 1.05)

        #Z distribution
        axsZ_dist = plt.subplot(gs[3:4,8:10])
        axsZ_dist.bar(optical_bin_edges[:-1], Z_hist_full, width=optical_bin_sizes, align = 'edge', zorder = 1, label = '5 by 5 arcmin')
        axsZ_dist.bar(optical_bin_edges[:-1], Z_hist_2b2,  width=optical_bin_sizes, align = 'edge', zorder = 2, label = '2 by 2 arcmin')
        axsZ_dist.bar(optical_bin_edges[:-1], Z_hist_1b1,  width=optical_bin_sizes, align = 'edge', zorder = 3, label = '1 by 1 arcmin')

        # Set plot properties
        axsZ_dist.set_xlim(optical_bin_edges[0],optical_bin_edges[-1])
        axsZ_dist.set_xscale('log')
        axsZ_dist.set_ylim(1, 1.01*optical_bin_max)
        axsZ_dist.set_yscale('log')
        axsZ_dist.set_title('Z', y = 0.45, x = 1.05)

    
    
    # figure details
    
    title = "Source {}; Ra {:.5f}, Dec {:.5f}".format(source_id, ra, dec)
    
    if(Related is not None):
        title += ", Related source {}".format(Related)
        
    if(Class is not None):
        title += ", Class {}".format(Class)
    
    if(EXT is not None):
        title += ", EXT {:.5f}".format(EXT)
    
    if(EXT_LIKE is not None):
        title += ", EXT_LIKE {:.5f}".format(EXT_LIKE)
        
    if(PNT_DET is not None):
        title += ", PNT_DET_ML {:.5f}".format(PNT_DET)
    
    if(Conf is not None):
        title += ", Confidence {:.5f}".format(Conf)
    
    if(z is not None):
        title += ", z {:.5f}".format(z)
    
    
    
    fig.suptitle(title)
    
    # Set plot size depending on if plotting optical histograms
    if(Optical_hists):
        fig.set_size_inches(18,7.8)
    else:
        fig.set_size_inches(14.8,7.8)
    
    # save fig and clear buffer
    fig.savefig("{}Spureous_images/{}.png".format(parent_folder, source_id), format='png', dpi = 300)
    plt.close(fig)
    
    #plt.show()
    
    
    
    
    
    
    
    
