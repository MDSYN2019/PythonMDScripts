#! /usr/bin/env python
import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy
from scipy import ndimage
from skimage import filters
import skimage
from skimage import feature

import matplotlib.pyplot as plt

def render_image(image,colorMap,colorRange,ticks,filename):

    plt.figure()
    a=plt.imshow(image.transpose(), cmap=colorMap,origin='lower',extent=[0,numberPixels[0],0,numberPixels[1]])
    plt.axis('off')    
    a.set_clim(colorRange)        
    plt.savefig(filename,dpi=dpi,transparent=True,bbox_inches='tight')
    plt.close()

# protect the runnable part script to allow e.g. testing etc
if __name__ == "__main__":
    
    outputStem="test"
    dpi=150
    pixelWidth=1
    canningSigma=2.0
    gaussianSigma=8.0
        
    u = MDAnalysis.Universe("analysis.pdb")
        
    # create an array with the maximum number of pixels (voxels)
    originalNumberPixels = numpy.array([int((u.dimensions[0]*1.1) * pixelWidth),int((u.dimensions[1]*1.1) * pixelWidth),int((u.dimensions[2]*1.1) * pixelWidth)])

    speciesList = ['unsaturated','saturated','cholesterol']

    # find out where the centre of the bilayer lies
    bilayer = u.select_atoms("name PO4 or name ROH")
    bilayerCentre = bilayer.center_of_geometry()[2]
    
    # use the MDAnalysis analysis function LeafletFinder to identify the two bilayers
    bilayerLeaflets = LeafletFinder(u,'name PO4 or name ROH')
    bilayerLeaflet = {}

    # check to see where the first leaflet lies
    if bilayerLeaflets.group(0).centroid()[2] > bilayerCentre:
        leafletList = ["upper","lower"]
    else:
        leafletList = ["lower","upper"]            
    
    imageLookup = {}
    imageLookup['unsaturated'] = {'colour': 'Blues', 'range': [0,0.0240], 'sigma': 4.1, 'ticks': [0,0.012,0.024]}
    imageLookup['saturated'] = {'colour': 'Reds', 'range': [0,0.0240], 'sigma': 4.1, 'ticks': [0,0.012,0.024]}
    imageLookup['cholesterol'] = {'colour': 'Purples', 'range': [0,0.0240], 'sigma': 4.1, 'ticks': [0,0.012,0.024]}
    imageLookup['mask'] = {'colour': 'RdBu', 'range': [-0.2,1.2], 'ticks': [0,1]}
    imageLookup['edges'] = {'colour': 'Greys', 'range': [0,1], 'ticks': [0,1]}
    
            
    # calculate how wide the bilayer is in pixels
    numberPixels = numpy.array([int((u.dimensions[0]+2) * pixelWidth),int((u.dimensions[1]+2) * pixelWidth),int((u.dimensions[2]+2) * pixelWidth)])

    # hence calculate the total number of pixels in each image
    totalPixels = numberPixels[0] * numberPixels[1]

    # create a pixel coordinate mesh grid
    P,Q = numpy.mgrid[0:(numberPixels[0]), 0:(numberPixels[1])]
                        
    # define the dictionaries we will use to store the arrays
    Coordinates = {}
    Pixels = {}
    PixelMesh = {}        
    DensityMesh = {}
    ThicknessMesh = {}

    # initialise the XYZCoordinates
    for species in ['unsaturated','saturated','cholesterol']:
        for leaflet in ["upper","lower"]:
            Coordinates[(species,leaflet)] = []
    
    # extract the coordinates of the lipid species present
    for i in [0,1]:

        # define saturated and unsaturated phases 
        bilayerLeaflet[('saturated',i)] = bilayerLeaflets.groups(i).select_atoms("resname DPPC and name PO4")
        bilayerLeaflet[('cholesterol',i)] = bilayerLeaflets.groups(i).select_atoms("resname CHOL and name ROH")
        bilayerLeaflet[('unsaturated',i)] = bilayerLeaflets.groups(i).select_atoms("resname DFPC and name PO4")
                    
        for species in ['unsaturated','saturated','cholesterol']:
            Coordinates[(species,leafletList[i])] = pixelWidth*bilayerLeaflet[(species,i)].positions

    # now process the coordinates
    for leaflet in leafletList:

        for species in speciesList:
            
            frontstem = "png/image-"  
            endstem = "-" + species + "-" + leaflet + "-" + outputStem + ".png"

            # create a list of the coordinates in pixel format
            Pixels[(species,leaflet)] = Coordinates[(species,leaflet)].astype(int)
            
            # wrap the values in the list using modulus to check no pixels are set out of range
            Pixels[(species,leaflet)] = numpy.mod(Pixels[(species,leaflet)],numberPixels)
            
            # create empty mesh arrays populated with zeros
            PixelMesh[(species,leaflet)] = numpy.zeros((numberPixels[0],numberPixels[1]), dtype=int)
            DensityMesh[(species,leaflet)] = numpy.zeros((numberPixels[0],numberPixels[1]))
            
            # set values to unity where a molecular species is present
            PixelMesh[(species,leaflet)][Pixels[(species,leaflet)][:,0],Pixels[(species,leaflet)][:,1]]=1    

            # render an image of the pixels to disc
            render_image(PixelMesh[(species,leaflet)],imageLookup[species]['colour'],imageLookup[species]['range'],imageLookup[species]['ticks'],frontstem + "pixels" + endstem)
            
            # convolve the pixel mesh array with a Gaussian with an appropriate width        
            DensityMesh[(species,leaflet)] = ndimage.gaussian_filter(PixelMesh[(species,leaflet)].astype(float),imageLookup[species]['sigma'])

            # render an image of the density to disc
            render_image(DensityMesh[(species,leaflet)],imageLookup[species]['colour'],imageLookup[species]['range'],imageLookup[species]['ticks'],frontstem + "density" + endstem)
        
        frontstem = "png/image-"  
        endstem = "-" + leaflet + "-" + outputStem + ".png"
        
        for i in ['difference','mask','edges']:
            DensityMesh[(i,leaflet)] = numpy.zeros((numberPixels[0],numberPixels[1]))        

        # create a image where the pixels of the saturated phase is subtracted from the unsaturated phase
        DensityMesh[('difference',leaflet)] = DensityMesh[('unsaturated',leaflet)] - DensityMesh[('saturated',leaflet)]
        
        # now we can create a binary mask using the average value of the above as a cutoff
        DensityMesh[('mask',leaflet)] = DensityMesh[('difference',leaflet)] > DensityMesh[('difference',leaflet)].mean()                    
        render_image(DensityMesh[('mask',leaflet)],imageLookup['mask']['colour'],imageLookup['mask']['range'],imageLookup['mask']['ticks'],frontstem + "mask" + endstem)

        # the edges (as single pixel lines) can be detected using a Canny filter
        # note: the input option canningSigma controls the "smoothness" of the edge detection
        DensityMesh[('edges',leaflet)] = feature.canny(DensityMesh[('mask',leaflet)].astype(float),sigma=canningSigma)
        render_image(DensityMesh[('edges',leaflet)],imageLookup['edges']['colour'],imageLookup['edges']['range'],imageLookup['edges']['ticks'],frontstem + "edges" + endstem)                    
            


