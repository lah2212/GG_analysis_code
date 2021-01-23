from PIL import Image
import numpy as np
from numpy import asarray
import itk
    
    
# Bundles a selection of images and saves them as an array
def bundle():
    
# initialize the imageset and select the number of images to be bundled:
    n = 5
    imageset = np.zeros((4096,4096,n))
# iterate through the images and fill the array; (n x (X x Y))
    for i in range(n):
       filename = 'Pt_94kx_Conical 5sec_1fs_70umObj_5frames_02_1.ser_'+str(96+i)+'.png'
       imageset[:,:,i] = asarray(Image.open(filename))
# save imageset as multichannel bundle
    np.save("Merged_image(temp)",imageset, allow_pickle=True)
    
    
# construct a 3D image out of the slices
    
    n = 5
    
    PixelType = itk.F
    InputImageType = itk.Image[PixelType, 2]
    OutputImageType = itk.Image[PixelType, 3]
    
    reader = itk.ImageFileReader[InputImageType].New()
    tileFilter = itk.TileImageFilter[InputImageType, OutputImageType].New()
    
    layout = [1, 1, 0]
    tileFilter.SetLayout(layout)
    
    for i in range(n):
        filename = 'Pt_94kx_Conical 5sec_1fs_70umObj_5frames_02_1.ser_'+str(96+i)+'.png'
        reader.SetFileName(filename)
        reader.Update()
        
        inputImage = reader.GetOutput()
        inputImage.DisconnectPipeline()
        
        tileFilter.SetInput(i, inputImage)
        
    defaultValue = 128
    tileFilter.SetDefaultPixelValue(defaultValue)
    tileFilter.Update()
    
    writer = itk.ImageFileWriter[OutputImageType].New()
    writer.SetFileName('Bundle(temp).jpg')
    writer.SetInput(tileFilter.GetOutput())
    writer.Update()
    image = tileFilter[1]

    return imageset, image