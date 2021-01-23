
def main():
    
    state = True
    while state:
        
        imagestack, imagebundle = bundle()
    
        # prompt user for smoothening operation and execute   
        x = input('Enter "B" for bm3d or "I" for itkdiffuse:\n\t')
        
        if x == 'B':
            refined = bm3d(imagestack, [0.1,0.1,0.1,0.1,0.1])
            state = False
            
            
            
            
            
            
            
            
        elif x == 'I':
            
            numberOfIterations = 1
            conductance = 0.5
            timeStep = 0.1
            
            inputImage = image
            
            InputPixelType = itk.D
            OutputPixelType = itk.UC
            Dimension = 3
            
            InputImageType = itk.Image[InputPixelType, Dimension]
            OutputImageType = itk.Image[OutputPixelType, Dimension]
            
            FilterType = itk.GradientAnisotropicDiffusionImageFilter[InputImageType, InputImageType]
           # gradientAnisotropicDiffusionFilter = gADF
            gADF = FilterType.New()
            gADF.SetInput(image)
            gADF.SetNumberOfIterations(numberOfIterations)
            gADF.SetTimeStep(timeStep)
            gADF.SetConductanceParameter(conductance)
    
            outputPixelTypeMinimum = itk.NumericTraits[OutputPixelType].min()
            outputPixelTypeMaximum = itk.NumericTraits[OutputPixelType].max()

            rescaler.SetOutputMinimum(outputPixelTypeMinimum)
            rescaler.SetOutputMaximum(outputPixelTypeMaximum)
            
            WriterType = itk.ImageFileWriter[OutputImageType]
            writer = WriterType.New()
            writer.SetFileName('Smoothed.TIF')
            writer.SetInput(rescaler.GetOutput())
            writer.Update()
            
            refined = rescaler.GetOutput()
    
            state = False
            
            
            
            
            
        else: print('input not recognized, please run again')
        
        return refined
    
 #  plt.figure(1)
 #  plt.imshow(itk.GetArrayViewFromImage(image), cmap='gray', vmin=20000)

 #  plt.imshow(disp_mat)
 #  plt.show()

if __name__ == '__main__':
    main()
