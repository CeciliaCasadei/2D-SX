# -*- coding: utf-8 -*-
import joblib
import numpy
selectedRun = '0127'

resolutionLimit = 7.0 # A

cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I    

lattices = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'
                           %(selectedRun, selectedRun, selectedRun))
print len(lattices)

lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'%(selectedRun, selectedRun, selectedRun), 'r')
lattices_names = list(lattices_names) 
print len(lattices_names)  
 
images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(selectedRun, selectedRun))
images_names = list(images_names)

fOpen = open('./Output_r%s/r%s_dataToEMC.txt'%(selectedRun, selectedRun), 'w')

n = 0
for index in range(0, len(lattices_names)):
    n = n+1
    lattice_name = lattices_names[index]
    runNumber = lattice_name[10:14]
    imageNumber = int(lattice_name[80:84])
    latticeNumber = lattice_name[92]
    latticeMatrix = lattices[index]                    # h k q_rod I_corrected flag=1 i_unassembled j_unassembled
    image_name = images_names[imageNumber-1].split()[1]
    print '*********************'
    print lattice_name
    print runNumber, imageNumber, latticeNumber
    print image_name
    fOpen.write('----- Begin chunk -----\n')
    fOpen.write('Image filename: \n')
    fOpen.write('%s\n'%image_name)
    fOpen.write('--- Begin crystal\n')
    fOpen.write('Crystal label: %s'%lattice_name)
    fOpen.write('Reflections measured after indexing\n')
    fOpen.write('  h   k   l          I    phase   sigma(I)  counts  fs/px  ss/px\n')
    for spot in latticeMatrix:
        h = int(spot[0])
        k = int(spot[1])
        I = float(spot[[3]])
        if numpy.isnan(I):
            continue
        
        reciprocalVector = [h, k]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
        resolution = 2* numpy.pi / q_2D     # A
        if resolution < resolutionLimit:
            continue
        
        fOpen.write('%3d%4d   0%11.2f        -          0       1      0      0\n'%(h, k, I))   
    fOpen.write('End of reflections\n')
    fOpen.write('--- End crystal\n')
    fOpen.write('----- End chunk -----\n')
    
fOpen.close()
print n