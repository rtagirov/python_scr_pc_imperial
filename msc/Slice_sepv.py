# program to extract data from FITS files into files to repressent a box of the stellar atmosphere, including turbulent velocity.
import pyfits
import numpy as np
import os.path
import os

# Currently using Slice_rot_3D, set vbswitch to 0 if v and B files are not required.

def Slice(locin, locout, star, mag, num, dz):
    # slice along one horizontal axis.
    # dz = height resolution is
    # 11.25 km for F3V
    # 10.03 km for G2V
    #  6 km for K0V
    #  4 km for M0V
    # for the width and length, if you wish dx=dy, the resolution is
    # 58.5938 km for F3V,
    # 17.5781 km for G2V,
    # 11.7188 km for K0V
    # 4.8828 km for M0V.
    folderin = locin+'/'+star+'/'+star+'_'+mag
    folderout = locout+'/'+ star+'/'+star+'_'+mag+'/1000/'+num+'/'
    pref = star+'_'+mag+'_'+num+'_'
    dataT = pyfits.getdata(folderin + '/eosT.' + num +  '.fits')
    dataP = pyfits.getdata(folderin + '/eosP.' + num +  '.fits')
    datarho = pyfits.getdata(folderin + '/result_0.' + num +  '.fits')
    datavz = pyfits.getdata(folderin + '/result_2.' + num +  '.fits')
    dataB = pyfits.getdata(folderin + '/result_6.' + num +  '.fits')
    """ for disk centre v line of sight = v2 (vz)
    therefore datav1 = pyfits.getdata(loc + '/result_1.' + num +  '.fits')
    datav3 = pyfits.getdata(loc + '/result_3.' + num +  '.fits') 
    and calculation for v line of sight not required
    """
    v = datavz/datarho
    maxind = dataT.shape[1]-1
    minind = 0 #300
    for xi in range(0,dataT.shape[0]):
        print(str(xi))
        f = open(folderout+pref + str(xi),'w')
        fv = open(folderout+'v/v_'+pref + str(xi),'w')
        fB = open(folderout+'B/B_'+pref + str(xi),'w')
        f.write(str(dataT.shape[2]) + '\n')
        f.write(str(dataT.shape[1]-minind) + '\n')
        fv.write(str(dataT.shape[2]) + '\n')
        fv.write(str(dataT.shape[1]-minind) + '\n')
        fB.write(str(dataT.shape[2]) + '\n')
        fB.write(str(dataT.shape[1]-minind) + '\n')
        for yi in range(0,dataT.shape[2]):
            for zi in range(maxind,minind-1,-1):
                #print(str(zi))
                f.write('{0:12.6}{1:12.6}{2:12.6}{3:12.6}\n'.format(dz*(maxind-zi), dataT[xi,zi,yi], dataP[xi,zi,yi], datarho[xi,zi,yi]))
                fv.write(str(v[xi,zi,yi])+"\n")
                fB.write(str(dataB[xi,zi,yi])+"\n")
        f.close()
        fv.close()
        fB.close()
        """print(str(dataT[0,zi,0]) + '   ' + str(dataP[0,zi,0]) + '   ' + str(datarho[0,zi,0]) + '   ' + str(vturb[0,zi,0]))"""
    return

def Slice_rot(locin, locout, star, mag, num, dz):
    # Slice along the other horizontal axis. Generally now used as matches with Kok Lengs slices.
    # dz = height resolution is 11.25 km for F3V
    # 10 km for G2V
    #  6 km for K0V
    #  4 km for M0V
    # for the width and length, if you wish dx=dy, the resolution is
    # 58.5938 km for F3V,
    # 17.5781 km for G2V,
    # 11.7188 km for K0V 4.8828 km for M0V.
    #star = 'G2'
    #mag = '500G'
    folderin = locin+'/'+star+'/'+star+'_'+mag
    folderout = locout+'/'+star+'/'+star+'_'+mag+'/1000/'+num+'/'
    if os.path.exists(folderout): 
        print(folderout)
    else:
        os.makedirs(folderout)
        os.makedirs(folderout+'/v')
        if mag!='hydro':
            os.makedirs(folderout+'/B')
        print('created:')
        print(folderout)
    pref = star+'_'+mag+'_'+num+'_'
    dataT = pyfits.getdata(folderin + '/eosT.' + num +  '.fits')
    dataP = pyfits.getdata(folderin + '/eosP.' + num +  '.fits')
    datarho = pyfits.getdata(folderin + '/result_0.' + num +  '.fits')
    datavz = pyfits.getdata(folderin + '/result_2.' + num +  '.fits')
    print(dataT.shape)
    if mag!='hydro':
        dataB = pyfits.getdata(folderin + '/result_6.' + num +  '.fits')
    """ for disk centre v line of sight = v2 (vz)
    therefore datav1 = pyfits.getdata(loc + '/result_1.' + num +  '.fits')
    datav3 = pyfits.getdata(loc + '/result_3.' + num +  '.fits') 
    and calculation for v line of sight not required
    """
    v = datavz/datarho
    maxind = dataT.shape[1]-1
    minind = 0
    for xi in range(0,dataT.shape[2]):#0,dataT.shape[2]):
        print(str(xi))
        f = open(folderout+pref + 'rot_'+str(xi),'w')
        fv = open(folderout+'v/v_'+pref + 'rot_' + str(xi),'w')
        if mag!='hydro':
            fB = open(folderout+'B/B_'+pref + 'rot_' + str(xi),'w')
        f.write(str(dataT.shape[0]) + '\n')
        f.write(str(dataT.shape[1]-minind) + '\n')
        fv.write(str(dataT.shape[0]) + '\n')
        fv.write(str(dataT.shape[1]-minind) + '\n')
        if mag!='hydro':
            fB.write(str(dataT.shape[0]) + '\n')
            fB.write(str(dataT.shape[1]-minind) + '\n')
        for yi in range(0,dataT.shape[0]):
            for zi in range(maxind,minind-1,-1):
                #print(str(zi))
                f.write('{0:12.6}{1:12.6}{2:12.6}{3:12.6}\n'.format(dz*(maxind-zi), dataT[yi,zi,xi], dataP[yi,zi,xi], datarho[yi,zi,xi]))
                fv.write(str(v[yi,zi,xi])+"\n")
                if mag!='hydro':
                    fB.write(str(dataB[yi,zi,xi])+"\n")
        f.close()
        fv.close()
        if mag!='hydro':
            fB.close()
        """print(str(dataT[0,zi,0]) + '   ' + str(dataP[0,zi,0]) + '   ' + str(datarho[0,zi,0]) + '   ' + str(vturb[0,zi,0]))"""
    return

def Slice_rot_3D(locin, locout, star, mag, num, dz, vbswitch, minind=0, endstr=''):
    """
     Extracts relevant infro from eosT, eosP and result MURaM files and writes out slice files ready to be used by modcon.
     e.g. ss.Slice_rot_3D('/media/BENJAMINS SOLAR CUBES/New_Solar_Cubes_untar','/media/Stellar3/slices_v','M0','hydro','030000',4.0,1)
     Same as Slice_rot, but slices v and B with all 3 dimensions, if vbswitch is 1.
     dz = height resolution: 
     1.25 km for F3V
     10 km for G2V
     6 km for K0V
     4 km for M0V
     3.2 km for M2V
     vbswitch  = 1 to extract 3D v and B files
     minind = 0 (for most stars) = 300 (for F3 stars) how many points to cut off the bottom.

     for the width and length, dx=dy, the resolution is
     58.5938 km for F3V,
     17.5781 km for G2V,
     11.7188 km for K0V, 
     4.8828 km for M0V,
     3.0469 km for M2V.
    """
    folderin = locin+'/'+star+'/'+star+'_'+mag
    folderout = locout+'/'+star+'/'+star+'_'+mag+'/1000/'+num+'/'
    if os.path.exists(folderout): 
        print(folderout)
    else:
        os.makedirs(folderout)
        print('created:')
        print(folderout)
    pref = star+'_'+mag+'_'+num+'_'
    dataT = pyfits.getdata(folderin + '/eosT.' + num +  '.fits')
    dataP = pyfits.getdata(folderin + '/eosP.' + num +  '.fits')
    datarho = pyfits.getdata(folderin + '/result_0.' + num +  '.fits')
    print(dataT.shape)
    maxind = dataT.shape[1]-1
    for xi in range(0,dataT.shape[2]):#0,dataT.shape[2]):
        print(str(xi))
        f = open(folderout+pref + 'rot_'+endstr+str(xi),'w')
        f.write(str(dataT.shape[0]) + '\n')
        f.write(str(dataT.shape[1]-minind) + '\n')
        #print(str(dataT.shape[1]-minind))
        for yi in range(0,dataT.shape[0]):
            for zi in range(maxind,minind-1,-1):
                #print(str(zi))
                f.write('{0:12.6}{1:12.6}{2:12.6}{3:12.6}\n'.format(dz*(maxind-zi), dataT[yi,zi,xi], dataP[yi,zi,xi], datarho[yi,zi,xi]))

        f.close()
    if vbswitch:
        Slice_rot_v(locin, locout, star, mag, num, datarho, minind, endstr)
        if mag!='hydro':
            Slice_rot_B(locin, locout, star, mag, num, minind, endstr)
    return

def Slice_rot_B(locin, locout, star, mag, num, minind, endstr=''):
    folderin = locin+'/'+star+'/'+star+'_'+mag
    folderout = locout+'/'+star+'/'+star+'_'+mag+'/1000/'+num+'/B'
    if os.path.exists(folderout): 
        print(folderout)
    else:
        os.makedirs(folderout)
        print('created:')
        print(folderout)
    pref = star+'_'+mag+'_'+num+'_'
    dataBx = pyfits.getdata(folderin + '/result_5.' + num +  '.fits')
    dataBz = pyfits.getdata(folderin + '/result_6.' + num +  '.fits')
    dataBy = pyfits.getdata(folderin + '/result_7.' + num +  '.fits')
    maxind = dataBx.shape[1]-1
    for xi in range(0,dataBx.shape[2]):
        fB = open(folderout+'/B_'+pref + 'rot_'+endstr + str(xi),'w')
        fB.write(str(dataBx.shape[0]) + '\n')
        fB.write(str(dataBx.shape[1]-minind) + '\n')
        fB.write('x \t y \t z \n')
        for yi in range(0,dataBx.shape[0]):
            for zi in range(maxind,minind-1,-1):
                fB.write(str(dataBx[yi,zi,xi])+'\t'+str(dataBy[yi,zi,xi])+'\t'+str(dataBz[yi,zi,xi])+"\n")
        fB.close()
    return

def Slice_rot_v(locin, locout, star, mag, num, datarho, minind, endstr=''):
    folderin = locin+'/'+star+'/'+star+'_'+mag
    folderout = locout+'/'+star+'/'+star+'_'+mag+'/1000/'+num+'/v'
    if os.path.exists(folderout): 
        print(folderout)
    else:
        os.makedirs(folderout)
        print('created:')
        print(folderout)
    pref = star+'_'+mag+'_'+num+'_'
    datavx = pyfits.getdata(folderin + '/result_1.' + num +  '.fits')
    datavz = pyfits.getdata(folderin + '/result_2.' + num +  '.fits')
    datavy = pyfits.getdata(folderin + '/result_3.' + num +  '.fits')
    maxind = datavx.shape[1]-1
    vx = datavx/datarho
    vz = datavz/datarho
    vy = datavy/datarho
    for xi in range(0,datavx.shape[2]):
        print(xi)
        fv = open(folderout+'/v_'+pref + 'rot_'+endstr + str(xi),'w')
        fv.write(str(datavx.shape[0]) + '\n')
        fv.write(str(datavx.shape[1]-minind) + '\n')
        fv.write('x \t y \t z \n')
        for yi in range(0,datavx.shape[0]):
            for zi in range(maxind,minind-1,-1):
                fv.write(str(vx[yi,zi,xi])+'\t'+str(vy[yi,zi,xi])+'\t'+str(vz[yi,zi,xi])+"\n")
        fv.close()
    return

def writefile(filenam, num, dz, dataT, dataP, datarho):
    """
    dz = height resolution is 
      11.25 km for F3V
      10 km for G2V
       6 km for K0V
       4 km for M0V
      for the width and length, if you wish dx=dy, the resolution is
      58.5938 km for F3V,
      17.5781 km for G2V,
      11.7188 km for K0V 4.8828 km for M0V.
    star = 'G2'
    mag = 'hydro'
    folder = './slices_v/'+star+'/'+star+'_'+mag+'/1000/'+num+'/'
    pref = star+'_'+mag+'_'+num+'_'
    for disk centre v line of sight = v2 (vz)
    therefore datav1 = pyfits.getdata(loc + '/result_1.' + num +  '.fits')
    datav3 = pyfits.getdata(loc + '/result_3.' + num +  '.fits') 
    and calculation for v line of sight not required
    """
    f = open(folder+pref+'avg2','w')
    f.write(str(1) + '\n')
    f.write(str(dataT.shape[0]) + '\n')
    for zi in range(0,dataT.shape[0]):
        #print(str(zi))
        f.write('{0:12.6}{1:12.6}{2:12.6}{3:12.6}\n'.format(dz*(zi), dataT[zi], dataP[zi], datarho[zi]))
    f.close()
    """print(str(dataT[0,zi,0]) + '   ' + str(dataP[0,zi,0]) + '   ' + str(datarho[0,zi,0]) + '   ' + str(vturb[0,zi,0]))"""
    return

def writefile_newres(filenam, num, dz, numx, numy, numdepths):
    # dz = height resolution is 11.25 km for F3V
    # 10 km for G2V
    #  6 km for K0V
    #  4 km for M0V
    # for the width and length, if you wish dx=dy, the resolution is
    # 58.5938 km for F3V,
    # 17.5781 km for G2V,
    # 11.7188 km for K0V 4.8828 km for M0V.
    star = 'G2'
    mag = '500G'
    folder = './slices_v/'+star+'/'+star+'_'+mag+'/1000/'+num+'/'
    pref = star+'_'+mag+'_'+num+'_'
    """ for disk centre v line of sight = v2 (vz)
    therefore datav1 = pyfits.getdata(loc + '/result_1.' + num +  '.fits')
    datav3 = pyfits.getdata(loc + '/result_3.' + num +  '.fits') 
    and calculation for v line of sight not required
    """
    for ix in range(0,numx):
        data = np.genfromtxt(folder+pref+str(ix),dtype = 'float',skip_header = 2)
        f = open(folder+pref+'_newres_'+str(ix),'w')
        f.write(str(numy) + '\n')
        f.write(str(numdepths) + '\n')
        for yi in range(0,numy):
            start = yi*numdepths
            for zi in range(0,numdepths):
                #print(str(zi))
                f.write('{0:12.6}{1:12.6}{2:12.6}{3:12.6}\n'.format(dz*(zi), data[start+zi,1], data[start+zi,2], data[start+zi,3]))
        f.close()
    """print(str(dataT[0,zi,0]) + '   ' + str(dataP[0,zi,0]) + '   ' + str(datarho[0,zi,0]) + '   ' + str(vturb[0,zi,0]))"""
    return
