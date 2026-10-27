import numpy as np
from glob import glob
import datetime
import matplotlib.pyplot as plt
import math
import sys
import re  #regular expression operations
import scipy.stats as ss
import scipy.interpolate as interp
import scipy.io as io
import os
from matplotlib.colors import LogNorm

def read_info_file(file):

    AllLines = [] #存所有的dat内容
    with open(file,'r') as f:
        for line in f:
            AllLines.append(line)
        f.close()
        
    nLines = len(AllLines)
    IsDone = 0
    iLine = 0

    Vars = {} #type=dict=key/value hash table structure: {key1:value1, key2:value2...}
    
    while (not IsDone):

        Line = AllLines[iLine]
        
        m = re.match(r'#GRID',Line) #run_info.dat里本来是单独一行#Grid 实际代码读出来是#Grid\n r'#Grid'只读Grid 不读new line的\n
        if m: #读取#Grid的内容
            iLine = iLine + 1
            s = AllLines[iLine].split() #s=124 nFieldLines iLine=2
            nFieldLines = int(s[0].strip())

            iLine = iLine + 1
            s = AllLines[iLine].split()
            nz0 = int(s[0].strip()) #nz0=nPointsAlongFieldLine0 iLine=3

            iLine = iLine + 1
            s = AllLines[iLine].split()
            nPoints = int(s[0].strip()) #nPoints=nPointsAlongFieldLine iLine=4

            iLine = iLine + 1
            s = AllLines[iLine].split()
            nLons = int(s[0].strip()) #nLons=# of longtitude  iLine=5

            iLine = iLine + 1
            s = AllLines[iLine].split()
            nProcs = int(s[0].strip()) #nProcs=# of workers   iLine=6

        m = re.match(r'#VARIABLE',Line)
        if m:
            iLine = iLine + 1 
            while (not IsDone):  #iLine=9
                cols = AllLines[iLine].split(',')
                Vars[cols[0]] = cols[1]+' ('+cols[2].strip()+')' #Vars='teu.dat': 'Electron Temperature (K)'
                iLine = iLine + 1
                if (iLine >= nLines):
                    IsDone = 1
            
        iLine = iLine + 1
        
        if (iLine >= nLines):
            IsDone = 1

    Grid = {'nFL' : nFieldLines,
            'Points' : nPoints,
            'nLons' : (nLons-2) * nProcs} #【这个nLongitude和nWorkers的关系是为什么啊？】-2=ghost cells

    return(Grid,Vars)




def read_time_file(file, basetime):
    # 【只是为了读time.dat有多少行 为了下面while loop加上basetime？】
    AllLines = []
    with open(file,'r') as f:
        for line in f:
            AllLines.append(line) #time.dat中全部行变成list
        f.close()
    
    nLines = len(AllLines)
    
#     print(AllLines)
    IsDone = 0
    iLine = 0
    times = []
     #【time.dat中的每一行time加上basetime的2017/9/6】
    while (not IsDone):
        Line = AllLines[iLine] #time.dat中每一行=list
        cols = AllLines[iLine].split() #time.dat每一行的list变成str
        times.append(basetime + datetime.timedelta(hours=float(cols[4]))) #time.dat中最后一个数字4.2251EE+01=42.251hour=1day,65703s
        # basetime=[2017,9,6,0,0]   times=[2017yr,9mon,7day,18hr,15min,3s,600000micros]
        iLine = iLine + 1
        if (iLine >= nLines): 
            IsDone = 1

    return times


def read_coordinates(ddir,filelon,filelat,filealt, nf, nz, nlt):

    ifile1 = ddir+filelon
    ifile2 = ddir+filelat
    ifile3 = ddir+filealt

    # lon,lat,alt
    with open(ifile1,'rb') as f: #【‘rb’是binary file】
        f.seek(4) #read from the 4th characters
        tmp = np.fromfile(f,dtype='float64',count=nf*nz*nlt)
    lon0=np.reshape(tmp,(nlt,nf,nz)) #lon0=ndarray(720,124,304)

    with open(ifile2,'rb') as f:
        f.seek(4)
        tmp = np.fromfile(f,dtype='float64',count=nf*nz*nlt) # or float 64
    lat0=np.reshape(tmp,(nlt,nf,nz))

    with open(ifile3,'rb') as f:
        f.seek(4)
        tmp = np.fromfile(f,dtype='float64',count=nf*nz*nlt)
    alt0=np.reshape(tmp,(nlt,nf,nz))

#     with open(ifile4,'rb') as f:
#         for k in np.arange(nt):
#             print('Reading time :',k)
#             # Fortran length:
#             dd = np.fromfile(f,dtype=np.int32,count=1)
#             # Read data:
#             dene0 = np.fromfile(f,dtype=np.double,count=nlt*nf*nz)
#             # Fortran length:
#             dd = np.fromfile(f,dtype=np.int32,count=1)
#             dene[:,:,:,k] = np.reshape(dene0,(nlt,nf,nz))

    return lon0,lat0,alt0


def find_alt_index(alt,sel_alt,nz): #找single altitude的alt_grid[a] --> sel_alt
    
# Find index of a certain altitude in SAMI3 coordinate 

    if(sel_alt < alt[0]): #alt先84变大到90再变到84
        k = np.nan
        l = np.nan
#         print('--smaller than alt(0)')
    elif (sel_alt > np.max(alt)):
        k = np.nan
        l = np.nan
#         print('--larger than maxval(alt)')
    else:
        mm = 1 #alt从头往后找
        while (sel_alt >= alt[mm]):
            mm = mm+1
        k = mm-1
            #print('South alt/k,k+1',alt[k],alt[k+1])

        mm = nz-1 #alt从尾往前找
        while (sel_alt >= alt[mm]):
            mm = mm-1
        l = mm+1
            #print('North alt/k,k-1',alt[k],alt[k-1])

    return k,l


def interp_alt(x1,x2,y1,y2,x):
    return (y2*(x-x1)+y1*(x2-x))/(x2-x1)



def contour_ialt(alt,data,blonst_g,blatst_g,altst_g,nf, nz, nlt):

    # Find grids (2D cut) at a single altitude
    
    ndata = np.ones((nlt,nf*2))*np.nan #nlt=6 nf=124
    lons = np.ones((nlt,nf*2))*np.nan
    lats = np.ones((nlt,nf*2))*np.nan

    for iLon in np.arange(nlt):
        for iFl in np.arange(nf):

            alts = altst_g[iLon,iFl,:]
            alt_index1,alt_index2 = find_alt_index(alts,alt,nz) # Find north and south ends

            if (alt_index1 != alt_index1):
                continue
            lontmp = blonst_g[iLon,iFl,alt_index1]
            lons[iLon,iFl] = lontmp
            lattmp = blatst_g[iLon,iFl,alt_index1]
            lats[iLon,iFl] = lattmp
            #interp_alt(x1,x2,y1,y2,x): (y2*(x-x1)+y1*(x2-x))/(x2-x1)
            ndata[iLon,iFl] = interp_alt(alts[alt_index1], alts[alt_index1+1], 
                                         data[iLon,iFl,alt_index1], data[iLon,iFl,alt_index1+1], alt)

            lontmp = blonst_g[iLon,iFl,alt_index2]
            lons[iLon,-iFl] = lontmp
            lattmp = blatst_g[iLon,iFl,alt_index2]
            lats[iLon,-iFl] = lattmp
            ndata[iLon,-iFl] = interp_alt(alts[alt_index2], alts[alt_index2-1], 
                                         data[iLon,iFl,alt_index2], data[iLon,iFl,alt_index2-1], alt)
            #print((altst_g[iLon,iFl,alt_index1]),
            #        (altst_g[iLon,iFl,alt_index2]),
            #        (blatst_g[iLon,iFl,alt_index1]),
            #        (blatst_g[iLon,iFl,alt_index2]))
    return ndata,lons,lats


def read_one_file_one_time(file,iTime,Grid, nf, nz, nlt):

    #iSize = nf*nz*nlt
    iSize = Grid["nFL"] * Grid["Points"] * Grid["nLons"]

    # Assume 8 byte floats, with 4 byte length header and footer:
    iPos = iTime * (4 + iSize * 8 + 4)

    print('Reading File : '+file, iPos)
    
    with open(file,'rb') as f:
        f.seek(iPos)
        dd = np.fromfile(f,dtype=np.int32,count=1)
        tmp = np.fromfile(f,dtype='float64',count=iSize) #每次读nf*nz*nlt grid的size
    data=np.reshape(tmp,(nlt,nf,nz))
    f.close
    
    return data



def calc_tec(data,blonst_g,blatst_g,altst_g, nf, nz, nlt, altindex, lonbins=90, latbins=90):
    
    tec=np.zeros((latbins,lonbins))
    #bAlt=each longitude(720) has 124 field lines, with each field line has nz(304) points
    alt_grid=np.max(altst_g, axis=2) #altst_g size=(96nlt,124nf,304nz)  alt_grid=(nlt,nf) at highest altitude
    alt_grid=alt_grid[0]-1 #size=(nf) 1st row(longitude)'s alt_grid
    alt_grid=alt_grid[:altindex] 
    #【为什么只要前一百的highest altitude grid ~15000 altitude？？-->靠近pole的fieldline很高～open field line 
    
    ne=np.zeros((latbins, lonbins, alt_grid.shape[-1]-1)) #ne(90,720,100-1alt_grid)
    
    # Find grids at a single altitude
    for k, a in enumerate(alt_grid[1:]): #loop each single altitude grid
        ndata = np.ones((nlt,nf*2))*np.nan #(720,248)
        lons = np.ones((nlt,nf*2))*np.nan
        lats = np.ones((nlt,nf*2))*np.nan

        for iLon in np.arange(nlt):
            for iFl in np.arange(nf):

                alts = altst_g[iLon,iFl,:] #每个longitude和fieldline point上的304nz=每个磁感线上的点
                alt_index1,alt_index2 = find_alt_index(alts,a,nz) #【Find north and south ends】: a在loop through之前的100alt_grid
                
                if (alt_index1 != alt_index1): #找到南北不同ends就继续往下，没找到(index相同)lons/lats/ndata就是nan
                    continue
                lontmp = blonst_g[iLon,iFl,alt_index1] #lons/lats/ndata前半段south alt
                lons[iLon,iFl] = lontmp
                lattmp = blatst_g[iLon,iFl,alt_index1]
                lats[iLon,iFl] = lattmp
                #interp_alt(x1,x2,y1,y2,x): (y2*(x-x1)+y1*(x2-x))/(x2-x1)
                ndata[iLon,iFl] = interp_alt(alts[alt_index1], alts[alt_index1+1], 
                                         data[iLon,iFl,alt_index1], data[iLon,iFl,alt_index1+1], a)

                lontmp = blonst_g[iLon,iFl,alt_index2] #lons/lats/ndata后半段north alt
                lons[iLon,-iFl] = lontmp
                lattmp = blatst_g[iLon,iFl,alt_index2]
                lats[iLon,-iFl] = lattmp
                ndata[iLon,-iFl] = interp_alt(alts[alt_index2], alts[alt_index2-1], 
                                         data[iLon,iFl,alt_index2], data[iLon,iFl,alt_index2-1], a)
               
        lp = (np.isnan(ndata))
        lp = np.invert(lp) #ndata不是nan的 (nlt,2*nf) lons=(nlt,2*nf) lons[lp]只留下不是nan的然后2darray变成1d--2nf dimension并到每一行
        #binned_statistic_2d:lo=x_edge la=y_edge
        st, lo, la, c=ss.binned_statistic_2d(lons[lp], lats[lp], \
                    ndata[lp], bins=[np.linspace(0,360,lonbins+1), \
                                     np.linspace(-85,85,latbins+1)])
        st=st.transpose() #st(720,90)=(721latbins,91lonbins)=Statistics value for each bin
        
        for i in range(lonbins): #loop st每一列把nan的值interp1d
            lp = (np.isnan(st[:,i]))
            lp = np.invert(lp) #loop每个latbins找不是nan的
            dl=170/latbins
            x=np.linspace(-85,85,latbins+1)
            x=x[:-1]+dl/2 #【为什么要这么定义interp的x？---binning center】
            f = interp.interp1d(x[lp], st[lp,i], fill_value='extrapolate') 
            xnew = np.linspace(-85,85,latbins+1)
            st[:,i] = f(xnew[:-1])
        
        ne[:,:,k]=st
    
    cnt = np.zeros((latbins, lonbins))
    for i in range(latbins): #开始积分TEC了
        for j in range(lonbins):
            for k in range(alt_grid.shape[-1]-1):
                if not np.isnan(ne[i,j,k]):
                    tec[i,j] = tec[i,j] + ne[i,j,k] * 1e6 * (alt_grid[k+1] - alt_grid[k-int(cnt[i,j])]) * 1e3 * 1e-16 
                    # 1e6:/cm^3->/m^3, 1e3:km->m 1e-16:TECU[el/m^2]
                    cnt[i,j] = 0
                else: #ne[i,j,k]不是nan的话 alt积分时就是alt[k+1]和[k]的grid, 如果ne[i,j,k]是nan的话 alt积分grid数量[k-cnt]倒回去几个grid
                    cnt[i,j] = cnt[i,j] + 1
        
    return la[:-1], lo[:-1], ne, tec



def calc_parameter(data,blonst_g,blatst_g,altst_g, nf, nz, nlt,lonbins=90, latbins=90):
    
    #bAlt=each longitude(720) has 124 field lines, with each field line has nz(304) points
    alt_grid=np.max(altst_g, axis=2) #altst_g size=(96nlt,124nf,304nz)  alt_grid=(nlt,nf) at highest altitude
    alt_grid=alt_grid[0]-1 #size=(nf) 1st row(longitude)'s alt_grid
    alt_grid=alt_grid[:100] 
    #【为什么只要前一百的highest altitude grid ~15000 altitude？？-->靠近pole的fieldline很高～open field line 
    
    out=np.zeros((latbins, lonbins, alt_grid.shape[-1]-1)) #ne(90,720,100-1alt_grid)
    
    # Find grids at a single altitude
    for k, a in enumerate(alt_grid[1:]): #loop each single altitude grid
        ndata = np.ones((nlt,nf*2))*np.nan #(720,248)
        lons = np.ones((nlt,nf*2))*np.nan
        lats = np.ones((nlt,nf*2))*np.nan

        for iLon in np.arange(nlt):
            for iFl in np.arange(nf):

                alts = altst_g[iLon,iFl,:] #每个longitude和fieldline point上的304nz=每个磁感线上的点
                alt_index1,alt_index2 = find_alt_index(alts,a,nz) #【Find north and south ends】: a在loop through之前的100alt_grid
                
                if (alt_index1 != alt_index1): #找到南北不同ends就继续往下，没找到(index相同)lons/lats/ndata就是nan
                    continue
                lontmp = blonst_g[iLon,iFl,alt_index1] #lons/lats/ndata前半段south alt
                lons[iLon,iFl] = lontmp
                lattmp = blatst_g[iLon,iFl,alt_index1]
                lats[iLon,iFl] = lattmp
                #interp_alt(x1,x2,y1,y2,x): (y2*(x-x1)+y1*(x2-x))/(x2-x1)
                ndata[iLon,iFl] = interp_alt(alts[alt_index1], alts[alt_index1+1], 
                                         data[iLon,iFl,alt_index1], data[iLon,iFl,alt_index1+1], a)

                lontmp = blonst_g[iLon,iFl,alt_index2] #lons/lats/ndata后半段north alt
                lons[iLon,-iFl] = lontmp
                lattmp = blatst_g[iLon,iFl,alt_index2]
                lats[iLon,-iFl] = lattmp
                ndata[iLon,-iFl] = interp_alt(alts[alt_index2], alts[alt_index2-1], 
                                         data[iLon,iFl,alt_index2], data[iLon,iFl,alt_index2-1], a)
               
        lp = (np.isnan(ndata))
        lp = np.invert(lp) #ndata不是nan的 (nlt,2*nf) lons=(nlt,2*nf) lons[lp]只留下不是nan的然后2darray变成1d--2nf dimension并到每一行
        #binned_statistic_2d:lo=x_edge la=y_edge
        st, lo, la, c=ss.binned_statistic_2d(lons[lp], lats[lp], \
                    ndata[lp], bins=[np.linspace(0,360,lonbins+1), \
                                     np.linspace(-85,85,latbins+1)])
        st=st.transpose() #st(720,90)=(721latbins,91lonbins)=Statistics value for each bin
        
        for i in range(lonbins): #loop st每一列把nan的值interp1d
            lp = (np.isnan(st[:,i]))
            lp = np.invert(lp) #loop每个latbins找不是nan的
            dl=170/latbins
            x=np.linspace(-85,85,latbins+1)
            x=x[:-1]+dl/2 #【为什么要这么定义interp的x？】
            f = interp.interp1d(x[lp], st[lp,i], fill_value='extrapolate') 
            xnew = np.linspace(-85,85,latbins+1)
            st[:,i] = f(xnew[:-1])
        
        out[:,:,k]=st
    
        
    return la[:-1], lo[:-1], out, alt_grid[:-1]


def calc_on2(oden, n2den, height):
    nlats = oden.shape[0]
    nlons = oden.shape[1]
    nalts = height.shape[0]
    MaxValN2 = 1.0e21 # 1e21/m^2
    o = np.zeros([nlats, nlons])
    n2 = np.zeros([nlats, nlons])
    for i in range(nlats):
        for j in range(nlons):
            for k in range(nalts-2, 1, -1):
                dAlt = (height[k] - height[k-1]) * 1e3 # baltu.dat [km]->[m]
                oMid = (oden[i,j,k] + oden[i,j,k-1]) / 2.0 * 1e6 # density [/cm^3] -> [/m^3] 
                n2Mid = (n2den[i,j,k] + n2den[i,j,k-1]) / 2.0 * 1e6
                if n2[i,j] + n2Mid*dAlt < MaxValN2:
                    n2[i,j,] = n2[i,j] + n2Mid * dAlt
                    o[i,j] = o[i,j] + oMid * dAlt
                else:
                    dAlt = (MaxValN2 - n2[i,j]) / n2Mid
                    n2[i,j] = n2[i,j] + n2Mid * dAlt
                    o[i,j]  =  o[i,j] +  oMid * dAlt
                    break
    on2ratio = o/n2
    return on2ratio



def read_growth(ddir,file,nlt,nf,nt):

    data = np.zeros((nlt,nf,nt))

    with open(ddir+file,'rb') as f:
        f.seek(4)
        for k in range(nt):
            print('Reading time :',k)
#             # Fortran length:
            # Read data:
            d0 = np.fromfile(f,dtype=np.float64,count=nlt*nf)
            # Fortran length:
#             dd = np.fromfile(f,dtype=np.int32,count=1)
            data[:,:,k] = np.reshape(d0,(nlt,nf))
        f.seek(4)
    return data

# grti=read_growth(EPB_dir,'gammaEE_cd.dat') #grti=[180Lon,124nf,96times]


def read_EPB_file_one_time(file,iTime,Grid,nf,nlt):

    iSize = Grid["nFL"]  * Grid["nLons"] #124*180

    # Assume 8 byte floats, with 4 byte length header and footer:
    iPos = iTime * (4 + iSize * 8 + 4)

    print('Reading File : '+file, iPos)
    
    with open(file,'rb') as f:
        f.seek(iPos)
        dd = np.fromfile(f,dtype=np.int32,count=1)
        tmp = np.fromfile(f,dtype='float64',count=iSize) #每次读nf*nz*nlt grid的size
    data=np.reshape(tmp,(nlt,nf))
    f.close
    
    return data