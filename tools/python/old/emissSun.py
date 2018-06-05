import read as read
import numpy as np
import os

#============================================================================
def emiss(data,ion,Teunit,nunit,Lunit,othick=0):
    '''get optically thin emission in specific wave length'''
    print '=== Getting emission ==='
# peel data:
    pointdata=data.get('pointdata')
# points:
    rho=read.extract(pointdata,'rho',attribute_mode='point')
    Te=read.extract(pointdata,'Te',attribute_mode='point')
    rho=rho*nunit
    Te=Te*Teunit
# read in emissivity table
    tabledic=loadGtable(ion,filenm='aia')
    logT=tabledic['logt']
    logn=tabledic['n_e_lg']
    gmat=tabledic['goft_mat']
    ce=findintable(ion,Te,rho,logT,logn,gmat)
# emission flux in unit of DN s^-1 per unit length
    emiss = rho**2 * ce * Lunit
    if othick != 0:
        emiss = np.where(rho>2.e10,0.,emiss)
# mask out the wind zone:
    print '=== Done with emission! ==='
    return {'emiss': emiss, 'points': data.get('points')}
#============================================================================
def integral(data,val,unit,los):
    '''get optically thin emission in specific wave length'''
    print '=== Getting emission ==='
# peel data:
    pointdata=data.get('pointdata')
# points:
    emiss=read.extract(pointdata,val,attribute_mode='point')*unit/los
# mask out the wind zone:
    print '=== Done with integral value ==='
    return {'emiss': emiss, 'points': data.get('points')}
#============================================================================
def viewangle(emissdic,theta,phiy,phi,delta=0,nu=1):
    '''set view angle for line of sight'''
    emissdic['phi']=phi
    emissdic['phiy']=phiy
    emissdic['theta']=theta
    emissdic['delta']=delta
    emissdic['nu']=nu
    emissdic['alpha']=0
    emissdic['recipe']=0
    return emissdic
#============================================================================
def loadGtable(ion,w0=None,filename=None,filenm=None,abund='co'):
    '''reads response function (n_e, T_e) 2D tables, requires ion'''

    if filename == None:
    # Construct filename:
        dirgot = os.environ['AMRVAC_DIR']+'/tools/python/sdoaia/'
        if filenm != None:
            if abund in ['co','Co','corona','Corona']:
                filegot = ''.join(['goft_table_',filenm,str(ion),'_abco_extro.dat'])
            else:
                filegot = ''.join(['goft_table_',filenm,str(ion),'_abph_extro.dat'])
        else:
            w0nm = w0.zfill(4)
            filegot = ''.join(['goft_table_',str(ion),'_',w0nm,'.dat'])

        filename = dirgot+filegot

    print'load G(n,T) table: ',filename   

    f = open(filename, 'r')

    ion   = float(f.readline())
    w00   = float(f.readline())
    watom = float(f.readline())
    (numn,numt) = f.readline().split()
    numn=int(float(numn));numt=int(numt)
    
    logt  = readArray(f,numt)

    goft_mat = np.empty((numn,numt))
    n_e_lg   = np.empty(numn)

    i=0
    while True:
        line     = f.readline()
        if not line: break
        n_e_0    = float(line)
        g_t      = readArray(f,numt)
        goft_mat[i,:] = g_t
        n_e_lg[i] = n_e_0
        i=i+1

    f.close()
    return {'w0': w0, 'ion': ion, 'n_e_lg': n_e_lg, 'logt': logt, 
            'goft_mat': goft_mat,'watom': watom}
#============================================================================
def readArray(f,num):
    '''read array from data file, require filehandle and length of the 
    array'''
    arr=np.empty(num)
    i=0
    while i < num:
        tmp=f.readline().split()
        for elem in tmp:
            arr[i] = float(elem)
            i=i+1

    return arr
#============================================================================
def findintable(ion,Te,ne,logT,logn,gmat):
    '''Find emissivity in perticular wave length or filter, given density and
       temperature'''
    Tel=np.log10(Te)
    nel=np.log10(ne)
    LT=findL(Tel,logT)
    Ln=findL(nel,logn)
    #bilinear interpolation in gmat(logT,logn) 2D table
    dlog=logT[LT[1]]-logT[LT[0]]
    dlog[dlog==0.]=1.
    frac1=(logT[LT[1]]-Tel)/dlog
    frac2=(Tel-logT[LT[0]])/dlog
    frac1[LT[1]==LT[0]]=0.5
    frac2[LT[1]==LT[0]]=0.5
    #linear interpolation in logT
    gr1=frac1*gmat[Ln[0],LT[0]]+frac2*gmat[Ln[0],LT[1]] 
    gr2=frac1*gmat[Ln[1],LT[0]]+frac2*gmat[Ln[1],LT[1]]
    dlog=logn[Ln[1]]-logn[Ln[0]]
    dlog[dlog==0.]=1.
    frac1=(logn[Ln[1]]-nel)/dlog
    frac2=(nel-logn[Ln[0]])/dlog
    frac1[Ln[1]==Ln[0]]=0.5
    frac2[Ln[1]==Ln[0]]=0.5
    corona_emiss=frac1*gr1+frac2*gr2 #linear interpolation in logn
    return corona_emiss
#============================================================================
def findL(parr,aseq):
    '''Find the location of a value in an ascending equi-distant sequence, 
    return list of two indices'''
    lenseq=len(aseq)
    dx=aseq[1]-aseq[0]
    jl=np.int32((parr-aseq[0])/dx)
    jh=jl+1
    jh[jl<0]=0
    jl[jl<0]=0
    jl[jh>lenseq-1]=lenseq-1
    jh[jh>lenseq-1]=lenseq-1
    return [jl,jh]
#============================================================================
def aiact(ion):
    '''reads AIA RGB color tables, requires ion'''
    # Construct filename:
    filename=os.environ['AMRVAC_DIR']+'/tools/python/sdoaia/aiacolartable'+str(ion)+'.txt'

    print'file name of color table: ',filename   

    f = open(filename, 'r')
    
    rr = readArray(f,256)
    gg = readArray(f,256)
    bb = readArray(f,256)
 
    f.close()

    rm= np.empty((256,3))
    rm[:,0]=np.arange(256.)/255.
    rm[:,1]=rr
    rm[:,2]=rr
    gm= np.empty((256,3))
    gm[:,0]=np.arange(256.)/255.
    gm[:,1]=gg
    gm[:,2]=gg
    bm= np.empty((256,3))
    bm[:,0]=np.arange(256.)/255.
    bm[:,1]=bb
    bm[:,2]=bb

    return {'red': tuple(map(tuple,rm)),'green': tuple(map(tuple,gm)),
            'blue': tuple(map(tuple,bm))}
#============================================================================
