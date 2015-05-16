import struct


### Provide extra information to read the .dat file:
# Grid sizes at first AMR level:
nxlone = [16,40]
# Physical domain:
xprobmin = [0.0,0.0]
xprobmax = [0.0315421,0.03]
filename = './KH0080.dat'




def main():

    dat = open(filename,'r')
    dat.seek(0,2)       # goto EOF
    offset = -36        # 36 = 7*size_int + size_double = 7*4 + 8
    dat.seek(offset,1)  # go back 36 bytes



    ### Read 7 footer ints + 1 double
    int_in = dat.read(4)
    nleafs = struct.unpack('i',int_in)[0]    # number of active tree leafs nleafs (= #blocks)
    int_in = dat.read(4)
    levmaxini = struct.unpack('i',int_in)[0] # maximal refinement level present levmax
    int_in = dat.read(4)
    ndimini = struct.unpack('i',int_in)[0]   # dimensionality NDIM
    int_in = dat.read(4)
    ndirini = struct.unpack('i',int_in)[0]   # number of vector components NDIR
    int_in = dat.read(4)
    nwini = struct.unpack('i',int_in)[0]     # number of variables nw
    int_in = dat.read(4)
    neqparini = struct.unpack('i',int_in)[0] # number of equation-specific variables neqpar+nspecialpar
    int_in = dat.read(4)
    it = struct.unpack('i',int_in)[0]        # integer time counter it
    flt_in = dat.read(8)
    t = struct.unpack('d',flt_in)[0]         # global time t




    ### read block size in each dimension and all eqpars
    offset = offset-int(ndimini*4 + neqparini*8) # int size = 4, double = 8
    dat.seek(offset,1)
    nx = []                     # Block size
    neqpar = []
    for i in range(ndimini):
        int_in = dat.read(4)
        nx.append(struct.unpack('i',int_in)[0])

    for i in range(neqparini):
        flt_in = dat.read(8)
        neqpar.append(struct.unpack('d',flt_in)[0])
        
    cells = 1
    for i in range(ndimini):
        cells *= nx[i]           # number of cells per block
        
    size_block = cells*nwini*8   # block size in bytes
    print neqpar





    ### Read forest
    dat.seek(nleafs*size_block,0)   # go to end of block stuctures, start of the forest
    ng = [1,1,1]
    for i in range(ndimini):
        ng[i] = nxlone[i]/nx[i]     # number of blocks in each direction

    block_info = []    
    igrid = 0
    refine = 0
    forest = []
    level = 1
    for k in range(1,ng[2]+1):
        for j in range(1,ng[1]+1):
            for i in range(1,ng[0]+1):
                (igrid,refine) = read_node(dat,forest,igrid,refine,ndimini,level,block_info,[i,j,k])


                
    
    # open file for output
    outname = ''
    for i in range(len(filename)-3):
        outname += filename[i]
    outname += 'csv'    
    outfile = open(outname,'w')
    outfile.write('X,Y,Z,rho,m1,m2,m3,e,rhod1,rhod2,m1d1,m1d2,m2d1,m2d2,m3d1\n') 
    #outfile.write(' ' + str(cells*nleafs) + '\n')     #total number of cells        
    
    


    # Calculate physical sizes of blocks and cells on level one
    dx  = []    # physical block size in each direction (on the first level)
    dxc = []    # physical cell size in each direction (on the first level) 
    for i in range(ndimini):
        dx.append((xprobmax[i] - xprobmin[i])/ng[i])
        dxc.append((xprobmax[i] - xprobmin[i])/nxlone[i])
     


    
    ### Start reading data blocks
    print 'Start reading data blocks...'
    dat.seek(0,0)  # go to start of file
    for i in range(nleafs):
        lb = []         # coordinates of bottom left corner of block 
        dxc_block = []  # local cell size
        for n in range(ndimini):
            lb.append(xprobmin[n]+(block_info[i][1][n] - 1)*dx[n]/(2**(block_info[i][0]-1)))
            dxc_block.append(dxc[n]/(2**(block_info[i][0]-1)))     
        values = [ [] for i in range(nwini) ]
        for nw in range(nwini):
            for c in range(cells):
                flt_in = dat.read(8)
                dbl = struct.unpack('d',flt_in)[0]
                values[nw].append(dbl)
        ### Output the block to ascii file
        for c in range(cells):
            xyz = calc_coord(c,nx,lb,dxc_block)
            for d in range(3):
                outfile.write(('%e' % xyz[d]) + ',' )
            for nw in range(nwini-1):
                outfile.write(('%e' % values[nw][c]) + ',')
            outfile.write(('%e' % values[nwini-1][c]) + '\n')                        



    dat.close()
    outfile.close()




def read_node(dat,forest,igrid,refine,ndimini,level,block_info,indices):
    # Recursive function that checks if the block is a leaf or not. 
    # If leaf, write block info. If not, run recursively for all 
    # possible child nodes

    b = dat.read(4)
    leaf = struct.unpack('i',b)[0]
    if (leaf):
        forest.append(1)
        igrid += 1
        block_info.append([level,indices])
    else:
        forest.append(0)
        refine += 1
        for i in range(2**ndimini):
            child_index = new_index(ndimini,i,indices)
            (igrid,refine) = read_node(dat,forest,igrid,refine,ndimini,level+1,block_info,child_index)
    return igrid,refine
    
    
    
    
def new_index(dims,i,ig):
    # Calculate new grid indices for child nodes if node is refined. 
    # See Keppens (2012) section 3.3 
    if (dims==1):
        new = [2*(ig[0]-1) + 1 + i,1,1]
    if (dims==2):
        new = [2*(ig[0]-1) + 1 + i%2,2*(ig[1]-1) + 1 + i/2,1]
    if (dims==3):
        new = [2*(ig[0]-1) + 1 + i%2,2*(ig[1]-1) + 1 + (i/2)%2,2*(ig[2]-1) + 1 + i/4]
    return new




def calc_coord(c,nx,lb,dxc_block):
    # calculate the cell coordinates from the cell number
    local_ig = [1,1,1]
    local_ig[0] = 1 + (c%nx[0])
    local_ig[1] = 1 + (int(c/nx[0])%nx[0])
    if(len(nx)==3):     #if 3D blocks
        local_ig[2] = 1 + int(c/(nx[0]*nx[1]))
    xyz = [0,0,0]
    for i in range(len(nx)):
        xyz[i] = lb[i] + (local_ig[i] - 0.5)*dxc_block[i]    
    return xyz


if __name__ == "__main__":
    main()
    















