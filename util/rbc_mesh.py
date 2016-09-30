# Generate RBC triangular mesh by optimizing a randomly initialized one

from numpy import *
from scipy import *
from scipy.spatial import *
from time import time
import sys
from getopt import getopt

opts, args = getopt( sys.argv[1:], 's:v:a:e:' )
opts       = dict( opts )
n_step     = int(opts['-s']) if '-s' in opts else 100
n_vert     = int(opts['-v']) if '-v' in opts else 500
seed       = int(opts['-e']) if '-e' in opts else 0
arx        = opts['-a'] if '-a' in opts else 'v%ds%d' % ( n_vert, n_step ) 
k          = 0.004 * n_vert ** (1/3.0);
eps        = 1e-7
cut        = 4 / n_vert ** 0.5

print 'Generating particles by Poisson disk'
t0    = time()
random.seed( seed )
vert  = random.randn( 1, 3 )
vert /= linalg.norm(vert)
cache = vert.copy()
for i in range( 0, n_vert - 1 ):
    if i % floor( n_vert / 10 ) == 0:
        sys.stdout.write( '%d...' % i )
        sys.stdout.flush()
    tree = cKDTree( vert )
    while True:
        nv   = random.randn( 1, 3 )
        nv  /= linalg.norm( nv )
        nv  *= 1 + eps * ( random.rand() - 0.5 )
        d, i = tree.query( nv )
        if d[0] >= cut * 0.65:
            break
    vert = append( vert, nv, axis = 0 )
sys.stdout.write('Done\n') 

# Reorder to improve locality
nbin = floor( 2 / cut ) / 2
binv = nbin / 2
bid  = sum( floor( ( vert + 1 ) * binv ) * array( [ 1, nbin, nbin * nbin ] )[newaxis,:], axis = 1 )
vert = vert[ argsort( bid ), : ]  

print '%d vertices initialization => %.2f seconds' % ( n_vert, time() - t0 )

t0 = time()

# Iteratively optimize the mesh
for step in range(n_step):
    # Repulse
    if step % floor(n_step/10) == 0:
        sys.stdout.write( '%d...' % step )
        sys.stdout.flush()
    dtri   = Delaunay( vert )
    face   = dtri.convex_hull.copy()
    edge   = concatenate( ( face[ :, [0,1] ], face[ :, [0,2] ], face[ :, [1,2] ] ) )
    edge   = sort( edge, axis = 1 )
    e2i    = unique( edge[:,0] * n_vert + edge[:,1] )
    bond   = vstack( ( e2i / n_vert, e2i % n_vert ) ).transpose()
    dx     = vert[ bond[:,0], : ] - vert[ bond[:,1], : ]
    dr     = linalg.norm( dx, axis = 1 )
    nr     = dx / tile( dr, (3,1) ).transpose()
    r0     = mean( dr )
    force  = nr * tile( -k * ( dr - r0 ), (3,1) ).transpose()
    for b, f in zip( bond, force ):
        vert[ b[0], : ] += f
        vert[ b[1], : ] -= f
    # Constrain
    dist   = linalg.norm( vert, axis = 1 )
    alpha  = 0.1
    centri = vert * tile( alpha / dist - alpha, ( 3, 1 ) ).transpose()
    vert   = vert + centri

# Post-processing: project to unit sphere
vert /= linalg.norm( vert, axis = 1 )[ :, newaxis ]

sys.stdout.write('Done\n')
print '%d vertices, %d iterations => %.2f seconds' % ( n_vert, n_step, time() - t0 )

savetxt( '%s.vert.txt' % arx, vert, fmt='% 2.8e', delimiter=' ' )
savetxt( '%s.face.txt' % arx, face, fmt='%-7d', delimiter=' ' )
savetxt( '%s.bond.txt' % arx, bond, fmt='%-7d', delimiter=' ' )

