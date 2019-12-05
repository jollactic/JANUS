# JANUS
--- Pyhon-based nanoparticle interface generator ---    
  
This small code can be used to create a set of interface geometries between two nanoparticles given by two xyz-file . The code is based on the quick hull algorithm by Barber et al. [ACM Trans Math Softw  22, 469 (1996)] which allows for, among other things, rapid identification of convex hulls using a set of points as input, in the current case the coordinates of atoms in a nanoparticle model. The procedure to construct interfaces is as follows:  
•	Facets are constructed by merging simplicies of the convex hull whose surface normal form an arcus cosine of 0.99 or larger.  This is done for each particle individually.  
•	The two particles are aligned such that the normal of their largest facet are pointing towards each other and are parallel to the x-axis.  
•	The geometrical mean of the two facets are placed at a distance 'offset_along_x' from each other along the x-axis.  
•	The procedure is repeated with the particles being rotated with respect to one another around the x-axis using angles between 0 and 180 in steps of 15. For each rotation 12 additional interface geometries, with the particles also being displaced in the yz-plane on a circle with radius 'offset_in_yz_plane', using  angles between 0 and 360 in steps of 30, are generated.  
  
This procedure leads to a total of 156 interface geometries for each particle pair. The strucutres are stored in ase-trajectory file called 'JANUS.traj'.
  
To run the code type:  
python3 JANUS.py particle_1.xyz particle_2.xyz offset_along_x offset_in_yz_plane  
  
where 'particle_1.xyz' and 'particle_2.xyz' are two user supplied xyz-files, and 'offset_along_x' and 'offset_in_yz_plane' are to floating numbers controlling the mutual displacements (in Ångströms) of the two particles as described in bullet points above.
