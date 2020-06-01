# raytrace
Shanti's Raytracer: fast, precise, and accurate optical ray propagation

Start with "examples.m" and read the instructions in "raytrace.m".

To see how one might instantiate a giant segmented telescope, and compute the effect on the pupil-referred wavefront if any of the optical elements were to move a microradian or micrometer in any direction, run "exampleHabEx.m".

The [MATLAB](matlab/) code was originally developed at JPL, and is released as open-source under an [Apache 2.0 license](LICENSE). US Government sponsorship is acknowledged.

## Tutorial
You will create an object for a ray bundle, calculate how they are affected by a sequence of surface interfaces. The ray propagation equations are coordinate-free, which means that every position and direction needs to be expressed in 3D.

Start with a source ray bundle with one of the built-in functions, **sourceColumn** for collimated beams, **sourcePoint** for spherical wavefronts, and **sourceFan** for a 1-dimensional section of a spherical wavefront.

    function sourceColumn(position, direction, x, radius, Nr, RefIndex)
   
Let's decide that _+z_ will be the optical axis, and define a collimated source centered 100 mm in the _+y_ direction. It will have 9 rays across the radius of 50 mm. It's good to tell it that the dimensions are in mm.

    source = sourceColumn([0, 100, 0],[0, 0, 1],[1, 0 0], 50, 9)
    source.units = 'mm';

Next, place a concave mirror some distance away, pointing in the _-z_ direction. We'll give it a radius of curvature of 500 mm, which means that the curvature, **cuy**, is 1/500.

    m1.position = [0 0 50];
    m1.direction = [0 0 -1];
    m1.cuy = 1/500;
    m1.type = 'reflect';
    
We'll put another surface, a focal plane, at prime focus. Let's tell it that we want wavefront errors in nanometers.

    fp.position = [0 0 -450];
    fp.direction = [0 0 1]; % direction doesn't actually matter for flat surfaces
    fp.displayUnits = 'nm';

Next, use the **raytrace** function to map the source through the surfaces.

    trace = raytrace(source,{m1,fp});

There are also some handy functions for displaying what happens. Let's plot the geometry, which you can spin around in 3D, and then the wavefront error map, a spot diagram, and a point-spread function for a 633 nm light source. This is why we set the **displayUnits** property for the focal plane -- so that wavefront error is printed in nm.

    figure(1);clf;
    plotSurfaces(trace); axis equal; sideview;
    hold on;plotRays(trace,'b');hold off;

    figure(2);clf;
    subplot(1,3,1); [pupil,mask,rmswfe] = displayWFE(trace{end});axis image;
    subplot(1,3,2); plotSpot(trace{end});axis image;
    subplot(1,3,3); 
    psf = pupilPSF(pupil, mask, .000633, 1, 1/2/m1.cuy, 19, source.samplingDistance, .010);
    imagesc(psf); axis image;

![Ray side view 1](/tutorial/matlabtutorial1.png)

![Ray side view 1](/tutorial/matlabtutorial2.png)

That wavefront looks awful, of course, and the PSF is meaningless. So let's change the mirror from a sphere to a parabola:

    m1.K = -1;
    
Now, run it all again

![Ray side view 1](/tutorial/matlabtutorial3.png)

![Ray side view 1](/tutorial/matlabtutorial4.png)

That looks better! Note the PSF has concentrated all the energy into the central pixel.

Next, what happens if the primary mirror is out of place? We'll create a perturbation, a 10 microradian rotation about the _x_ axis.
    
    pert.axis = normr([1 0 0]);
    pert.angle = 1e-5;
    
Raytrace the same source through the perturbed mirror, and what happens? There's a 70 nm wavefront error that looks like tilt. Don't believe too much in the PSF, though -- the PSF calculation is quite susceptible to computation artifacts.
 
    trace = raytrace(source,{perturb(m1,pert),fp});

![Ray side view 1](/tutorial/matlabtutorial6.png)

## References

[C J Mitchell 1981 J. Opt. 12 301](https://iopscience.iop.org/article/10.1088/0150-536X/12/5/003) Generalized ray-tracing for diffraction gratings of arbitrary form

[Richard J. Mathar Serb. Astr. J. 179 (2009) 107-120](https://arxiv.org/abs/0809.2368) Zernike Basis to Cartesian Transformations
