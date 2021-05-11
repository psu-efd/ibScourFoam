# README #
This repository is for the code and cases of ibScourFoam, a solver based on OpenFOAM for 3D scour modeling around structures/objects with complex shape. The scour process is simulated with an improved immersed boundary method which produces smooth wall shear stresses. 

The code is developed with OpenFOAM v5. To use this code, it is required that OpenFOAM has been properly installed. Current code has only been used in Linux. Porting to Windows and Mac OS have not been done, but possible.

This code utilizes some data structure and implementation of an immersed boundary method in the OpenFOAM-extend project. Details can be found in the citations and references of the following paper:

Y. Xu and X. Liu (2021). An immersed boundary method with y+-adaptation wall function for smooth wall shear. International Journal of Numerical Methods in Fluids. https://doi.org/10.1002/fld.4960

## Example simulation cases ##

This repository comes with some examples showing how to use the code and demonstrating its capability to deal with scour around complex structures. Animations of the simulated scour process can be viewed at our group's YouTube channel:

- The classic case of scour around vertical cylinder: [YouTube link](https://www.youtube.com/watch?v=JIKRwKyth5s&list=PLPt0QqiTKrmv6r52cZyp9l11i6auKGSMY&index=4).

- Scour around a short cylinder on a supporting cradle: [YouTube link](https://www.youtube.com/watch?v=tf6YIUMZle0&list=PLPt0QqiTKrmv6r52cZyp9l11i6auKGSMY&index=5)

- Scour around an object with complex shape on a supporting cradle: [YouTube link](https://www.youtube.com/watch?v=VU2N0s7i3nI&list=PLPt0QqiTKrmv6r52cZyp9l11i6auKGSMY&index=6)

If you use our code to simulation interesting cases, please send figures or link to animations to <psuefd@gmail.com>.

## Acknowledgements ##
This work is supported by the Strategic Environmental Research and Development Program (SERDP, Award Number W74RDV70063408). 

## Disclaimer ##
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via <www.openfoam.com>, and owner of the OPENFOAM&copy;  and OpenCFD&copy; trade marks.

OPENFOAM&copy; is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via <www.openfoam.com>.
 
## Authors and contributors: ##
Xiaofeng Liu, Ph.D., P.E.  
Associate Professor  
Department of Civil and Environmental Engineering  
Institute of Computational and Data Sciences  
Penn State University

Yalan Song  
Graduate Research Assistant  
Department of Civil and Environmental Engineering  
Penn State University  

Yuncheng (Cloud) Xu (Former Ph.D. student at Penn State)  
Lecturer  
College of Water Resources and Civil Engineering  
China Agricultural University

## License ##
GPL v3 License



