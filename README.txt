spatialnde
----------

The purpose of this package is to facilitate data projection from
imaging NDE modalities (such as thermography) onto surface
parameterizations of CAD models.

It leverages the OpenCV package to do much of the real work
such as camera distortion correction and pose estimation

It is intended as a prototype. The successor package spatialnde2
should generally be used in most cases.

Primary components/directory structure
--------------------------------------

spatialnde/         The Python package
pt_steps            Spatial processing steps for performinag
                    flat perspective corrections and
		    distortion correction
scripts/            Scripts to facilitate OpenCV camera calibration,
                    coming up with fiducial (landmark) names,
doc/                Documentation of some of the projection calculations
dataguzzler/        Dataguzzler configuration for real time
                    perspective correction, real-time UV projection
		    of camera data, definition of landmarks, and
		    recording of UV image "stacks" (data cubes"
demos/              Various demonstration scripts.

SpatialNDE Python package structure
-----------------------------------
spatialnde                    Main package; performs camera spatial
                              calibration and distortio correction
			      using OpenCV; defines landmarks, objects,
			      rays, and parameterizations; accelerates
			      geometry calculations including 3D ray
			      tracing
spatialnde.cadpart            Ability to represent CAD objects by 3D
                              surface and parameterization meshes.
spatialnde.cadpart.loaders    Ability to read .x3d and .stl geometry files
spatialnde.dataguzzler        Creation of 3D dataguzzler channels;
                              extraction of landmarks from dataguzzler channels
spatialnde.opencascade        Loading of curvature data via OpenCascade
                              geometry kernel
spatialnde.exporters          Ability to write VRML and X3D geometry files
