# triNURBS
## Experimental software to test the use of volumetric textures in a trivariate NURBS context

The main idea is to construct a trivariate NURBS representation that could be used both as a representation of geometry and for analysis but with the key difference over existing IGA approaches that the volume is composed of a cellular or lattice geometry.

The main application of the discretisatoin is for additive manufacturing and for novel metamaterial devices where such multiscale geometries are common.  A real challenge is how to represent fine scale features efficiently and one approach is to use techniques from the CAGD community for texture mapping.  In texture mapping the aim is to create realistic surface effects (e.g. textiles, waves, materials) which are important for realistic rendered images and much effort has been focussed on how to produce realistic effects that can be rendered in reasonable timeframes.  This code extends the main ideas of texture mapping into a volumetric context (as opposed to surface based texture maps for the majority of rendering scenarios).

The main advantage of a texture mapping approach is that we do not explicity store a mesh/discretisation for the entire multiscale structure.  We only generate the fine scale geometry features when we need them, much like in rendering where fine scale features are only required if they are in the rendered view.  From the perspective of additive manufacturing (AM) this is highly advantageous where complex multiscale geometries are common and using the present approach we theoretically only need to generate fine scale geometry features when creating 'slicing' images for relevant layers in AM processes. 
