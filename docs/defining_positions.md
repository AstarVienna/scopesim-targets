# Defining positional information
Positional information is handles internally by `astropy.coordinates.SkyCoord`.
This pages gives a quick overview on how to supply positional information.

## Possible way to specify positions
In short:
* `astropy.coordinates.SkyCoord` instance (with or without a `distance`)
* A 2-tuple of x and y offsets (in arcsec) relative to the target center
* A dictionary of more specific information (see below), optionally including a `distance`

## Offsets
The `position` and `offset` arguments are mutually exclusive.
The exception to this is the `Binary` class, where `position` implicitly applies to the primary and `offset` to the secondary.
Binaries with explicit positions for both components can be defined as a `StarField`.
If neither is specified, a (relative) position of `(0, 0)` is assumed as a default.

## Radial velocity
If given as a `astropy.coordinates.SkyCoord` instance (directly in Python or via YAML constructor), it is possible to also include a `radial_velocity` as a keyword argument (see Astropy docs).
If used, it should be supplied as a Quantity and is used to shift the spectrum accordingly.
For cosmological distances, the `distance` attribute is used to apply cosmological redshift to the spectrum.
Note that this feature is currently experimental and not all target subclasses fully support it.
