# Defining positional information
Introduction TBA.

## Possible way to specify positions
In short:
* `astropy.coordinates.SkyCoord` instance
* A 2-tuple of x and y offsets (in arcsec) relative to the target center
* A dictionary of more specific information (see below)

## Offsets
The `position` and `offset` arguments are mutually exclusive.
The exception to this is the `Binary` class, where `position` implicitly applies to the primary and `offset` to the secondary.
Binaries with explicit positions for both components can be defined as a `StarField`.
If neither is specified, a (relative) position of `(0, 0)` is assumed as a default.

More TBA.
