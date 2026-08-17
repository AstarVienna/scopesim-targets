# Defining positional information
Positional information is handled internally by `astropy.coordinates.SkyCoord`.
This page gives a quick overview on how to supply positional information.

## Possible ways to specify positions
In short:
* `astropy.coordinates.SkyCoord` instance (with or without a `distance`) — an **absolute** position.
* A 2-tuple of x and y offsets (in arcsec) — a **relative** position, see below.
* A dictionary of more specific information, optionally including a `distance`.

## Absolute vs. relative positions
These are distinct concepts and are kept as distinct types internally:
a relative offset is *never* silently promoted to an absolute coordinate.

* `position` is always absolute (world coordinates).
* `offset` is always relative to the **parent frame**: for components of a
  `TargetScene`, this is the scene's `position` (the *frame center*).
  Offsets follow the usual on-sky convention: +x = East, +y = North.

Note that the scene frame center is deliberately *not* the telescope pointing;
pointing, dither patterns and rotator angles are observation configuration and
live in `ScopeSim`, not in the target definition.

The `position` and `offset` arguments are mutually exclusive; specifying both
is a validation error.
The exception to this is the `Binary` class, where `position` implicitly applies
to the primary and `offset` to the secondary (i.e. `position` is the primary's
position, *not* the barycenter or photocenter).
Binaries with explicit positions for both components can be defined as a `StarField`.
If neither is specified, a relative offset of `(0, 0)` is assumed as a default.

## Position angles
Wherever a `position_angle` is accepted (e.g. `Binary` separations), it follows
the astronomical convention: measured **from North, through East**.

## Separations in physical units
Angular separations may be given in physical units (e.g. `0.1 AU`), in which
case a `distance` is required for the conversion.
A missing distance is an error; no default distance is assumed.

## Radial velocity and redshift
Line-of-sight motion is specified by exactly one of two mutually exclusive
attributes of the (absolute) position:

* `radial_velocity`: a velocity `Quantity`, for nearby objects.
  Used to Doppler-shift the spectrum.
* `z`: the (unitless) cosmological redshift, for distant objects.

If, for any reason, both are present, `radial_velocity` takes priority
(matching the ESO target specification).
The `distance` attribute is used **only** for geometric conversions
(parallax, physical separations); it is *never* used to infer a redshift,
so no cosmology is silently assumed.
Note that this feature is currently experimental and not all target subclasses
fully support it.
