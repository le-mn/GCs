# Tabulated Powerspectra

Most of these tabulate P(k) in [(Mpc/h)^3] against k in [h/Mpc]. The
normaliztion is arbitrary because most codes that read these will renormalize
to a given sigma8.

TODO: fill in cosmoogical paramters for these.

## `pk_EAGLE_norm.dat`

EAGLE CDM (Eagle project via Galform/Lan et al.)

Seems to have a 10% shift to higher k compared to the others? Different h?

## `pk_MillGas_norm.dat`

MillGas CDM (rom Galform)

## `pk_Mill_pch.dat`

Millennium CDM, from PCH 2008 code.

## `delta2_WDM33_COCO_WMAP7.dat`, `pk_WDM33_COCO_WMAP7.dat`

Renamed from `/cosma6/data/dp004/PROJECTS/COCO/COCO_WDM1/ICS/PS_TR33_0WMAP7.txt`

In the first file, power is Delta^2, not P(k). I created the second file,
`pk_WDM33_COCO_WMAP7.dat` with the same linear k, P(k) units as the other
powerspectra in this directory.

This has a somewhat strange powerlaw extension to high k, compared to `pk_WDMDove.dat`.

## `pk_WDMDove.dat`

EAGLE WDM (Eagle project via Galform/Lan et al., seemingly made by APC).

