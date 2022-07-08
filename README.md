TO DO: ADD EXPLANATION

TO DO: HIDE UNNECESSARY FOLDERS

TO DO: ADD GIFS

# Script to analyze force harmonics in a segmented PM machine

This script is used for my Master Thesis at the TU Delft, called *Forces and Vibrations in a Modular HVDC machine*.
It aims to set up an analytical model of the magnetic flux density distribution in the air gap, and apply Maxwell's stress tensor to obtain the resulting forces.

The method is based on *Noise in Polyphase Electric Motors* By J.F. Gieras, 2D magnetic flux models by Z.Q. Zhu, and work by Jean Le Besnerais and Mostafa Valavi. 
The complex relative permeance function is based on the papers by Damir Žarko, and other references can be found in the thesis report.

Casper Klop, July 2022

---
The specific machine type assumed in this script is a Surface-mounted PM machine with Fractional-Slot Concentrated Windings, with the stator segmented in multiple 3-phase sections.
The model is tested with the parameters set in [Set_params.m](Set_params.m), but effort is put in parametrizing all equations so other models can be analyzed as well. Implicit, underlying assumptions in some of the equations might result in errors if another machine is analyzed.

It is assumed that all segments are equal, and all calculations are done for one segment only.

## Flux density calculations
The magnetic flux density is calculated close to the stator inner radius, as the thesis is focused on calculating the forces acting on the stator segments.
First, the magnetic flux density resulting from the Permanent Magnets is calculated, assuming a slotless stator, which only has a radial component:

![Radial flux density due to rotor](GIFs/FluxPMRadial.gif)

Second, the magnetic flux density resulting from the armature reaction is calculated, assuming a slotless stator, which has both a radial and a tangential component:

![Radial flux density due to stator](GIFs/FluxArmRadial.gif)

Neglecting saturation, the total magnetic flux density for a slotless machine is calculated by the linear superposition of the two magnetic fields.

## Influence of slots and segment gaps


