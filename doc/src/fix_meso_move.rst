.. index:: fix meso/move

fix meso/move command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID meso/move style args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* meso/move = style name of this fix command
* style = *linear* or *wiggle* or *rotate* or *variable*

  .. parsed-literal::

       *linear* args = Vx Vy Vz
         Vx,Vy,Vz = components of velocity vector (velocity units), any component can be specified as NULL
       *wiggle* args = Ax Ay Az period
         Ax,Ay,Az = components of amplitude vector (distance units), any component can be specified as NULL
         period = period of oscillation (time units)
       *rotate* args = Px Py Pz Rx Ry Rz period
         Px,Py,Pz = origin point of axis of rotation (distance units)
         Rx,Ry,Rz = axis of rotation vector
         period = period of rotation (time units)
       *variable* args = v_dx v_dy v_dz v_vx v_vy v_vz
         v_dx,v_dy,v_dz = 3 variable names that calculate x,y,z displacement as function of time, any component can be specified as NULL
         v_vx,v_vy,v_vz = 3 variable names that calculate x,y,z velocity as function of time, any component can be specified as NULL

* zero or more keyword/value pairs may be appended
* keyword = *units*

  .. parsed-literal::

       *units* value = *box* or *lattice*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 boundary meso/move wiggle 3.0 0.0 0.0 1.0 units box
   fix 2 boundary meso/move rotate 0.0 0.0 0.0 0.0 0.0 1.0 5.0
   fix 2 boundary meso/move variable v_myx v_myy NULL v_VX v_VY NULL

Description
"""""""""""

Perform updates of position, velocity, internal energy and local
density for mesoscopic particles in the group each timestep using the
specified settings or formulas, without regard to forces on the
particles. This can be useful for boundary, solid bodies or other
particles, whose movement can influence nearby particles.

The operation of this fix is exactly like that described by the
:doc:`fix move <fix_move>` command, except that particles' density,
internal energy and extrapolated velocity are also updated.

.. note::

   The particles affected by this fix should not be time integrated by
   other fixes (e.g. :doc:`fix sph <fix_sph>`, :doc:`fix
   sph/stationary <fix_sph_stationary>`), since that will change their
   positions and velocities twice.

.. note::

   As particles move due to this fix, they will pass through periodic
   boundaries and be remapped to the other side of the simulation box,
   just as they would during normal time integration (e.g. via the
   :doc:`fix sph <fix_sph>` command).  It is up to you to decide
   whether periodic boundaries are appropriate with the kind of
   particle motion you are prescribing with this fix.

.. note::

   As discussed below, particles are moved relative to their initial
   position at the time the fix is specified.  These initial coordinates
   are stored by the fix in "unwrapped" form, by using the image flags
   associated with each particle.  See the :doc:`dump custom <dump>` command
   for a discussion of "unwrapped" coordinates.  See the Atoms section of
   the :doc:`read_data <read_data>` command for a discussion of image flags
   and how they are set for each particle.  You can reset the image flags
   (e.g. to 0) before invoking this fix by using the :doc:`set image <set>`
   command.

----------

The *linear* style moves particles at a constant velocity, so that their
position *X* = (x,y,z) as a function of time is given in vector
notation as

.. parsed-literal::

   X(t) = X0 + V \* delta

where *X0* = (x0,y0,z0) is their position at the time the fix is
specified, *V* is the specified velocity vector with components
(Vx,Vy,Vz), and *delta* is the time elapsed since the fix was
specified.  This style also sets the velocity of each particle to V =
(Vx,Vy,Vz).  If any of the velocity components is specified as NULL,
then the position and velocity of that component is time integrated
the same as the :doc:`fix sph <fix_sph>` command would perform, using
the corresponding force component on the particle.

Note that the *linear* style is identical to using the *variable*
style with an :doc:`equal-style variable <variable>` that uses the
vdisplace() function.  E.g.

.. code-block:: LAMMPS

   variable V equal 10.0
   variable x equal vdisplace(0.0,$V)
   fix 1 boundary move variable v_x NULL NULL v_V NULL NULL

The *wiggle* style moves particles in an oscillatory fashion, so that
their position *X* = (x,y,z) as a function of time is given in vector
notation as

.. parsed-literal::

   X(t) = X0 + A sin(omega\*delta)

where *X0* = (x0,y0,z0) is their position at the time the fix is
specified, *A* is the specified amplitude vector with components
(Ax,Ay,Az), *omega* is 2 PI / *period*, and *delta* is the time
elapsed since the fix was specified.  This style also sets the
velocity of each particle to the time derivative of this expression.
If any of the amplitude components is specified as NULL, then the
position and velocity of that component is time integrated the same as
the :doc:`fix sph <fix_sph>` command would perform, using the
corresponding force component on the particle.

Note that the *wiggle* style is identical to using the *variable*
style with :doc:`equal-style variables <variable>` that use the
swiggle() and cwiggle() functions.  E.g.

.. code-block:: LAMMPS

   variable A equal 10.0
   variable T equal 5.0
   variable omega equal 2.0*PI/$T
   variable x equal swiggle(0.0,$A,$T)
   variable v equal v_omega*($A-cwiggle(0.0,$A,$T))
   fix 1 boundary move variable v_x NULL NULL v_v NULL NULL

The *rotate* style rotates particles around a rotation axis *R* =
(Rx,Ry,Rz) that goes through a point *P* = (Px,Py,Pz).  The *period* of
the rotation is also specified.  The direction of rotation for the
particles around the rotation axis is consistent with the right-hand
rule: if your right-hand thumb points along *R*, then your fingers wrap
around the axis in the direction of rotation.

This style also sets the velocity of each particle to (omega cross
Rperp) where omega is its angular velocity around the rotation axis and
Rperp is a perpendicular vector from the rotation axis to the particle.

The *variable* style allows the position and velocity components of
each particle to be set by formulas specified via the
:doc:`variable <variable>` command.  Each of the 6 variables is
specified as an argument to the fix as v_name, where name is the
variable name that is defined elsewhere in the input script.

Each variable must be of either the *equal* or *atom* style.
*Equal*\ -style variables compute a single numeric quantity, that can be
a function of the timestep as well as of other simulation values.
*Atom*\ -style variables compute a numeric quantity for each particle, that
can be a function per-atom quantities, such as the particle's position, as
well as of the timestep and other simulation values.  Note that this
fix stores the original coordinates of each particle (see note below) so
that per-atom quantity can be used in an atom-style variable formula.
See the :doc:`variable <variable>` command for details.

The first 3 variables (v_dx,v_dy,v_dz) specified for the *variable*
style are used to calculate a displacement from the particle's original
position at the time the fix was specified.  The second 3 variables
(v_vx,v_vy,v_vz) specified are used to compute a velocity for each
particle.

Any of the 6 variables can be specified as NULL.  If both the
displacement and velocity variables for a particular x,y,z component
are specified as NULL, then the position and velocity of that
component is time integrated the same as the :doc:`fix sph <fix_sph>`
command would perform, using the corresponding force component on the
particle.  If only the velocity variable for a component is specified
as NULL, then the displacement variable will be used to set the
position of the particle, and its velocity component will not be
changed. If only the displacement variable for a component is
specified as NULL, then the velocity variable will be used to set the
velocity of the particle, and the position of the particle will be
time integrated using that velocity.

The *units* keyword determines the meaning of the distance units used
to define the *linear* velocity and *wiggle* amplitude and *rotate*
origin.  This setting is ignored for the *variable* style.  A *box*
value selects standard units as defined by the :doc:`units <units>`
command, e.g. velocity in Angstroms/fs and amplitude and position
in Angstroms for units = real.  A *lattice* value means the velocity
units are in lattice spacings per time and the amplitude and position
are in lattice spacings.  The :doc:`lattice <lattice>` command must have
been previously used to define the lattice spacing.  Each of these 3
quantities may be dependent on the x,y,z dimension, since the lattice
spacings can be different in x,y,z.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the original coordinates of moving particles to :doc:`binary restart files <restart>`, as well as the initial timestep, so that
the motion can be continuous in a restarted simulation.  See the
:doc:`read_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file, so that the
operation of the fix continues in an uninterrupted fashion.

.. note::

   Because the move positions are a function of the current
   timestep and the initial timestep, you cannot reset the timestep to a
   different value after reading a restart file, if you expect a fix move
   command to work in an uninterrupted fashion.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix produces a per-atom array which can be accessed by various
:doc:`output commands <Howto_output>`.  The number of columns for each
atom is 3, and the columns store the original unwrapped x,y,z coords
of each particle.  The per-atom values can be accessed on any timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the DPD-SMOOTH package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This fix requires that atoms store density and internal energy as
defined by the :doc:`atom_style sph <atom_style>` command.

All particles in the group must be mesoscopic SPH/SDPD particles.

.. versionchanged:: TBD

This fix is incompatible with deformation controls that remap velocity,
for instance the *remap v* option of :doc:`fix deform <fix_deform>`.

Related commands
""""""""""""""""

:doc:`fix move <fix_move>`, :doc:`fix sph <fix_sph>`,
:doc:`displace_atoms <displace_atoms>`

Default
"""""""

The option default is units = lattice.
