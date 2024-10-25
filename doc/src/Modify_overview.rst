Overview
========

The best way to add a new feature to LAMMPS is to find a similar
feature and look at the corresponding source and header files to
figure out what it does.  You will need some knowledge of C++ to
understand the high-level structure of LAMMPS and its class
organization.  Functions (class methods) that do actual computations
are mostly written in C-style code and operate on simple C-style data
structures (vectors and arrays).  A high-level overview of the
programming style choices in LAMMPS is :doc:`given elsewhere
<Developer_code_design>`.

Most of the new features described on the :doc:`Modify <Modify>` doc
page require you to write a new C++ derived class (excluding exceptions
described below, this can often be done by making small edits to existing
files).  Creating a new class requires 2 files, a source code file
(\*.cpp) and a header file (\*.h).  The derived class must provide
certain methods to work as a new option.  Depending on how different
your new feature is compared to existing features, you can either
derive from the base class itself, or from a derived class that
already exists.  Enabling LAMMPS to invoke the new class is as simple
as putting the two source files in the src directory and re-building
LAMMPS.

The advantage of C++ and its object-orientation is that all the code
and variables needed to define the new feature are in the 2 files you
write.  Thus, it should not make the rest of LAMMPS more complex or
cause bugs through unwanted side effects.

Here is a concrete example.  Suppose you write 2 files
``pair_foo.cpp`` and ``pair_foo.h`` that define a new class
``PairFoo`` which computes pairwise potentials described in the
classic 1997 :ref:`paper <Foo>` by Foo, *et al.* If you wish to invoke
those potentials in a LAMMPS input script with a command like:

.. code-block:: LAMMPS

   pair_style foo 0.1 3.5

then your ``pair_foo.h`` file should be structured as follows:

.. code-block:: c++

   #ifdef PAIR_CLASS
   // clang-format off
   PairStyle(foo,PairFoo);
   #else
   // clang-format on
   ...
   (class definition for PairFoo)
   ...
   #endif

where "foo" is the style keyword in the pair_style command, and
``PairFoo`` is the class name defined in your ``pair_foo.cpp`` and
``pair_foo.h`` files.

When you re-build LAMMPS, your new pairwise potential becomes part of
the executable and can be invoked with a pair_style command like the
example above.  Arguments like 0.1 and 3.5 can be defined and
processed by your new class.

As illustrated by this example, many features referred to in the
LAMMPS documentation are called a "style" of a particular command.

The :doc:`Modify page <Modify>` lists all the common styles in LAMMPS,
and discusses the header file for the base class that these styles
derive from.  Public variables in that file can be used and set by the
derived classes, and may also be used by the base class.  Sometimes
they are also accessed by the rest of LAMMPS.  Pure functions, which
means functions declared as virtual in the base class header file and
which are also set to 0, are functions you **must** implement in your
new derived class to give it the functionality LAMMPS expects. Virtual
functions that are not set to 0 are functions you may override or not.
Those are usually defined with an empty function body.

Additionally, new output options can be added directly to the
thermo.cpp, dump_custom.cpp, and variable.cpp files.  These are also
listed on the :doc:`Modify page <Modify>`.

Here are additional guidelines for modifying LAMMPS and adding new
functionality:

* Think about whether what you want to do would be better as a pre- or
  post-processing step.  Many computations are more easily and more
  quickly done that way.
* Do not try to do anything within the timestepping of a run that is not
  parallel.  For example, do not accumulate a bunch of data on a single
  processor and analyze it.  That would run the risk of seriously degrading
  the parallel efficiency.
* If your new feature reads arguments or writes output, make sure you
  follow the unit conventions discussed by the :doc:`units <units>`
  command.

----------

.. _Foo:

**(Foo)** Foo, Morefoo, and Maxfoo, J of Classic Potentials, 75, 345 (1997).
