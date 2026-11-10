Technical notes
===============

Short notes (about 4 pages) on methods and features of NGSolve. Each note
lives in ``docs/technotes/<number>_<name>/<name>.tex`` together with its figures.
The PDFs are not stored in the repository, they are built with the
documentation (``make technotes``, needs ``latexmk``).

**2026-01** `High-order de Rham sequences on mixed meshes <2026-01_derham_sequences.pdf>`_
    *Joachim Schöberl.*
    The spaces ``H1``, ``HCurl``, ``HDiv`` and ``L2`` on triangles, quads,
    tets, prisms and hexes: order convention, the ``type1`` and ``RT``
    flags, degrees of freedom, and the exact sequences they form.

**2026-02** `The TDNNS method for elasticity <2026-02_tdnns.pdf>`_
    *Astrid Pechstein, Michael Neunteufel, Joachim Schöberl.*
    Tangential-displacement normal-normal-stress elements: formulation,
    hybridization, and the space-creation syntax (``CreateStressSpace``,
    ``CreateMultiplierSpace``).

**2026-03** `Mass conserving mixed stress (MCS) methods <2026-03_mcs.pdf>`_
    *Philip Lederer, Xaver Mooslechner, Jay Gopalakrishnan, Christopher
    Lackner, Joachim Schöberl.*
    Spaces and hybrid formulation of the MCS method for Stokes flow, and
    the space-creation syntax (``CreateStressSpace``,
    ``CreateMultiplierSpace``, ``CreateDivergenceSpace``).

**2026-04** `Nonlinear HHJ shell with frictional contact <2026-04_ghost_method.pdf>`_
    *Michael Neunteufel, Joachim Schöberl.*
    A falling sheet draped over a moving sphere: nonlinear
    Hellan-Herrmann-Johnson shell, contact and friction, time stepping.

Adding a note
-------------

Take the next number ``<year>-<running number>``, create
``docs/technotes/<number>_<name>/<name>.tex``, put figures next to it, write
the number into its ``\date`` line, and add an entry above that links to
``<number>_<name>.pdf``.
