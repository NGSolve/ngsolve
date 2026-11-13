Technical notes
===============

Short notes (about 4 pages) on methods and features of NGSolve. Each note
lives in ``docs/technotes/<number>_<name>/<name>.tex`` together with its figures.
The PDF and HTML versions are not stored in the repository, they are built
with the documentation (``make technotes``, needs ``latexmk`` for the PDF and
``pandoc`` for the HTML page).

**2026-01** `High-order de Rham sequences on mixed meshes <2026-01_derham_sequences.html>`_ (`PDF <2026-01_derham_sequences.pdf>`__)
    *Joachim Schöberl.*
    The spaces ``H1``, ``HCurl``, ``HDiv`` and ``L2`` on triangles, quads,
    tets, prisms and hexes: order convention, the ``type1`` and ``RT``
    flags, degrees of freedom, and the exact sequences they form.

**2026-02** `The TDNNS method for elasticity <2026-02_tdnns.html>`_ (`PDF <2026-02_tdnns.pdf>`__)
    *Astrid Pechstein, Michael Neunteufel, Joachim Schöberl.*
    Tangential-displacement normal-normal-stress elements: formulation,
    hybridization, and the space-creation syntax (``CreateStressSpace``,
    ``CreateMultiplierSpace``).

**2026-03** `Mass conserving mixed stress (MCS) methods <2026-03_mcs.html>`_ (`PDF <2026-03_mcs.pdf>`__)
    *Philip Lederer, Xaver Mooslechner, Jay Gopalakrishnan, Christopher
    Lackner, Joachim Schöberl.*
    Spaces and hybrid formulation of the MCS method for Stokes flow, and
    the space-creation syntax (``CreateStressSpace``,
    ``CreateMultiplierSpace``, ``CreateDivergenceSpace``).

**2026-04** `Nonlinear HHJ shell with frictional contact <2026-04_ghost_method.html>`_ (`PDF <2026-04_ghost_method.pdf>`__)
    *Michael Neunteufel, Joachim Schöberl.*
    A falling sheet draped over a moving sphere: nonlinear
    Hellan-Herrmann-Johnson shell, contact and friction, time stepping.

**2026-05** `Ordering for saddle point matrices in SparseCholesky <2026-05_sparsecholesky_ordering.html>`_ (`PDF <2026-05_sparsecholesky_ordering.pdf>`__)
    *Joachim Schöberl.*
    Zero diagonal entries and dense rows (``NumberSpace``) in the minimum
    degree ordering: the wait rule, the placement of constraints, the flags
    ``late``, ``densefac``, ``lastdofs`` and ``waitfrac`` of ``Inverse``,
    and the measured fill.

Adding a note
-------------

Take the next number ``<year>-<running number>``, create
``docs/technotes/<number>_<name>/<name>.tex``, put figures next to it, write
the number into its ``\date`` line, and add an entry above that links to
``<number>_<name>.html`` and ``<number>_<name>.pdf``.

The HTML page is made from the same source by ``tex2html.py`` and pandoc.
Check it in a browser, pandoc knows less LaTeX than LaTeX:

* figures included as PDF need a PNG or SVG of the same name next to them;
* lines starting with ``%html`` are read for the HTML page only, use them to
  redefine macros pandoc cannot handle (TikZ pictures, ``\\`` in a macro);
* put braces around macro arguments in formulas (``\tilde{\Pp}``).
