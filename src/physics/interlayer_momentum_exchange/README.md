# interlayer momentum exchange

## Purpose and role in IMEX-SfloW2D

Pair-wide interfacial momentum drag without exchanging mass, species or thermal content.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`interlayer_drag_chezy_2d.f90`](./interlayer_drag_chezy_2d.f90) | Exact quadratic Chezy drag impulse |
| [`interlayer_momentum_exchange_2d.f90`](./interlayer_momentum_exchange_2d.f90) | Pair-wide config, method dispatch and atomic momenta |

## Physical equations

With $\boldsymbol w=\boldsymbol u_2-\boldsymbol u_1$, the Chézy drag is

$$
\boldsymbol D=f\rho_2|\boldsymbol w|\boldsymbol w,\quad
\dot{\boldsymbol m}_1=\boldsymbol D,\quad
\dot{\boldsymbol m}_2=-\boldsymbol D,\quad \boldsymbol m_\ell=M_\ell\boldsymbol u_\ell.
$$

For frozen layer masses and L2 density, $k=f\rho_2(1/M_1+1/M_2)$, hence the exact analytic relative-velocity flow is

$$
\boldsymbol w(t+\Delta t)=\frac{\boldsymbol w(t)}{1+k|\boldsymbol w(t)|\Delta t}.
$$

The equal and opposite impulses conserve the total horizontal momentum. The implemented pair update requires both layer depths to exceed the 1 mm contact threshold. The analytic method is supported; other methods explicitly reject use.

## Architectural and numerical contracts

Layer-pair configuration has a single owner. Paired momentum updates are atomic. This is interfacial, not basal, drag.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
