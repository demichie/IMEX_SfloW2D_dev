"""Explicit production policies over UNMODIFIED checksum-pinned spatial cores.

The historical SSPRK2 prototype projects cells at dry_h_tol and uses their
original desingularized velocity in the normal HP momentum bounds. Production
instead retains positive cell mass down to machine epsilon. The explicitly
authorized correction uses zero auxiliary reconstruction velocity at
h<=dry_h_tol, matching the compatible normal HP momentum interval. A separate
authorized guard omits only local curvature below that existing threshold.
Original volumetric reconstruction candidates remain unchanged. A separately
authorized final accepted-state projection removes only unresolved conservative
momenta; check_cases performs it after recording/checking the raw assembly.
Q1 bed, eta, blend weights, quadrature, path terms and flux functions remain
the clean kernels. This is NOT a new clean reference, nor a claim of identical
historical semantics; historical SSPRK2 is compared separately.

On compact domains production also copies ALL LS derivatives across the
two-cell boundary strip, whereas the prototype copies pure derivatives only
along their differentiation direction. Explicitly match this boundary policy
before forming G; never change the authoritative Q1 face elevations.
"""
from dataclasses import replace

import numpy as np


def production_curvature_rhs(q, bed, dx, dy, time, params, cache, core):
    """Assemble the unchanged spatial RHS plus curvature only at resolved depth.

    Evaluate the quadratic velocity contraction AFTER selecting h>dry_h_tol.
    Merely multiplying a full-array contraction by zero would still overflow
    for unresolved dry velocities. This is the explicitly user-authorized
    source guard, not a modification of either checksum-pinned clean core.
    The historical SSPRK2 oracle continues to use the original unguarded core.
    """
    spatial_params = replace(params, curvature_term=False)
    rhs = core.rhs2d(q, bed, dx, dy, time, spatial_params, stage_cache=cache)
    if not params.curvature_term:
        return rhs
    geometry = cache["slope_geometry"]
    physical = core.primitive2d(cache["q"], params)
    wet = physical["h"] > params.dry_h_tol
    u, v = physical["u"][wet], physical["v"][wet]
    acceleration = (geometry.Bxx_center[wet]*u**2
                    + 2*geometry.Bxy_center[wet]*u*v
                    + geometry.Byy_center[wet]*v**2)
    gravity = geometry.G_center[wet] if params.slope_correction else 1.0
    coefficient = -gravity*physical["M"][wet]*acceleration
    rhs[..., 1][wet] += coefficient*geometry.Bx_center[wet]
    rhs[..., 2][wet] += coefficient*geometry.By_center[wet]
    return rhs


def production_slope_geometry(bed, dx, dy, core):
    """Use the clean LS kernels with production's tensor boundary-strip copying."""
    geometry = core.build_slope_geometry_2d(bed, dx, dy)
    derivatives = [getattr(geometry, name).copy() for name in
        ("Bx_center", "By_center", "Bxx_center", "Bxy_center", "Byy_center")]
    ny, nx = bed.shape
    if nx < 5 or ny < 5:
        raise AssertionError("2D comparator requires resolved axes of at least five cells")
    for field in derivatives:
        field[:, :2] = field[:, 2:3]
        field[:, -2:] = field[:, -3:-2]
        field[:2] = field[2:3]
        field[-2:] = field[-3:-2]
    bx, by, bxx, bxy, byy = derivatives
    G = 1/(1+bx*bx+by*by)
    Gx = np.column_stack((G[:, 0], 0.5*(G[:, :-1]+G[:, 1:]), G[:, -1]))
    Gy = np.vstack((G[0], 0.5*(G[:-1]+G[1:]), G[-1]))
    return core.SlopeGeometry2D(G, Gx, Gy, bx, by, bxx, bxy, byy)


def production_dry_momentum_policy(cache, params, core):
    """Reassemble normal momentum bounds with zero velocity in dry cell averages.

    Direct auxiliary velocity candidates and the final admissibility interval
    both use zero at the EXISTING dry threshold, as explicitly authorized by
    the user. Volumetric candidates and conservative cell states are retained.
    """
    count = 0
    for axis in ("x", "y"):
        for entry in cache[axis]:
            line = entry["cache1d"]
            data = line["data"]
            h, u = data["prim_q"][3], data["prim_q"][1]
            dry = h <= params.dry_h_tol
            if not np.any(dry & (u != 0)):
                continue
            count += int(np.count_nonzero(dry & (u != 0)))
            _, _, _, _, _, _, hu_minus, hu_plus, um, up = core.reconstruct_state(
                line["q"], params, return_physical=True, primitive_values=data["prim_q"])
            # Authorized correction: direct auxiliary velocity candidates use
            # the same dry centre convention as final normal bounds. Retain
            # the original h*u volumetric candidates and conservative cells.
            safe = np.where(dry, 0, u)
            um, up = core.reconstruct_array(safe, params)
            hm, hp = core.reconstruct_array(h, params)
            hm = np.maximum(hm, 0); hp = np.maximum(hp, 0)
            bed = entry["bed"]
            heta_minus = np.maximum(data["etaL_cell"]-bed.B_face[:-1], 0)
            heta_plus = np.maximum(data["etaR_cell"]-bed.B_face[1:], 0)
            weight = data["w_eta"]
            hm = (1-weight)*hm + weight*heta_minus
            hp = (1-weight)*hp + weight*heta_plus
            weight_momentum = weight*(np.abs(np.diff(bed.B_face)) > params.bed_step_tol)
            hu_minus = (1-weight_momentum)*hu_minus + weight_momentum*heta_minus*um
            hu_plus = (1-weight_momentum)*hu_plus + weight_momentum*heta_plus*up
            neighbors = np.stack((safe, np.r_[safe[0], safe[:-1]], np.r_[safe[1:], safe[-1]]))
            lo, hi = neighbors.min(axis=0), neighbors.max(axis=0)
            target = 0.5*(hu_plus-hu_minus)-safe*0.5*(hp-hm)
            lower = np.maximum((safe-hi)*hm, (lo-safe)*hp)
            upper = np.minimum((safe-lo)*hm, (hi-safe)*hp)
            delta = np.minimum(np.maximum(target, lower), upper)
            # Explicit intermediate arrays avoid FMA in the zero-bound cancellation.
            mean_minus, mean_plus = safe*hm, safe*hp
            hu_minus, hu_plus = mean_minus-delta, mean_plus+delta
            am, ap = core.reconstruct_array(data["prim_q"][5], params)
            qm = core.state_from_h_alpha_hu(hm, am, hu_minus, params)
            qp = core.state_from_h_alpha_hu(hp, ap, hu_plus, params)
            vm, vp = core.primitive(qm, params)[1], core.primitive(qp, params)[1]
            data.update(qL_cell=qm, qR_cell=qp, uL_cell=vm, uR_cell=vp,
                        huL_cell=hu_minus, huR_cell=hu_plus)
            # Reassemble interior faces first. The final closure pass below
            # copies reconstructed interior states into the external ghosts.
            for key, first, interior, last in (
                ("qL_face", data["qL_face"][:1], qp, None),
                ("qR_face", None, qm, data["qR_face"][-1:]),
                ("uL_face", data["uL_face"][:1], vp, None),
                ("uR_face", None, vm, data["uR_face"][-1:]),
                ("huL_face", data["huL_face"][:1], hu_plus, None),
                ("huR_face", None, hu_minus, data["huR_face"][-1:])):
                data[key] = np.concatenate(tuple(v for v in (first, interior, last) if v is not None), axis=0)
            data["primL_face"] = core.primitive(data["qL_face"], params)
            data["primR_face"] = core.primitive(data["qR_face"], params)
            data["BhydL_face"] = data["etaL_face"]-data["primL_face"][3]
            data["BhydR_face"] = data["etaR_face"]-data["primR_face"][3]
            H, fL, fR, minus, plus = core.cu_flux_from_states(
                data["qL_face"], data["qR_face"], params,
                uL_aux=data["uL_face"], uR_aux=data["uR_face"],
                primitive_L=data["primL_face"], primitive_R=data["primR_face"],
                grav_coeff=line["G_face"])
            line.update(H0=H, fL=fL, fR=fR, a_minus=minus, a_plus=plus)
    return count


def auxiliary_tangential_faces(qslice, normal, params, core):
    """Reconstruct only auxiliary tangential velocity from dry-safe centres."""
    velocity = core._tangential_velocity_slice(qslice, normal, params)
    h = core.raw_thickness2d(qslice, params)
    safe = np.where(h <= params.dry_h_tol, 0, velocity)
    return core._reconstruct_scalar_faces(safe, params)


def correct_auxiliary_tangential_transport(rhs, cache, dx, dy, params, core):
    """Replace only tangential CU fluxes affected by dry auxiliary neighbours.

    The pinned core reconstructs tangential raw velocities. Express the
    authorized production correction as a conservative flux DIFFERENCE;
    retain all clean normal flux, path and source kernels unchanged.
    """
    q = cache['q']
    for axis, spacing, component in (('x', dx, 2), ('y', dy, 1)):
        for index, entry in enumerate(cache[axis]):
            row = q[index] if axis == 'x' else q[:, index]
            velocity = core._tangential_velocity_slice(row, axis, params)
            h = core.raw_thickness2d(row, params)
            if not np.any((h <= params.dry_h_tol) & (velocity != 0)):
                continue
            line = entry['cache1d']; data = line['data']
            diagnostic = {**data, 'a_minus': line['a_minus'], 'a_plus': line['a_plus']}
            original = core._tangential_flux_from_diag(row, axis, params, diagnostic)
            left, right = auxiliary_tangential_faces(row, axis, params, core)
            ml, mr = data['qL_face'][:, 0]*left, data['qR_face'][:, 0]*right
            am, ap = line['a_minus'], line['a_plus']
            denominator = ap-am
            corrected = (ap*data['uL_face']*ml-am*data['uR_face']*mr+
                         ap*am*(mr-ml))/np.where(denominator > params.speed_eps, denominator, 1)
            degenerate = denominator <= params.speed_eps
            corrected[degenerate] = 0.5*(data['uL_face'][degenerate]*ml[degenerate]+
                                         data['uR_face'][degenerate]*mr[degenerate])
            difference = corrected-original
            change = -(difference[1:]-difference[:-1])/spacing
            if axis == 'x':
                rhs[index, :, component] += change
            else:
                rhs[:, index, component] += change


def final_dry_face_momenta(cache, params, core):
    """Close dry diagnostics and copy final interior traces into external ghosts.

    Production copies the FINAL reconstructed state, not the cell average,
    at transmissive domain faces (reconstruction_2d's external ghost pass).
    This matters when a nonzero near-dry average has a dry reconstructed face.
    Preserve this explicitly instead of assuming every boundary stays quiet.
    """
    for axis in ("x", "y"):
        for entry in cache[axis]:
            data = entry["cache1d"]["data"]
            for side in ("L", "R"):
                for location in ("cell", "face"):
                    dry = data[f"q{side}_{location}"][..., 0] == 0
                    data[f"hu{side}_{location}"][dry] = 0
            for prefix in ("q", "u", "hu", "eta"):
                data[f"{prefix}L_face"][0] = data[f"{prefix}R_face"][0]
                data[f"{prefix}R_face"][-1] = data[f"{prefix}L_face"][-1]
            for side in ("L", "R"):
                data[f"prim{side}_face"] = core.primitive(data[f"q{side}_face"], params)
                data[f"Bhyd{side}_face"] = data[f"eta{side}_face"]-data[f"prim{side}_face"][3]
            line = entry["cache1d"]
            H, fL, fR, minus, plus = core.cu_flux_from_states(
                data["qL_face"], data["qR_face"], params,
                uL_aux=data["uL_face"], uR_aux=data["uR_face"],
                primitive_L=data["primL_face"], primitive_R=data["primR_face"],
                grav_coeff=line["G_face"])
            line.update(H0=H, fL=fL, fR=fR, a_minus=minus, a_plus=plus)


def stationary_roundoff_policy(rhs, cache, dx, dy, params, core):
    """Independently apply the existing 64-epsilon exact-rest pressure allowance.

    No moving cell/stencil qualifies, even at arbitrarily small nonzero speed.
    Scalars are untouched. This is the a8148e4 production policy, absent from
    the historical clean prototype, not a general momentum filter.
    """
    q = cache["q"]
    resting = np.all(q[..., 1:3] == 0, axis=-1)
    scales = []
    for axis, normal in (("x", 0), ("y", 1)):
        line_masks, line_scales = [], []
        for entry in cache[axis]:
            line = entry["cache1d"]; data = line["data"]
            left, right = auxiliary_tangential_faces(
                q[len(line_masks)] if axis == 'x' else q[:, len(line_masks)], axis, params, core)
            mt_left = data["qL_face"][:, 0]*left
            mt_right = data["qR_face"][:, 0]*right
            face_rest = (data["qL_face"][:, 1] == 0) & (data["qR_face"][:, 1] == 0)
            face_rest &= (mt_left == 0) & (mt_right == 0)
            line_masks.append(face_rest[:-1] & face_rest[1:])
            face_scale = np.zeros(len(face_rest))
            for side in ("L", "R"):
                primitive = data[f"prim{side}_face"]
                h, gamma = primitive[3], primitive[8]
                value = abs(line["G_face"]*gamma)*h*np.maximum(h, abs(data[f"eta{side}_face"]))
                face_scale = np.maximum(face_scale, value)
            line_scales.append(np.maximum(face_scale[:-1], face_scale[1:]))
        resting &= np.stack(line_masks, axis=normal)
        scales.append(np.stack(line_scales, axis=normal))
    for component, scale, spacing in ((1, scales[0], dx), (2, scales[1], dy)):
        quiet = resting & (np.abs(rhs[..., component]) <= 64*np.finfo(float).eps*scale/spacing)
        rhs[..., component][quiet] = 0
