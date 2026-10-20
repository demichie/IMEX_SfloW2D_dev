"""Independent NumPy implementation of the current Fortran KT spatial operator.

Q=(mixture mass, x momentum, y momentum, solid mass). The supported closure is
an isothermal liquid with one solid species (including pure liquid). This core
uses no Fortran output as input to its fluxes. Historical HP-PCCU routines are
used only for the unchanged velocity regularization and bed geometry helpers.
External faces use the production copied-ghost/outflow-only convention. Source
boundaries, gas closures and rheology are outside this reference's scope.
"""
from dataclasses import dataclass
import numpy as np
import hp_pccu_2d_core as geometry

EPS = np.finfo(float).eps


@dataclass
class Parameters:
    rho_c: float = 1000.0
    rho_s: float = 2500.0
    rho_a: float = 0.0
    g: float = 9.81
    eps_sing: float = 1e-10
    dry_h_tol: float = 1e-10
    speed_eps: float = 1e-14
    theta: float = 1.3
    limiter: int = 3
    reconstr_coeff: float = 1.0
    reconstruction: str = 'froude_blend'
    slope_correction: bool = False
    curvature_term: bool = False
    cfl: float = 0.24


class VertexBed(geometry.BedGeometry2D):
    """Q1 geometry with the production four-vertex cell-mean operation order."""
    def __init__(self, vertices):
        v = np.asarray(vertices,float)
        super().__init__(0.5*(v[:-1]+v[1:]),0.5*(v[:,:-1]+v[:,1:]))
        object.__setattr__(self,'vertices',v)

    @property
    def B_center(self):
        v = self.vertices
        return 0.25*(v[:-1,:-1]+v[:-1,1:]+v[1:,:-1]+v[1:,1:])


def primitive(q, p):
    mass = q[..., 0]
    xs = np.divide(q[..., 3], mass, out=np.zeros_like(mass), where=mass > 0)
    if np.any(xs < -4096*EPS) or np.any(xs > 1+4096*EPS):
        raise ValueError('inadmissible component mass fraction')
    xs = np.clip(xs, 0, 1)
    density = 1 / ((1-xs)/p.rho_c + xs/p.rho_s)
    h = np.maximum(mass/density, 0)
    u = geometry.desingularized_velocity(mass, q[..., 1], p.eps_sing)
    v = geometry.desingularized_velocity(mass, q[..., 2], p.eps_sing)
    return h, u, v, xs


def minmod(a, b):
    return np.where(a*b > 0, np.sign(a)*np.minimum(abs(a), abs(b)), 0)


def slope(values, axis, p):
    """Limited increment across a cell on a uniform grid, before coefficient."""
    v = np.moveaxis(np.asarray(values), axis, -1)
    out = np.zeros_like(v)
    if v.shape[-1] < 3 or p.limiter == 0:
        return np.moveaxis(out, -1, axis)
    forward, backward = v[..., 2:]-v[..., 1:-1], v[..., 1:-1]-v[..., :-2]
    centered = 0.5*(v[..., 2:]-v[..., :-2])
    s = limited(forward, backward, centered, p)
    out[..., 1:-1] = s
    return np.moveaxis(out, -1, axis)


def limited(forward, backward, centered, p):
    if p.limiter == 0:
        return np.zeros_like(centered)
    if p.limiter == 1:
        return minmod(forward, backward)
    if p.limiter == 2:
        a, b = minmod(forward, 2*backward), minmod(2*forward, backward)
        return np.where(abs(a) >= abs(b), a, b)
    if p.limiter in (3, 4):
        return minmod(centered, (p.theta if p.limiter == 3 else 2)*minmod(forward, backward))
    if p.limiter == 5:
        return centered
    if p.limiter == 6:
        return forward
    if p.limiter == 7:
        return backward
    raise ValueError('unsupported limiter')


def velocity_traces(q, h, value, component, axis, p):
    """Physical boundary auxiliary velocities use a half-cell ghost stencil."""
    safe = np.where(h <= p.dry_h_tol, 0, value)
    increment = slope(safe, axis, p)
    values = np.moveaxis(safe, axis, -1)
    increments = np.moveaxis(increment, axis, -1)
    raw = np.divide(q[..., component+1], q[..., 0], out=np.zeros_like(h), where=q[..., 0] > 0)
    raw = np.moveaxis(raw, axis, -1)
    exterior_minus, exterior_plus = raw[..., 0].copy(), raw[..., -1].copy()
    if component == (0 if axis == 1 else 1):
        exterior_minus = np.minimum(exterior_minus, 0)
        exterior_plus = np.maximum(exterior_plus, 0)
    if values.shape[-1] > 1:
        increments[..., 0] = limited(values[..., 1]-values[..., 0],
            2*(values[..., 0]-exterior_minus), (values[..., 1]-exterior_minus)/1.5, p)
        increments[..., -1] = limited(2*(exterior_plus-values[..., -1]),
            values[..., -1]-values[..., -2], (exterior_plus-values[..., -2])/1.5, p)
    minus, plus = values-0.5*p.reconstr_coeff*increments, values+0.5*p.reconstr_coeff*increments
    minus[..., 0], plus[..., -1] = exterior_minus, exterior_plus
    return np.moveaxis(minus, -1, axis), np.moveaxis(plus, -1, axis)


def face_pair(minus, plus, axis):
    """Map cell endpoints to faces, copying each final interior boundary trace."""
    m, r = np.moveaxis(minus, axis, 0), np.moveaxis(plus, axis, 0)
    left = np.concatenate((m[:1], r), axis=0)
    right = np.concatenate((m, r[-1:]), axis=0)
    return np.moveaxis(left, 0, axis), np.moveaxis(right, 0, axis)


def positive_eta(eta, increment, bed_minus, bed_plus):
    m, r = eta-0.5*increment, eta+0.5*increment
    high, low = r < bed_plus, m < bed_minus
    r = np.where(high, bed_plus, r)
    m = np.where(high, 2*eta-r, m)
    m = np.where(~high & low, bed_minus, m)
    r = np.where(~high & low, 2*eta-m, r)
    return m, r


def hydrostatic_path(hm, hp, gm, gp, em, ep, Gm, Gp, p):
    """Three-point Gauss integral of -G*(Gamma*h*deta + h²/2*dGamma)."""
    gm, gp = np.array(gm, copy=True), np.array(gp, copy=True)
    left_dry, right_dry = hm <= p.dry_h_tol, hp <= p.dry_h_tol
    gm = np.where(left_dry & ~right_dry, gp, gm)
    gp = np.where(right_dry & ~left_dry, gm, gp)
    result = np.zeros_like(hm)
    abscissae = (0.5-0.5*np.sqrt(3/5), 0.5, 0.5+0.5*np.sqrt(3/5))
    for s, w in zip(abscissae, (5/18, 4/9, 5/18)):
        h, gamma = hm+s*(hp-hm), gm+s*(gp-gm)
        G = Gm+s*(Gp-Gm)
        result += w*(-G*(gamma*h*(ep-em)+0.5*h*h*(gp-gm)))
    return result


def direction(q, bed, spacing, axis, p, Gface):
    h, u, v, xs = primitive(q, p)
    bed_center = bed.B_center
    bface = bed.Bx_face if axis == 1 else bed.By_face
    bm = np.take(bface, range(h.shape[axis]), axis=axis)
    bp = np.take(bface, range(1, h.shape[axis]+1), axis=axis)
    eta = h+bed_center
    de = p.reconstr_coeff*slope(eta, axis, p)
    dh = p.reconstr_coeff*slope(h, axis, p)
    em, ep = positive_eta(eta, de, bm, bp)
    if p.reconstruction == 'h_x':
        hm, hp = np.maximum(h-0.5*dh, 0), np.maximum(h+0.5*dh, 0)
    else:
        increment = de
        if p.reconstruction == 'minmod':
            increment = minmod(de, dh+bp-bm)
        elif p.reconstruction == 'froude_blend':
            un = u if axis == 1 else v
            fr = np.divide(abs(un), np.sqrt(p.g*h), out=np.zeros_like(h), where=h > p.dry_h_tol)
            z = np.clip((fr-0.8)/0.4, 0, 1)
            weight = z*z*(3-2*z)
            increment = (1-weight)*de+weight*(dh+bp-bm)
        elif p.reconstruction != 'eta_x':
            raise ValueError('unknown trace reconstruction')
        tm, tp = positive_eta(eta, increment, bm, bp)
        hm, hp = np.maximum(tm-bm, 0), np.maximum(tp-bp, 0)
    hL, hR = face_pair(hm, hp, axis)
    eL, eR = face_pair(em, ep, axis)
    states, velocities, gammas, densities = [], [], [], []
    for side in (0, 1):
        depth = hL if side == 0 else hR
        traces = []
        for component, value in enumerate((u, v, xs)):
            if component < 2:
                minus, plus = velocity_traces(q,h,value,component,axis,p)
            else:
                delta = 0.5*p.reconstr_coeff*slope(value,axis,p)
                minus, plus = value-delta, value+delta
            traces.append(face_pair(minus,plus,axis)[side])
        uf, vf, fraction = traces
        fraction = np.clip(fraction, 0, 1)
        dry = depth <= p.dry_h_tol
        uf, vf = np.where(dry, 0, uf), np.where(dry, 0, vf)
        depth = np.where(dry, 0, depth)
        rho = 1/((1-fraction)/p.rho_c+fraction/p.rho_s)
        mass = rho*depth
        states.append(np.stack((mass, mass*uf, mass*vf, mass*fraction), axis=-1))
        velocities.append(uf if axis == 1 else vf)
        gammas.append(rho*(p.g*(rho-p.rho_a)/rho))
        densities.append(rho)
    qL, qR = states
    uL, uR = velocities
    gL, gR = gammas
    # Speeds use the selected depth, including its dry zero-velocity convention.
    rhoL, rhoR = densities
    cL = np.sqrt(np.maximum(Gface*gL/rhoL*hL, 0))
    cR = np.sqrt(np.maximum(Gface*gR/rhoR*hR, 0))
    am, ap = np.minimum(np.minimum(uL-cL, uR-cR), 0), np.maximum(np.maximum(uL+cL, uR+cR), 0)
    den = ap-am
    safe = np.where(den > p.speed_eps, den, 1)
    transport = (ap[..., None]*(uL-am)[..., None]*qL-am[..., None]*(uR-ap)[..., None]*qR)/safe[..., None]
    transport = np.where((den > p.speed_eps)[..., None], transport, 0.5*(uL[...,None]*qL+uR[...,None]*qR))
    limited_component = ((transport[...,0] > 0) & (transport[...,3] > transport[...,0])) | ((transport[...,0] < 0) & (transport[...,3] < transport[...,0]))
    transport[...,3] = np.where(limited_component,transport[...,0],transport[...,3])
    donorL, donorR = face_pair(h, h, axis)
    blocked = ((transport[...,0] > 0) & ((donorL <= p.dry_h_tol) | (hL <= p.dry_h_tol))) | ((transport[...,0] < 0) & ((donorR <= p.dry_h_tol) | (hR <= p.dry_h_tol)))
    transport[blocked] = 0
    rest = (uL == 0) & (uR == 0) & (abs(eL-eR) <= 8*EPS*np.maximum(1,np.maximum(abs(eL),abs(eR))))
    transport[...,0] = np.where(rest, 0, transport[...,0])
    transport[...,3] = np.where(rest, 0, transport[...,3])
    normal = 1 if axis == 1 else 2
    # Outflow-only physical boundaries remove mass-carried transport, retaining normal momentum.
    for index, sign in ((0, 1), (-1, -1)):
        section = [slice(None)]*2;section[axis] = index;section=tuple(section)
        incoming = sign*transport[section][...,0] > 0
        for component in (0, 3, 3-normal):
            plane = transport[...,component][section]
            plane[incoming] = 0
    # Face path uses the selected transport depths with independent eta_x levels.
    hpathL, hpathR = hL, hR
    epathL, epathR = np.where(hL <= p.dry_h_tol,bface,eL), np.where(hR <= p.dry_h_tol,bface,eR)
    hpathL, hpathR = np.where(hL <= p.dry_h_tol,0,hpathL), np.where(hR <= p.dry_h_tol,0,hpathR)
    path = hydrostatic_path(hpathL,hpathR,gL,gR,epathL,epathR,Gface,Gface,p)
    left, right = transport.copy(), transport.copy()
    left[...,normal] += np.where(hL > p.dry_h_tol,am/safe*path,0)
    right[...,normal] += np.where(hR > p.dry_h_tol,ap/safe*path,0)
    left[...,normal] = np.where(den > p.speed_eps,left[...,normal],transport[...,normal])
    right[...,normal] = np.where(den > p.speed_eps,right[...,normal],transport[...,normal])
    # Cell-internal path uses independent eta_x endpoints and directional face gravity.
    gm = np.take(gR,range(h.shape[axis]),axis=axis)
    gp = np.take(gL,range(1,h.shape[axis]+1),axis=axis)
    cellpath = hydrostatic_path(np.maximum(em-bm,0),np.maximum(ep-bp,0),gm,gp,em,ep,
                              np.take(Gface,range(h.shape[axis]),axis=axis),np.take(Gface,range(1,h.shape[axis]+1),axis=axis),p)
    cellpath = np.where(h > p.dry_h_tol,cellpath,0)
    rhs = -(np.take(left,range(1,h.shape[axis]+1),axis=axis)-np.take(right,range(h.shape[axis]),axis=axis))/spacing
    rhs[...,normal] += cellpath/spacing
    return rhs, dict(a_minus=am,a_plus=ap,hL=hL,hR=hR,etaL=eL,etaR=eR,
                     mass_flux=transport[...,0],
                     rest_faces=np.all(qL[...,1:3] == 0,axis=-1) & np.all(qR[...,1:3] == 0,axis=-1),
                     arithmetic_scale=np.maximum(abs(Gface*gL)*hL*np.maximum(hL,abs(eL)),
                                                 abs(Gface*gR)*hR*np.maximum(hR,abs(eR))))


def rhs(q, vertices, dx, dy, p, return_diagnostics=False):
    q = np.asarray(q,float)
    if q.ndim != 3 or q.shape[-1] != 4 or not np.all(np.isfinite(q)) or np.any(q[...,0] < 0):
        raise ValueError('invalid raw conservative state')
    bed = VertexBed(vertices)
    sg = geometry.build_slope_geometry_2d(bed,dx,dy) if p.slope_correction or p.curvature_term else None
    if sg is not None:
        # geometry_2d copies both derivative boundary strips after its LS fit.
        derivatives=[np.array(getattr(sg,key),copy=True) for key in ('Bx_center','By_center','Bxx_center','Bxy_center','Byy_center')]
        for field in derivatives:
            if q.shape[1] >= 5:field[:,:2]=field[:,2:3];field[:,-2:]=field[:,-3:-2]
            if q.shape[0] >= 5:field[:2]=field[2:3];field[-2:]=field[-3:-2]
        bx,by,bxx,bxy,byy=derivatives
        G=1/(1+bx*bx+by*by)
        Gx=np.column_stack((G[:,0],0.5*(G[:,:-1]+G[:,1:]),G[:,-1]))
        Gy=np.vstack((G[0],0.5*(G[:-1]+G[1:]),G[-1]))
    else:
        G=np.ones(q.shape[:2]);Gx=np.ones(bed.Bx_face.shape);Gy=np.ones(bed.By_face.shape)
    result=np.zeros_like(q);diagnostics={}
    for axis,spacing,gravity in ((1,dx,Gx if p.slope_correction else np.ones_like(Gx)),(0,dy,Gy if p.slope_correction else np.ones_like(Gy))):
        if q.shape[axis] <= 1:continue
        r,d=direction(q,bed,spacing,axis,p,gravity);result+=r;diagnostics['x' if axis==1 else 'y']=d
    # Remove only pressure-scaled arithmetic noise in an exactly resting stencil.
    resting = np.all(q[...,1:3] == 0,axis=-1)
    for name,axis in (('x',1),('y',0)):
        if name not in diagnostics:
            continue
        face_rest = diagnostics[name]['rest_faces']
        resting &= np.take(face_rest,range(q.shape[axis]),axis=axis)
        resting &= np.take(face_rest,range(1,q.shape[axis]+1),axis=axis)
    for name,axis,component,spacing in (('x',1,1,dx),('y',0,2,dy)):
        if name not in diagnostics:
            continue
        scale = diagnostics[name]['arithmetic_scale']
        scale = np.maximum(np.take(scale,range(q.shape[axis]),axis=axis),
                           np.take(scale,range(1,q.shape[axis]+1),axis=axis))
        quiet = resting & (abs(result[...,component]) <= 64*EPS*scale/spacing)
        result[...,component][quiet] = 0
    if p.curvature_term:
        h,u,v,_=primitive(q,p)
        acceleration=bxx*u*u+2*bxy*u*v+byy*v*v
        scale=np.where(h > p.dry_h_tol,-q[...,0]*acceleration*(G if p.slope_correction else 1),0)
        result[...,1]+=scale*bx;result[...,2]+=scale*by
    return (result,diagnostics) if return_diagnostics else result


def advance(initial, vertices, dx, dy, dt_values, p, return_diagnostics=False):
    """Current explicit part of the Fortran n_RK=2 IMEX tableau: forward Euler."""
    q=np.asarray(initial,float).copy();minimum=float(np.min(q[...,0]/p.rho_c))
    boundary_transfer = 0.0
    for dt in dt_values:
        residual, diagnostics = rhs(q,vertices,dx,dy,p,return_diagnostics=True)
        boundary_rate = 0.0
        if 'x' in diagnostics:
            flux = diagnostics['x']['mass_flux']
            boundary_rate -= dy*(np.sum(flux[:,-1])-np.sum(flux[:,0]))
        if 'y' in diagnostics:
            flux = diagnostics['y']['mass_flux']
            boundary_rate -= dx*(np.sum(flux[-1])-np.sum(flux[0]))
        boundary_transfer += dt*boundary_rate
        q += dt*residual
        if not np.all(np.isfinite(q)) or np.min(q[...,0]) < -EPS:
            raise AssertionError('negative/nonfinite raw reference final state')
        raw_h = (q[...,0]-q[...,3])/p.rho_c+q[...,3]/p.rho_s
        minimum = min(minimum,float(np.min(raw_h)))
        tolerance = 4096*EPS*np.maximum(abs(q[...,0]),np.finfo(float).tiny)
        small_negative = (q[...,3] < 0) & (q[...,3] >= -tolerance)
        q[...,3][small_negative] = 0
        q[q[...,0] < EPS]=0
        h,*_=primitive(q,p)
        q[...,1:3][h <= p.dry_h_tol]=0
    if return_diagnostics:
        return q,minimum,dict(boundary_mass_transfer=float(boundary_transfer))
    return q,minimum
