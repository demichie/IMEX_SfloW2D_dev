#!/usr/bin/env python3
"""Find the first full-precision conservative state divergence in Gate-H."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from compare_fields import CONTRACT, common_input, accepted_times, kt


def diagnose(case_dir, mode):
    geometry = np.load(case_dir/"gate_h_geometry.npz")
    shape = geometry["bed"].shape
    dx = float(geometry["dx"])
    evidence = json.loads((case_dir/"accepted_time_evidence.json").read_text())
    states_path = case_dir/"accepted_states.bin"
    if "accepted_states_sha256" not in evidence:
        raise AssertionError("observer missing accepted-state stream attestation")
    if hashlib.sha256(states_path.read_bytes()).hexdigest() != evidence["accepted_states_sha256"]:
        raise AssertionError("accepted state stream hash mismatch")
    dt = accepted_times(np.fromfile(case_dir/"accepted_times.bin",dtype=np.float64).reshape(-1,2),
                        CONTRACT["time"])
    values = np.fromfile(states_path,dtype=np.float64)
    expected_size = len(dt)*shape[0]*shape[1]*3
    if len(values) != expected_size:
        raise AssertionError(f"state stream has {len(values)} values, expected {expected_size}")
    states = values.reshape((len(dt),*shape,3))
    initial, vertices = common_input(case_dir,shape,dx)
    settings=CONTRACT["current_IMEX"]
    flags=mode=="slope_curvature"
    params=kt.Parameters(rho_a=101300/(287.051*300),eps_sing=settings["eps_sing"],
                         reconstruction=settings["KT_TRACE_RECONSTRUCTION"],
                         slope_correction=flags,curvature_term=flags)
    mass_scale=max(1,float(initial[...,0].max()))
    momentum_scale=mass_scale*np.sqrt(params.g*max(1,float(initial[...,0].max())/params.rho_c))
    scales=np.array([mass_scale,momentum_scale,momentum_scale])
    eps=settings["roundoff_epsilon_multiplier"]*np.finfo(float).eps
    reference=initial.copy()
    first=None
    max_errors=np.zeros(3)
    for index, step_dt in enumerate(dt):
        reference,_,_=kt.advance(reference,vertices,dx,dx,[step_dt],params,
                                 return_diagnostics=True)
        difference=states[index]-reference[...,:3]
        scaled=np.max(np.abs(difference),axis=(0,1))/scales
        max_errors=np.maximum(max_errors,scaled)
        if first is None and np.any(scaled > eps*(index+1)):
            component=int(np.argmax(scaled))
            ij=np.unravel_index(np.argmax(np.abs(difference[...,component])),shape)
            first={
                "step_number":index+1, "start_time":float(np.sum(dt[:index])),
                "dt":float(step_dt), "cell_ij":[int(ij[1]),int(ij[0])],
                "component_index":component,
                "scaled_error":scaled.tolist(), "limit":float(eps*(index+1)),
                "fortran_q":states[index][ij].tolist(),
                "python_q":reference[ij][:3].tolist(),
                "signed_difference":difference[ij].tolist()
            }
    # Re-evaluate the first discrepant step from the *same* Fortran
    # pre-step state. This removes pre-existing trajectory drift and
    # exposes a genuinely different spatial update at that step.
    if first is not None:
        step=first["step_number"]
        previous=initial.copy() if step==1 else np.zeros_like(initial)
        if step>1:
            previous[...,:3]=states[step-2]
            previous[...,3]=0.0  # Gate-H has no solid species.
        predicted,_,_=kt.advance(previous,vertices,dx,dx,[dt[step-1]],params,
                                  return_diagnostics=True)
        local_error=states[step-1]-predicted[...,:3]
        local_scaled=np.max(np.abs(local_error),axis=(0,1))/scales
        component=int(np.argmax(local_scaled))
        j,i=np.unravel_index(np.argmax(np.abs(local_error[...,component])),shape)
        jj=slice(max(0,j-1),min(shape[0],j+2))
        ii=slice(max(0,i-1),min(shape[1],i+2))
        first["same_input_one_step"]={
            "max_scaled_error":local_scaled.tolist(),
            "component_index":component,
            "worst_cell_ij":[int(i),int(j)],
            "fortran_previous_q":previous[j,i,:3].tolist(),
            "fortran_accepted_q":states[step-1,j,i,:].tolist(),
            "python_update_from_fortran_q":predicted[j,i,:3].tolist(),
            "fortran_applied_rate":((states[step-1,j,i,:]-previous[j,i,:3])/dt[step-1]).tolist(),
            "python_rate":((predicted[j,i,:3]-previous[j,i,:3])/dt[step-1]).tolist(),
            "mass_error_patch":local_error[jj,ii,0].tolist(),
            "momentum_x_error_patch":local_error[jj,ii,1].tolist(),
            "momentum_y_error_patch":local_error[jj,ii,2].tolist(),
            "previous_mass_patch":previous[jj,ii,0].tolist(),
        }
        print(f"TRAJECTORY SAME-INPUT {mode} dx={dx:g} "
              f"step={step} worst={(i,j)} scaled={local_scaled.tolist()}",flush=True)
        print(f"TRAJECTORY SAME-INPUT local={first['same_input_one_step']}",flush=True)
    report={"mode":mode,"dx":dx,"steps":len(dt),
            "first_out_of_contract_step":first,
            "max_scaled_errors_over_trajectory":max_errors.tolist(),
            "accepted_states_sha256":evidence["accepted_states_sha256"],
            "note":"All comparisons use Fortran binary64 states, not formatted ASCII output."}
    (case_dir/"trajectory_evidence.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"TRAJECTORY {mode} dx={dx:g} first out-of-contract={first}",flush=True)
    print(f"TRAJECTORY largest scaled errors={max_errors.tolist()}",flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("case_dir",type=Path)
    cli.add_argument("mode",choices=CONTRACT["modes"])
    args=cli.parse_args()
    diagnose(args.case_dir,args.mode)
