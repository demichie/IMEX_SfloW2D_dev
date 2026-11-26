#!/usr/bin/env python3
"""Enforce all-step same-input KT checks and retain independent trajectory drift."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from compare_fields import CONTRACT, common_input, accepted_times, kt



def check_same_input_step(previous, observed, vertices, dx, dt, params, scales, epsilon):
    """Check every conservative component after one update from the same state."""
    if previous.ndim != 3 or previous.shape[-1] != 4 or observed.shape != previous.shape[:2]+(3,):
        raise AssertionError("invalid same-input conservative field shapes")
    if not np.all(np.isfinite(previous)) or not np.all(np.isfinite(observed)):
        raise AssertionError("nonfinite same-input conservative field")
    predicted, _ = kt.advance(previous, vertices, dx, dx, [dt], params)
    difference = observed-predicted[..., :3]
    scaled = np.max(np.abs(difference), axis=(0, 1))/scales
    if not np.all(np.isfinite(scaled)):
        raise AssertionError("nonfinite same-input error")
    component = int(np.argmax(scaled))
    j, i = np.unravel_index(np.argmax(np.abs(difference[..., component])), difference.shape[:2])
    result = {"max_scaled_error": scaled.tolist(),
              "worst_component": component, "worst_cell_ij": [int(i), int(j)],
              "fortran_q": observed[j, i].tolist(),
              "python_one_step_q": predicted[j, i, :3].tolist(),
              "signed_error": difference[j, i].tolist()}
    return result, bool(np.all(scaled <= epsilon))


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
    max_same_input=np.zeros(3)
    first_operator_failure=None
    worst_operator_step=None
    for index, step_dt in enumerate(dt):
        previous=initial.copy() if index==0 else np.zeros_like(initial)
        if index:
            previous[...,:3]=states[index-1]
        direct, passed=check_same_input_step(previous,states[index],vertices,dx,
                                             step_dt,params,scales,eps)
        max_same_input=np.maximum(max_same_input,np.asarray(direct["max_scaled_error"]))
        if worst_operator_step is None or max(direct["max_scaled_error"]) > max(worst_operator_step["detail"]["max_scaled_error"]):
            worst_operator_step={"step":index+1,"detail":direct}
        if not passed and first_operator_failure is None:
            first_operator_failure={"step":index+1,"start_time":float(np.sum(dt[:index])),
                                    "dt":float(step_dt),"detail":direct,"limit":float(eps)}
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
    report={"mode":mode,"dx":dx,"steps":len(dt),
            "first_out_of_contract_step":first,
            "max_scaled_errors_over_trajectory":max_errors.tolist(),
            "same_input_all_steps": {
                "status": "passed" if first_operator_failure is None else "failed",
                "steps_checked":len(dt),
                "scaled_limit_per_step":float(eps),
                "maximum_scaled_errors":max_same_input.tolist(),
                "worst_step":worst_operator_step,
                "first_failure":first_operator_failure,
            },
            "accepted_states_sha256":evidence["accepted_states_sha256"],
            "note":"All comparisons use Fortran binary64 states, not formatted ASCII output."}
    (case_dir/"trajectory_evidence.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"TRAJECTORY {mode} dx={dx:g} first out-of-contract={first}",flush=True)
    print(f"TRAJECTORY largest scaled errors={max_errors.tolist()}",flush=True)
    print(f"SAME-INPUT ALL-STEPS {mode} dx={dx:g} "
          f"status={report['same_input_all_steps']['status']} "
          f"maximum_scaled_errors={max_same_input.tolist()} limit={eps:.17g}",flush=True)
    if first_operator_failure is not None:
        raise AssertionError(f"Gate-H same-input operator mismatch: {first_operator_failure}")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("case_dir",type=Path)
    cli.add_argument("mode",choices=CONTRACT["modes"])
    args=cli.parse_args()
    diagnose(args.case_dir,args.mode)
