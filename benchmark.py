"""Usage:
  benchmark [options] <folderin>

Options:
  -c, --cuvec  : Whether to run new CuVec version.

Arguments:
  <folderin>  : Input directory [default: Ab_PET_mMR_test].
"""
import logging
from pathlib import Path
from time import time

import numpy as np
from argopt import argopt
from tqdm.auto import trange

from data_prep import cu, prepare_data
from niftypet import nimpa, nipet


def run_noncu(recdat, nitr=4):
    # reconstructed image
    y = np.ones_like(recdat['isen'][0])
    tic = time()
    for _ in trange(nitr, desc="OSEM"):
        for n in trange(recdat['Sn'], unit="subset", leave=False):
            # forward-project current reconstruction to sinogram space
            Xy = nipet.frwd_prj(y, recdat['params'], isub=recdat['sidx'][n], attenuation=False,
                                dev_out=True, fullsino_out=False)
            # add the normalised sum of randoms and scatter
            Xy += recdat['rsn_sub'][n]
            # corrections
            crr = recdat['msub'][n] / Xy
            # back-project corrections to image space
            bim = nipet.back_prj(crr, recdat['params'], isub=recdat['sidx'][n], dev_out=True)
            # update image
            y *= bim * recdat['isen'][n]
    return y, time() - tic


def run_cuvec(recdat, nitr=4, sync=False):
    # reconstructed image
    y = cu.ones_like(recdat['isen'][0])
    tic = time()
    # temporary variables
    Xy, crr, bim, mul = (None,) * 4
    for _ in trange(nitr, desc="OSEM"):
        for n in trange(recdat['Sn'], unit="subset", leave=False):
            # forward-project current reconstruction to sinogram space
            Xy = nipet.frwd_prj(y, recdat['params'], isub=recdat['sidx'][n], attenuation=False,
                                dev_out=True, fullsino_out=False, output=Xy, sync=sync)
            # add the normalised sum of randoms and scatter
            Xy = nimpa.add(Xy, recdat['rsn_sub'][n], output=Xy, sync=sync)
            # corrections
            crr = nimpa.div(recdat['msub'][n], Xy, default=1, output=crr, sync=sync)
            # back-project corrections to image space
            bim = nipet.back_prj(crr, recdat['params'], isub=recdat['sidx'][n], dev_out=True,
                                 output=bim, sync=sync)
            # apply FOV mask and normalise with scanner sensitivity
            mul = nimpa.mul(bim, recdat['isen'][n], output=mul, sync=sync)
            # update reconstructed image
            y = nimpa.mul(y, mul, output=y, sync=sync)
    cu.dev_sync()
    return y, time() - tic


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    args = argopt(__doc__).parse_args()
    folderin = Path(args.folderin).expanduser() # input path

    recdat = prepare_data(folderin, folderin / ('_output_cuvec' if args.cuvec else '_output_b4cu'),
                          cuvec=args.cuvec)
    # import matplotlib.pyplot as plt
    # plt.matshow(np.sum(recdat['isen'], axis=(0,3)))
    # plt.matshow(recdat['eim'][60,...])
    # plt.matshow(recdat['mu'][60,...])
    y, dt = (run_cuvec if args.cuvec else run_noncu)(recdat)
    print(f"Time elapsed:{dt:.3f}s")
    # import cProfile
    # cProfile.run('ycu, dt = run_cuvec()')
    # cProfile.run('y, dt = run_noncu()')
