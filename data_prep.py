import numpy as np
from tqdm.auto import trange

from niftypet import nimpa, nipet

try:
    import cuvec as cu
except ImportError:
    cu = None


def prepare_data(folderin, opth, cuvec=True):
    if cuvec:
        assert cu is not None
    mMRpars = nipet.get_mmrparams()                  # scanner parameters
    datain = nipet.classify_input(folderin, mMRpars) # categorise input data

    # hardware and object mu-maps
    muh = nipet.hdw_mumap(datain, [1, 2, 4], mMRpars, outpath=opth, use_stored=True)
    muo = nipet.obj_mumap(datain, mMRpars, outpath=opth, store=True)
    mu = muh['im'] + muo['im'] # combine mu-maps

    # histogram into measured sinogram
    hst = nipet.mmrhist(datain, mMRpars, outpath=opth, store=True, use_stored=True)
    m = nipet.mmraux.remgaps(hst['psino'], mMRpars['txLUT'], mMRpars['Cnt']).astype(np.float32)
    if cuvec:
        m = cu.asarray(m)

    # obtain the attenuation and normalisation sinograms
    A = nipet.frwd_prj(mu, mMRpars, attenuation=True, dev_out=True)
    N = nipet.mmrnorm.get_sinog(datain, hst, mMRpars['axLUT'], mMRpars['txLUT'], mMRpars['Cnt'])

    # combine attenuation and normalisation (diagonal of AN matrix)
    AN = A * N
    if cuvec:
        AN = cu.asarray(AN)

    # perform default reconstruction to obtain scatter
    rec_ = nipet.mmrchain(datain, mMRpars, mu_h=muh, mu_o=muo, itr=2, histo=hst, outpath=opth,
                          ret_sinos=True)
    eim = rec_['im']
    s = rec_['sinos']['ssino'] # scatter
    r = rec_['sinos']['rsino'] # randoms

    # normalise in GPU dimensions
    r = nipet.mmraux.remgaps(r, mMRpars['txLUT'], mMRpars['Cnt'])
    s = nipet.mmraux.remgaps(s, mMRpars['txLUT'], mMRpars['Cnt'])
    print("Randoms: %.3g%%" % (r.sum() / m.sum() * 100))
    print("Scatter: %.3g%%" % (s.sum() / m.sum() * 100))

    # image FoV mask
    msk = nipet.img.mmrimg.get_cylinder(mMRpars['Cnt'], rad=29., xo=0., yo=0., unival=1,
                                        gpu_dim=True) > 0.9

    # subsets and sensitivity image for each subset
    Sn = 14                   # number of subsets
    sinoTIdx = [None] * Sn    # subset indices
    sen = [None] * Sn         # sensitivity images for each subset
    sen_inv_msk = [None] * Sn # masked inverse sensitivity image

    # obtain sensitivity image for each subset, inverted and masked
    for n in trange(Sn, unit="subset"):
        sinoTIdx[n] = cu.asarray(nipet.prj.mmrrec.get_subsets14(n, mMRpars)[0], np.int32)
        if cuvec:
            sinoTIdx[n] = cu.asarray(sinoTIdx[n], np.int32)
            sen[n] = nipet.back_prj(nimpa.isub(AN, sinoTIdx[n], sync=False), mMRpars,
                                    isub=sinoTIdx[n], dev_out=True)
            sen_inv_msk[n] = cu.zeros_like(sen[n])
        else:
            sen[n] = nipet.back_prj(AN[sinoTIdx[n]], mMRpars, isub=sinoTIdx[n], dev_out=True)
            sen_inv_msk[n] = np.zeros_like(sen[n])
        sen_inv_msk[n][msk] = np.float32(1) / sen[n][msk]
        assert sen[n].shape == (mMRpars['Cnt']['SZ_IMX'], mMRpars['Cnt']['SZ_IMY'],
                                mMRpars['Cnt']['SZ_IMZ'])

    # get the measured and randoms/scatter data in subset form
    if cuvec:
        rs_AN = nimpa.div(r + s, AN, default=1)
        m_sub = [nimpa.isub(m, idx) for idx in sinoTIdx]
        rs_AN_sub = [nimpa.isub(rs_AN, idx) for idx in sinoTIdx]
    else:
        rs_AN = (r+s) / AN
        m_sub = [m[idx, :] for idx in sinoTIdx]
        rs_AN_sub = [rs_AN[idx, :] for idx in sinoTIdx]

    return {
        'params': mMRpars, 'datain': datain, 'mu': mu, 'msub': m_sub, 'rsn_sub': rs_AN_sub,
        'msk': msk, 'isen': sen_inv_msk, 'sidx': sinoTIdx, 'senim': sen, 'Sn': Sn, 'eim': eim}
