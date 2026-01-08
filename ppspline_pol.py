#!/usr/bin/env python

############
# ppspline #
############

# ppspline is a command-line program to make a frequency-parameterized model of
#    wideband profile evolution.  The parameterization is good within the range
#    of data frequencies, provided there are not huge gaps in frequency.  For
#    an input nchan x nbin average portrait of aligned profiles, the model is a
#    B-spline representation of the curve traced out by the nchan profile
#    amplitude vectors in an nbin vector space.  Since the profile shapes are
#    highly correlated, the B-spline representation can be reduced to << nbin
#    dimensions, and is limited to ten dimensions.  Therefore, the profile
#    variability is decomposed using PCA and a small number of eigenprofiles
#    that encompass most of the profile evolution are selected.  The input data
#    should be high S/N, averaged, and "aligned" (e.g. output from
#    ppalign.py).  Pre-normalization of the input is encouraged (specifically
#    using normalize_portrait('prof')), as is using smooth=True in the
#    make_model function.

# Written by Timothy T. Pennucci (TTP; tim.pennucci@nanograv.org).
# Updated by Gosia Curylo (gosia.curylo@monash.edu).

from __future__ import division
from __future__ import print_function

from past.utils import old_div
from pplib_pol import *


class DataPortrait(DataPortrait):
    """
    DataPortrait is a class that contains the data to which a model is fit.

    This class adds methods and attributes to the parent class specific to
        modeling profile evolution with a B-spline curve.
    """

    def make_spline_model(self, state='Intensity', max_ncomp=10, smooth=True, snr_cutoff=150.0,
                          rchi2_tol=0.1, k=3, sfac=1.0, max_nbreak=None, model_name=None,
                          quiet=False, **kwargs):
        """
        Make a model based on PCA and B-spline interpolation.

        state is the polarization state; either 'Intensity' (default), 'Stokes' (I, Q, U, V) or 'Stokes_L' (I, L, V).
        max_ncomp is the maximum number of PCA components to use in the
            B-spline parameterization; max_ncomp <= 10.
            Can be different for each polarisation, use a list with npol params.
        smooth=True will smooth the eigenvectors and mean profile using
            a reduced chi-squared figure-of-merit.
        snr_cutoff is the S/N ratio value above or equal to which an
            eigenvector is deemed "significant".  Setting it equal to np.inf
            would ensure only a mean profile model is returned.
        rchi2_tol is the tolerance parameter that will allow greater deviations
            in the smoothed profile from the input profiles' shapes.
        k is the polynomial degree of the spline; cubic splines (k=3)
            recommended; 1 <= k <= 5.  NB: polynomial order = degree + 1.
        sfac is a multiplicative smoothing factor passed to si.splprep; greater
            values result in more smoothing.  sfac=0 will make an interpolating
            model anchored on the input data profiles.
        max_nbreak is the maximum number of breakpoints (unique knots) to
            allow.  If provided, this may override sfac and enforce smoothing
            based on max_nbreak breakpoints.  That is, if the fit returns n >
            max_nbreak breakpoints, it will refit using maximum max_nbreak
            breakpoints, irrespective of the other smoothing condition.  The
            corresponding maximum number of B-splines will be max_nspline =
            max_nbreak + k - 1.  max_nbreak should be >= 2.
        model_name is the name of the model; defaults to self.source +
            '.spl'
        quiet=True suppresses output.
        **kwargs get passed to find_significant_eigvec(...).
        """

        if state=='Stokes' or state=='Stokes_L':
            print("********************************************")
            print("Let's begin to make the spline models...")
            print("********************************************")
            print("")
            states = ['I', 'Q', 'U', 'V'] if state=='Stokes' else ['I', 'L', 'V']
            nsub   = self.nsub
            npol   = self.npol if state=='Stokes' else 3
            nchan  = self.nchan
            nchanx = len(self.freqsxs[0])
            nbin   = self.nbin

            # Initiate new attributes
            list_attributes = ['ieig', 'proj_port', 'model_name', 'tck', 'u', 'fp', 'ier','msg']  #these need to be 2D lists of arrays
            for attr in list_attributes:
                setattr(self, attr, [[None for _ in range(npol)] for _ in range(nsub)])
            
            self.ncomp     = np.zeros([nsub,npol]) 
            self.eigvec    = np.zeros([nsub,npol,nbin,nbin])
            self.eigval    = np.zeros([nsub,npol,nbin])
            self.mean_prof = np.zeros([nsub,npol,nbin])
            if smooth:
                self.smooth_mean_prof = np.zeros([nsub,npol,nbin])
                self.smooth_eigvec    = np.zeros([nsub,npol,nbin,nbin])
            self.reconst_port = np.zeros([nsub,npol,nchanx, nbin])
            
            self.model        = np.zeros([nsub,npol,nchan,nbin])
            self.modelx       = np.zeros([nsub,npol,nchanx,nbin])
            self.model_masked = np.zeros([nsub,npol,nchan,nbin])

            #Same for each polarization
            freqs = self.freqsxs[0]
            nu_lo = freqs.min()
            nu_hi = freqs.max()

            #Save the user defined model params
            max_ncomp0  = max_ncomp
            snr_cutoff0 = snr_cutoff
            rchi2_tol0  = rchi2_tol
            sfac0       = sfac
            max_nbreak0 = max_nbreak
            k0          = k  
            for isub in range(nsub):
                for ipol in range(npol):

                    #Modeling params for each polarisation state
                    max_ncomp  = max_ncomp0[ipol] if isinstance(max_ncomp0, list) else max_ncomp0
                    snr_cutoff = snr_cutoff0[ipol] if isinstance(snr_cutoff0, list) else snr_cutoff0
                    rchi2_tol  = rchi2_tol0[ipol] if isinstance(rchi2_tol0, list) else rchi2_tol0
                    sfac       = sfac0[ipol] if isinstance(sfac0, list) else sfac0
                    max_nbreak = max_nbreak0[ipol] if isinstance(max_nbreak0, list) else max_nbreak0
                    k          = k0[ipol] if isinstance(k0, list) else k0
                    
                    print("")
                    print("Model parameters for state %s:"%states[ipol])
                    print("Max_ncomp = %s"%max_ncomp)
                    print("Snr_cutoff = %s"%snr_cutoff)
                    print("rchi2_tol = %s"%rchi2_tol)
                    print("Sfac = %s"%sfac)
                    print("Max_nbreak = %s"%max_nbreak)
                    print("k = %s"%k)
                    print("")

                    # Definitions
                    port        = self.portx[isub,ipol] if state=='Stokes' else self.portLx[isub,ipol]      #np.shape(nchanx,nbin), same I and V
                    snrs        = self.SNRsxs[isub,ipol] if state=='Stokes' else self.LSNRsxs[isub,ipol]    
                    noisexs     = self.noise_stdsxs[isub,ipol] if state=='Stokes' else self.Lnoise_stdsxs[isub,ipol]
                    pca_weights = old_div(snrs, np.sum(snrs))                                               #np.shape(nchanx)
                    mean_prof   = old_div((port.T * pca_weights).T.sum(axis=0), pca_weights.sum())          #np.shape(nbin)
                
                    # Check nbin
                    nbin = port.shape[1]
                    if nbin % 2 != 0:
                        print("nbin = %d is odd; cannot wavelet_smooth.\n" % nbin)
                        smooth = False

                    elif np.modf(np.log2(nbin))[0] != 0.0:
                        print(
                            "nbin = %d is not a power of two; can only try wavelet_smooth to one level; recommend resampling to a power-of-two number of phase bins.\n" % nbin)
                    
                    # Do principal component analysis
                    eigval, eigvec = pca(port, mean_prof, pca_weights, quiet=quiet)

                    # Get "significant" eigenvectors
                    if max_ncomp is None:
                        return_max = 10
                    else:
                        return_max = min(max_ncomp, 10)

                    if smooth:
                        if 'pywt' not in sys.modules:
                            raise ImportError("You failed to import pywt and need PyWavelets to use smooth=True!")
                        ieig, smooth_eigvec = find_significant_eigvec(eigvec, check_max=10,
                                                                    return_max=return_max, snr_cutoff=snr_cutoff,
                                                                    return_smooth=True, rchi2_tol=rchi2_tol, **kwargs)

                    else:
                        ieig = find_significant_eigvec(eigvec, check_max=10,
                                                    return_max=return_max, snr_cutoff=snr_cutoff,
                                                    return_smooth=False, rchi2_tol=rchi2_tol, **kwargs)
                    ncomp = len(ieig)

                    if smooth:
                        smooth_mean_prof = smart_smooth(mean_prof,
                                                        rchi2_tol=rchi2_tol)

                    if ncomp == 0:  # Will make model with constant average port
                        proj_port = port[:, :ncomp]
                        if smooth:
                            modelx = reconst_port = np.tile(smooth_mean_prof,
                                                            len(freqs)).reshape(len(freqs), port.shape[1])
                            model = np.tile(smooth_mean_prof,
                                            len(self.freqs[0])).reshape(len(self.freqs[0]),
                                                                        port.shape[1])
                        else:
                            modelx = reconst_port = np.tile(mean_prof,
                                                            len(freqs)).reshape(len(freqs), port.shape[1])
                            model = np.tile(mean_prof,
                                            len(self.freqs[0])).reshape(len(self.freqs[0]),
                                                                        port.shape[1])
                    else:
                        delta_port = port - mean_prof
                        if smooth:
                            reconst_port = reconstruct_portrait(port, mean_prof,
                                                                smooth_eigvec[:, ieig])
                            # Find the projections of the profiles onto the basis components
                            proj_port = np.dot(delta_port, smooth_eigvec[:, ieig])
                        else:
                            reconst_port = reconstruct_portrait(port, mean_prof,
                                                                eigvec[:, ieig])
                            # Find the projections of the profiles onto the basis components
                            proj_port = np.dot(delta_port, eigvec[:, ieig])

                    if ncomp == 0:
                        (tck, u) = [np.array([]), np.array([]), 0], np.array([])
                        fp, ier, msg = None, None, None
                    else:

                        spl_weights = pca_weights
                        #s = sfac * len(proj_port) * \
                        #        np.sum((self.SNRsxs[isub,ipol] * self.noise_stdsxs[isub,ipol])**2) / \
                        #        sum(self.SNRsxs[isub,ipol])**2
                        s = sfac * len(proj_port) * \
                                np.sum((snrs * noisexs)**2) / \
                                sum(snrs)**2
                        
                        if self.bw < 0: flip = -1   #u in si.splprep has to be increasing...
                        else: flip = 1

                        #Find the B-spline curve traced by the projected vectors,
                        #parameterized by frequency
                        (tck,u), fp, ier, msg = si.splprep(proj_port[::flip].T,
                                w=spl_weights[::flip], u=freqs[::flip], ub=nu_lo, ue=nu_hi,
                                k=k, task=0, s=s, t=None, full_output=1, nest=None, per=0,
                                quiet=int(quiet))

                        if max_nbreak is not None and len(np.unique(tck[0])) > max_nbreak:
                            if max_nbreak < 2:
                                print("max_nbreak not >= 2; setting max_nbreak = 2...")
                                max_nbreak = 2
                            if max_nbreak == 2: s = np.inf
                            (tck, u), fp, ier, msg = si.splprep(proj_port[::flip].T,
                                                                w=spl_weights[::flip], u=freqs[::flip], ub=nu_lo,
                                                                ue=nu_hi, k=k, task=0, s=s, t=None, full_output=1,
                                                                nest=max_nbreak + (k * 2), per=0, quiet=int(quiet))

                        if ier > 1:  # Will also catch when ier == "unknown"
                            print("Something went wrong in si.splprep for %s:\n%s" % (
                                self.source, msg))
                    # Build model
                    if ncomp != 0:
                        if smooth:
                            modelx = gen_spline_portrait(smooth_mean_prof, freqs,
                                                        smooth_eigvec[:, ieig], tck)
                            model = gen_spline_portrait(smooth_mean_prof, self.freqs[0],
                                                        smooth_eigvec[:, ieig], tck)
                        else:
                            modelx = gen_spline_portrait(mean_prof, freqs, eigvec[:, ieig],
                                                        tck)
                            model = gen_spline_portrait(mean_prof, self.freqs[0],
                                                        eigvec[:, ieig], tck)

                    # Assign new attributes
                    self.ieig[isub][ipol] = ieig
                    self.ncomp[isub,ipol] = ncomp
                    self.eigvec[isub,ipol] = eigvec
                    self.eigval[isub,ipol] = eigval
                    self.mean_prof[isub,ipol] = mean_prof
                    if smooth:
                        self.smooth_mean_prof[isub,ipol] = smooth_mean_prof
                        self.smooth_eigvec[isub,ipol] = smooth_eigvec
                    self.proj_port[isub][ipol] = proj_port
                    self.reconst_port[isub,ipol] = reconst_port
                    # tck contains the knot locations t, B-spline coefficients c, and
                    # polynomial degree k -- end knots will have multiplicity k+1, interior
                    # breakpoints will have multiplicity 1 for maximum continuity.  The
                    # number of B-splines will be n = l + k, where l is the number of
                    # intervals.  l = number of breakpoints - 1 = number of unique knots - 1
                    # = len(tck[0]) - 2*tck[2] - 1.
                    self.tck[isub][ipol] = tck
                    self.u[isub][ipol]   = u
                    self.fp[isub][ipol]  = fp
                    self.ier[isub][ipol] = ier
                    self.msg[isub][ipol] = msg

                    if model_name is None:
                        self.model_name[isub][ipol] = self.source + '_%s.spl'%(states[ipol])
                    else:
                        self.model_name[isub][ipol] = model_name + '_%s.spl'%(states[ipol])

                    self.model[isub,ipol] = model
                    self.modelx[isub,ipol] = modelx

                    #self.masks are np.shape(nsub, 4, nchan, nbin) but masks are the same for each pol so it doesn't matter
                    self.model_masked[isub,ipol] = self.model[isub,ipol] * self.masks[isub,ipol]

                    if not quiet:
                        if proj_port.sum():
                            print(
                                "B-spline interpolation model %s uses %d basis profile components and %d breakpoints (%d B-splines with k=%d)." % (
                                self.model_name[isub][ipol], 
                                ncomp,
                                len(np.unique(tck[0])),
                                len(tck[0]) - tck[2] - 1, tck[2]))
                        else:
                            print("Polstate: %s"%states[ipol])
                            print(
                                "B-spline interpolation model %s uses 0 basis profile components; it returns the average profile." % (
                                    self.model_name[isub][ipol]))

        #####################################
        ## Intensity = same as original code
        #####################################
        if state=='Intensity':
            # Definitions
            port = self.portx
            pca_weights = old_div(self.SNRsxs, np.sum(self.SNRsxs))
            mean_prof = old_div((port.T * pca_weights).T.sum(axis=0), pca_weights.sum())
            freqs = self.freqsxs[0]
            nu_lo = freqs.min()
            nu_hi = freqs.max()
            # Check nbin
            nbin = port.shape[1]
            if nbin % 2 != 0:
                print("nbin = %d is odd; cannot wavelet_smooth.\n" % nbin)
                smooth = False
            elif np.modf(np.log2(nbin))[0] != 0.0:
                print(
                    "nbin = %d is not a power of two; can only try wavelet_smooth to one level; recommend resampling to a power-of-two number of phase bins.\n" % nbin)
            # Do principal component analysis
            eigval, eigvec = pca(port, mean_prof, pca_weights, quiet=quiet)
            # Get "significant" eigenvectors
            if max_ncomp is None:
                return_max = 10
            else:
                return_max = min(max_ncomp, 10)
            if smooth:
                if 'pywt' not in sys.modules:
                    raise ImportError("You failed to import pywt and need PyWavelets to use smooth=True!")
                ieig, smooth_eigvec = find_significant_eigvec(eigvec, check_max=10,
                                                            return_max=return_max, snr_cutoff=snr_cutoff,
                                                            return_smooth=True, rchi2_tol=rchi2_tol, **kwargs)
            else:
                ieig = find_significant_eigvec(eigvec, check_max=10,
                                            return_max=return_max, snr_cutoff=snr_cutoff,
                                            return_smooth=False, rchi2_tol=rchi2_tol, **kwargs)
            ncomp = len(ieig) #Why is this a float?

            if smooth:
                smooth_mean_prof = smart_smooth(mean_prof,
                                                rchi2_tol=rchi2_tol)

            if ncomp == 0:  # Will make model with constant average port
                proj_port = port[:, :ncomp]
                if smooth:
                    modelx = reconst_port = np.tile(smooth_mean_prof,
                                                    len(freqs)).reshape(len(freqs), port.shape[1])
                    model = np.tile(smooth_mean_prof,
                                    len(self.freqs[0])).reshape(len(self.freqs[0]),
                                                                port.shape[1])
                else:
                    modelx = reconst_port = np.tile(mean_prof,
                                                    len(freqs)).reshape(len(freqs), port.shape[1])
                    model = np.tile(mean_prof,
                                    len(self.freqs[0])).reshape(len(self.freqs[0]),
                                                                port.shape[1])
            else:
                delta_port = port - mean_prof
                if smooth:
                    reconst_port = reconstruct_portrait(port, mean_prof,
                                                        smooth_eigvec[:, ieig])
                    # Find the projections of the profiles onto the basis components
                    proj_port = np.dot(delta_port, smooth_eigvec[:, ieig])
                else:
                    reconst_port = reconstruct_portrait(port, mean_prof,
                                                        eigvec[:, ieig])
                    # Find the projections of the profiles onto the basis components
                    proj_port = np.dot(delta_port, eigvec[:, ieig])

            if ncomp == 0:
                (tck, u) = [np.array([]), np.array([]), 0], np.array([])
                fp, ier, msg = None, None, None
            else:

                spl_weights = pca_weights
                s = sfac * len(proj_port) * \
                        np.sum((self.SNRsxs * self.noise_stdsxs)**2) / \
                        sum(self.SNRsxs)**2
                if self.bw < 0: flip = -1   #u in si.splprep has to be increasing...
                else: flip = 1
                #Find the B-spline curve traced by the projected vectors,
                #parameterized by frequency
                (tck,u), fp, ier, msg = si.splprep(proj_port[::flip].T,
                        w=spl_weights[::flip], u=freqs[::flip], ub=nu_lo, ue=nu_hi,
                        k=k, task=0, s=s, t=None, full_output=1, nest=None, per=0,
                        quiet=int(quiet))

                if max_nbreak is not None and len(np.unique(tck[0])) > max_nbreak:
                    if max_nbreak < 2:
                        print("max_nbreak not >= 2; setting max_nbreak = 2...")
                        max_nbreak = 2
                    if max_nbreak == 2: s = np.inf
                    (tck, u), fp, ier, msg = si.splprep(proj_port[::flip].T,
                                                        w=spl_weights[::flip], u=freqs[::flip], ub=nu_lo,
                                                        ue=nu_hi, k=k, task=0, s=s, t=None, full_output=1,
                                                        nest=max_nbreak + (k * 2), per=0, quiet=int(quiet))

                if ier > 1:  # Will also catch when ier == "unknown"
                    print("Something went wrong in si.splprep for %s:\n%s" % (
                        self.source, msg))

            # Build model
            if ncomp != 0:
                if smooth:
                    modelx = gen_spline_portrait(smooth_mean_prof, freqs,
                                                smooth_eigvec[:, ieig], tck)
                    model = gen_spline_portrait(smooth_mean_prof, self.freqs[0],
                                                smooth_eigvec[:, ieig], tck)
                else:
                    modelx = gen_spline_portrait(mean_prof, freqs, eigvec[:, ieig],
                                                tck)
                    model = gen_spline_portrait(mean_prof, self.freqs[0],
                                                eigvec[:, ieig], tck)

            # Assign new attributes
            self.ieig = ieig
            self.ncomp = ncomp
            self.eigvec = eigvec
            self.eigval = eigval
            self.mean_prof = mean_prof
            if smooth:
                self.smooth_mean_prof = smooth_mean_prof
                self.smooth_eigvec = smooth_eigvec
            self.proj_port = proj_port
            self.reconst_port = reconst_port
            # tck contains the knot locations t, B-spline coefficients c, and
            # polynomial degree k -- end knots will have multiplicity k+1, interior
            # breakpoints will have multiplicity 1 for maximum continuity.  The
            # number of B-splines will be n = l + k, where l is the number of
            # intervals.  l = number of breakpoints - 1 = number of unique knots - 1
            # = len(tck[0]) - 2*tck[2] - 1.
            self.tck, self.u, self.fp, self.ier, self.msg = tck, u, fp, ier, msg

            if model_name is None:
                self.model_name = self.source + '.spl'
            else:
                self.model_name = model_name

            self.model = model
            self.modelx = modelx
            self.model_masked = self.model * self.masks[0, 0]

            if not quiet:
                if proj_port.sum():
                    print(
                        "B-spline interpolation model %s uses %d basis profile components and %d breakpoints (%d B-splines with k=%d)." % (
                        self.model_name, ncomp,
                        len(np.unique(self.tck[0])),
                        len(self.tck[0]) - self.tck[2] - 1, self.tck[2]))
                else:
                    print(
                        "B-spline interpolation model %s uses 0 basis profile components; it returns the average profile." % (
                            self.model_name))

    def write_model(self, state='Intensity', resultdir=None, quiet=False):
        """
        Write the output (pickle file) model to outfile.
        If resuldir=None, file will be saved in current directory
        """
        
        if state=='Intensity':
            outfile = self.model_name if resultdir is None else resultdir + '/' + self.model_name
            of = open(outfile, "wb")
            if hasattr(self, "smooth_eigvec"):
                if len(self.ieig):
                    
                    pickle.dump([self.model_name, self.source, self.datafile,
                                self.smooth_mean_prof, self.smooth_eigvec[:, self.ieig],
                                self.tck], of, protocol=2)
                else:
                    pickle.dump([self.model_name, self.source, self.datafile,
                                self.smooth_mean_prof, self.smooth_eigvec[:, []],
                                self.tck], of, protocol=2)
            else:
                if len(self.ieig):
                    pickle.dump([self.model_name, self.source, self.datafile,
                                self.mean_prof, self.eigvec[:, self.ieig], self.tck], of, protocol=2)
                else:
                    pickle.dump([self.model_name, self.source, self.datafile,
                                self.mean_prof, self.eigvec[:, []], self.tck], of,
                                protocol=2)

            of.close()
            if not quiet:
                print("Wrote modelfile %s." % outfile)
            
        if state=='Stokes' or state=='Stokes_L':
            nsub = self.nsub
            npol = self.npol if state=='Stokes' else 3
            for isub in range(nsub):
                for ipol in range(npol):
                    outfile = self.model_name[isub][ipol] if resultdir is None else resultdir + '/' + self.model_name[isub][ipol]
                    of = open(outfile, "wb")
                    if hasattr(self, "smooth_eigvec"):
                        if len(self.ieig[isub][ipol]):
                            ieig = self.ieig[isub][ipol] #number of significang eig for this isub,ipol
                            
                            pickle.dump([self.model_name[isub][ipol], self.source, self.datafile,
                                        self.smooth_mean_prof[isub,ipol], self.smooth_eigvec[isub,ipol][:, ieig],
                                        self.tck[isub][ipol]], of, protocol=2)
                        else:
                            pickle.dump([self.model_name[isub][ipol], self.source, self.datafile,
                                        self.smooth_mean_prof[isub,ipol], self.smooth_eigvec[isub,ipol][:, []],
                                        self.tck[isub][ipol]], of, protocol=2)
                    else:
                        if len(self.ieig):
                            ieig = self.ieig[isub][ipol] #number of significang eig for this isub,ipol
                            pickle.dump([self.model_name[isub,ipol], self.source, self.datafile,
                                        self.mean_prof[isub,ipol], self.eigvec[isub,ipol][:, self.ieig], 
                                        self.tck[isub][ipol]], of, protocol=2)
                        else:
                            pickle.dump([self.model_name[isub,ipol], self.source, self.datafile,
                                        self.mean_prof[isub,ipol], self.eigvec[:, []], 
                                        self.tck[isub][ipol]], of,
                                        protocol=2)

                    of.close()
                    if not quiet:
                        print("Wrote modelfile %s." % outfile)

    def show_eigenprofiles(self, state='Intensity', ncomp=None, resultdir='', savefig=False, **kwargs):
        """
        Calls show_eigenprofiles(...) to make plots of mean/eigen profiles.

        see show_eigenprofiles(...) for details.

        ncomp=None plots self.ncomp PCA components, otherwise plots the number
            of components specified.
        resultdir on default results are saved in current directory unless specified
        **kwargs get passed to show_eigenprofiles(...).
        """
        if state=='Intensity':

            title = self.model_name[:-4]

            if ncomp is None: ncomp = self.ncomp
            if hasattr(self, "smooth_eigvec"):
                if ncomp:
                    eigvec = self.eigvec[:, self.ieig[:ncomp]].T          #np.shape(ncomp, nbin)
                    seigvec = self.smooth_eigvec[:, self.ieig[:ncomp]].T  #np.shape(ncomp, nbin)
                else:
                    eigvec = None
                    seigvec = None
                show_eigenprofiles(eigvec, seigvec, self.mean_prof,
                                self.smooth_mean_prof, title=title, resultdir=resultdir, savefig=savefig, 
                                **kwargs)
            else:
                if ncomp:
                    eigvec = self.eigvec[:, self.ieig[:ncomp]].T
                else:
                    eigvec = None
                show_eigenprofiles(eigvec, None, self.mean_prof, None, title=title,
                                   resultdir=resultdir, savefig=savefig,
                                **kwargs)
                
        if state=='Stokes' or state=="Stokes_L":
            nsub = self.nsub 
            npol = self.npol if state=='Stokes' else 3
            for isub in range(nsub):
                for ipol in range(npol):
                    title = self.model_name[isub][ipol][:-4]

                    if ncomp is None: ncompp = int(self.ncomp[isub,ipol])

                    eigvec  = self.eigvec[isub,ipol]
                    ieig    = self.ieig[isub][ipol]
                
                    if hasattr(self, "smooth_eigvec"):
                        seigvec = self.smooth_eigvec[isub,ipol]
                        if ncompp:
                            eigvec  = eigvec[:, ieig[:ncompp]].T
                            seigvec = seigvec[:, ieig[:ncompp]].T
                        else:
                            eigvec = None
                            seigvec = None

                        show_eigenprofiles(eigvec, seigvec, self.mean_prof[isub,ipol],
                                        self.smooth_mean_prof[isub,ipol], title=title,
                                        resultdir=resultdir, savefig=savefig,
                                        **kwargs)
                    else:
                        if ncomp:
                            eigvec = eigvec[:, ieig[:ncomp]].T
                        else:
                            eigvec = None

                        show_eigenprofiles(eigvec, None, self.mean_prof[isub,ipol], None, title=title,
                                           resultdir=resultdir, savefig=savefig,
                                        **kwargs)

    def show_spline_curve_projections(self, state='Intensity', ncomp=None, resultdir='', savefig=False, **kwargs):
        """
        Calls show_spline_curve_projections(...) to make plots of the model.

        see show_spline_curve_projections(...) for details.

        ncomp=None plots self.ncomp PCA components, otherwise plots the number
            of components specified.

        xaxis: can be freqs (default) or phase bin numbers (useful for cleaning data)
        **kwargs get passed to show_spline_curve_projections(...).
        """
        if state=='Intensity':
            title = self.model_name[:-4]
            if ncomp is None: ncomp = self.ncomp
            if ncomp:
                show_spline_curve_projections(self.proj_port, self.tck,
                                            self.freqsxs[0], old_div(self.SNRsxs, np.sum(self.SNRsxs)),
                                            ncoord=ncomp, title=title, resultdir=resultdir, savefig=savefig, **kwargs)
        if state=='Stokes' or state=="Stokes_L":
            nsub = self.nsub
            if state=='Stokes':
                npol = self.npol
            if state=='Stokes_L':
                npol = 3
            for isub in range(nsub):
                for ipol in range(npol):
                    title = self.model_name[isub][ipol][:-4]
                    ncompp = int(self.ncomp[isub,ipol])
                    show_spline_curve_projections(self.proj_port[isub][ipol], self.tck[isub][ipol],
                                                self.freqsxs[0], old_div(self.SNRsxs[isub][ipol], np.sum(self.SNRsxs[isub,ipol])),
                                                ncoord=ncompp, title=title, resultdir=resultdir, savefig=savefig, **kwargs)
        
    def get_outliers(self, state='Intensity', ncomp=None, resultdir='.', savefig=False, stokes=0,**kwargs):
    
        if state=='Stokes' or state=="Stokes_L":
            nsub = self.nsub
            if state=='Stokes':
                npol = self.npol
            if state=='Stokes_L':
                npol = 3
            for isub in range(nsub):
                for ipol in range(npol):
                    ncompp = int(self.ncomp[isub,ipol])
                    get_outliers(self.proj_port[isub][ipol], self.tck[isub][ipol],
                                                self.freqsxs[0], old_div(self.SNRsxs[isub][ipol], np.sum(self.SNRsxs[isub,ipol])),
                                                ncoord=ncompp, stokes=stokes)

if __name__ == "__main__":

    from optparse import OptionParser

    usage = "Usage: %prog -d <datafile> [options]"
    parser = OptionParser(usage)
    # parser.add_option("-h", "--help",
    #                  action="store_true", dest="help", default=False,
    #                  help="Show this help message and exit.")
    parser.add_option("-d", "--datafile",
                      action="store", metavar="archive", dest="datafile",
                      help="PSRCHIVE archive from which to make model, or a metafile listing multiple archives (i.e., from different bands).  If providing a metafile, the achives must already be aligned.")
    parser.add_option("-o", "--resultdir",
                      action="store", metavar="resultdir", dest="resultdir",
                      help="Specify where to save the results. Default current directory.")
    parser.add_option("-l", "--model_name",
                      action="store", metavar="model_name", dest="model_name",
                      default=None,
                      help="Optional name for model [default=datafile.spl]. If state='Stokes' extension is fixed and its a list: [model_name_npol.spl]")
    parser.add_option("-a", "--archive",
                      action="store", metavar="archive", dest="archive",
                      default=False,
                      help="If True, will save the spline model to as an archive 2D template.")
    parser.add_option("-N", "--norm",
                      action="store", metavar="normalization", dest="norm",
                      default="prof",
                      help="Normalize the input data by channel ('None', 'mean', 'max' (not recommended), 'rms' (off-pulse noise), 'prof' (mean profile flux) [default], or 'abs' (sqrt{vector modulus})).")
    parser.add_option("-st", "--state",
                      action="store", metavar="polstate", dest="state",
                      default="Intensity",
                      help="Choose either 'Intensity' or 'Stokes'." )
    parser.add_option("-s", "--smooth",
                      action="store_true", metavar="smooth", dest="smooth",
                      default=False,
                      help="Smooth the eigenvectors and mean profile [recommended] using default wavelet_smooth options and smart_smooth.")
    parser.add_option("-n", "--max_ncomp",
                      action="store", metavar="max_ncomp", dest="max_ncomp",
                      default=10,
                      help="Maximum number of principal components to use in PCA reconstruction of the data.  max_ncomp is limited to a maximum of 10 by the B-spline representation in scipy.interpolate.")
    parser.add_option("-S", "--snr",
                      action="store", metavar="snr_cutoff", dest="snr_cutoff",
                      default=150.0,
                      help="S/N ratio cutoff for determining 'significant' eigenprofiles.  A value somewhere over 100.0 should be good. [default=150.0].")
    parser.add_option("-T", "--rchi2_tol",
                      action="store", metavar="tolerance", dest="rchi2_tol",
                      default=0.1,
                      help="Tweak this between 0.0 and 0.1 [default] if the returned eigenprofiles are not smooth enough.")
    parser.add_option("-k", "--degree",
                      action="store", metavar="degree", dest="k", default=3,
                      help="Degree of the spline.  Cubic splines (k=3) are recommended [default]. 1 <= k <=5.")
    parser.add_option("-f", "--sfac",
                      action="store", metavar="smooth_factor", dest="sfac",
                      default=1.0,
                      help="To change the smoothness of the B-spline model, tweak this between 0.0 (interpolating spline that passes through all data points) and a large number (guarantees maximum two breakpoints = maximum smoothness).  Alternatively, use -t.")
    parser.add_option("-t", "--knots",
                      action="store", metavar="max_knots", dest="max_nbreak",
                      default=None,
                      help="The maximum number of unique knots.  This functions esentially as an ignorant smoothing condition in case the default settings return a fit with more than max_knots number of unique knots in the spline model.  e.g., 10 unique knots are more than usually necessary.")
    parser.add_option("--plots",
                      action="store_true", dest="make_plots", default=False,
                      help="Save some plots related to the model with basename model_name (-l).")
    parser.add_option("--quiet",
                      action="store_true", dest="quiet", default=False,
                      help="Suppresses output.")

    (options, args) = parser.parse_args()

    if (options.datafile is None):
        print("\nppspline.py - make a pulse portrait model using PCA & B-spline interpolation\n")
        parser.print_help()
        print("")
        parser.exit()

    datafile = options.datafile
    resultdir = options.resultdir
    model_name = options.model_name
    archive = options.archive
    norm = options.norm
    state = options.state
    smooth = options.smooth
    max_ncomp = int(options.max_ncomp)
    snr_cutoff = float(options.snr_cutoff)
    rchi2_tol = float(options.rchi2_tol)
    k = int(options.k)
    sfac = float(options.sfac)
    if options.max_nbreak is not None:
        max_nbreak = int(options.max_nbreak)
    else:
        max_nbreak = None
    make_plots = options.make_plots
    quiet = options.quiet

    dp = DataPortrait(datafile, state=state, quiet=quiet)

    if norm in ("mean", "max", "prof", "rms", "abs"):
        #Check if the portrait was already normalized
        dp.normalize_portrait(norm)

    dp.make_spline_model(state=state, max_ncomp=max_ncomp, smooth=smooth,
                         snr_cutoff=snr_cutoff, rchi2_tol=rchi2_tol, k=k, sfac=sfac,
                         max_nbreak=max_nbreak, model_name=model_name, quiet=quiet)

    #If Stokes, there will be 4 models saved for each state
    dp.write_model(state=state, resultdir=resultdir, quiet=quiet)

    if archive is not None and len(dp.datafiles) == 1:
        #Even if Stokes, it will be save in only one file
        dp.write_model_archive(state=state, resultdir=resultdir, quiet=quiet)

    if make_plots:
        dp.show_eigenprofiles(state=state, resultdir=resultdir, savefig=True)
        dp.show_spline_curve_projections(state=state, resultdir=resultdir, savefig=True)
        dp.show_model_fit(state=state, resultdir=resultdir, savefig=True)
