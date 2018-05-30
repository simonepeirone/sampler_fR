import os, math, sys

here = os.path.dirname(os.path.abspath(__file__))
getdist_python_path = here+'/../python/'
sys.path.insert(0, os.path.normpath(getdist_python_path))

import getdist.plots as gplot
import getdist.chains as chains
import getdist.mcsamples
import getdist.densities as dens
import planckStyle
import numpy as np

cache = False

results_dir = here
chains_dir = here+'/../chains'

stat_results_dir = results_dir+'/stats/'

analysis_settings = {'max_corr_2D': u'0.99',
		     'boundary_correction_order': u'0',
		     'converge_test_limit': u'0.95',
		     'smooth_scale_2D': u'0.5',
		     'credible_interval_threshold': u'0.05',
		     'contours': u'0.68 0.95 0.997',
		     'fine_bins_2D': u'256',
		     'num_bins': u'200',
		     'mult_bias_correction_order': u'0',
		     'fine_bins': u'1024',
		     'num_bins_2D': u'80',
		     'max_scatter_points': u'2000',
		     'range_ND_contour': u'0',
		     'range_confidence': u'0.001',
		     'smooth_scale_1D': u'0.3',
		     'ignore_rows': u'0.3'}

chains_roots = [
				'test_ph'
                ]


# create the results folder if needed:
if not os.path.isdir(results_dir):
    os.mkdir(results_dir)

if not os.path.isdir(stat_results_dir):
    os.mkdir(stat_results_dir)

# get the marginal constraints without post-processing:
for root in chains_roots:
    # get the samples:
    sample = getdist.loadMCSamples(chains_dir+'/'+root, settings=analysis_settings, no_cache=cache)
    # digest input:
    sample.out_dir = stat_results_dir
    sample.rootdirname = os.path.join(stat_results_dir, root)

    # write covariance:
    sample.writeCovMatrix()
    # write correlation:
    sample.writeCorrelationMatrix()
    # write marginal statistics:
    sample.getMargeStats().saveAsText( sample.rootdirname + '.margestats' )

# post process and plot the relevant chains:
sample_1 = getdist.loadMCSamples(chains_dir+'/test_ph', settings=analysis_settings, no_cache=cache)
# sample_2 = getdist.loadMCSamples(chains_dir+'/new_LCDM_kids', settings=analysis_settings, no_cache=cache)

# # test with the triangle before resampling:
# g1=gplot.getSubplotPlotter()
# params = [u'omegabh2',u'omegach2',u'theta',u'tau',u'logA',u'ns',u'sigma8']
# label=['LCDM', 'CPL']
# g1.add_legend(label,fontsize=10)
# g1.triangle_plot( [sample_1, sample_2] ,params, filled=True, shaded=False)
# g1.export(results_dir+'1_test_pre_kids.pdf')
#
# # test with the triangle before resampling:
# g1=gplot.getSubplotPlotter()
# params = [u'omegabh2',u'omegach2',u'theta',u'tau',u'logA',u'ns',u'sigma8']
# label=['LCDM', 'CPL']
# g1.add_legend(label,fontsize=10)
# g1.triangle_plot( [sample_planck_LCDM, sample_planck_CPL] ,params, filled=True, shaded=False)
# g1.export(results_dir+'1_test_pre_planck.pdf')
#
# # test with the triangle before resampling:
# g1=gplot.getSubplotPlotter()
# params = [u'omegabh2',u'omegach2',u'theta',u'tau',u'logA',u'ns',u'sigma8']
# label=['LCDM', 'CPL']
# g1.add_legend(label,fontsize=10)
# g1.triangle_plot( [sample_kidsplanck_LCDM, sample_kidsplanck_CPL] ,params, filled=True, shaded=False)
# g1.export(results_dir+'1_test_pre_kidsplanck.pdf')


# get the density of the first chain:
params_list = ['omegabh2', 'omegach2']
params_list = ['tau', 'logA']
params_list = ['omegabh2', 'omegach2', 'theta', 'tau', 'logA', 'ns']
par_numbers = [ sample_1.getParamNames().numberOfName(par) for par in params_list]

mean        = sample_1.getMeans()[par_numbers]
covariance  = sample_1.cov()[:,par_numbers][par_numbers,:]

inv_covariance = np.linalg.inv( covariance )

def gaussian( mean, x, inv_cov):
	return 0.1*np.dot(np.dot((mean-x),inv_cov),(mean-x))

par_numbers_1 = [ sample_1.getParamNames().numberOfName(par) for par in params_list]

new_weight = []
for weight, point in zip( sample_1.weights, sample_1.samples):
	x = point[par_numbers_1]
	new_weight.append( gaussian( mean, x, inv_covariance) )

new_weight = np.array( new_weight )
sample_1.reweightAddingLogLikes(new_weight)
sample_1.updateBaseStatistics()

# # get the density of the second chain:
# params_list = ['omegabh2', 'omegach2']
# params_list = ['tau', 'logA']
# params_list = ['omegabh2', 'omegach2', 'theta', 'tau', 'logA', 'ns']
# par_numbers = [ sample_2.getParamNames().numberOfName(par) for par in params_list]
#
# mean        = sample_2.getMeans()[par_numbers]
# covariance  = sample_2.cov()[:,par_numbers][par_numbers,:]
#
# inv_covariance = np.linalg.inv( covariance )
#
# def gaussian( mean, x, inv_cov):
# 	return 0.1*np.dot(np.dot((mean-x),inv_cov),(mean-x))
#
# par_numbers_2 = [ sample_2.getParamNames().numberOfName(par) for par in params_list]
#
# new_weight = []
# for weight, point in zip( sample_2.weights, sample_2.samples):
# 	x = point[par_numbers_2]
# 	new_weight.append( gaussian( mean, x, inv_covariance) )
#
# new_weight = np.array( new_weight )
# sample_2.reweightAddingLogLikes(new_weight)
# sample_2.updateBaseStatistics()

# # test with the triangle after resampling:
# g2=gplot.getSubplotPlotter( )
# params = [u'omegabh2',u'omegach2',u'theta',u'tau',u'logA',u'ns',u'sigma8']
# label=['LCDM', 'CPL']
# g2.add_legend(label,fontsize=10)
# g2.triangle_plot( [sample_1, sample_2] ,params, filled=True, shaded=False)
# g2.export(results_dir+'2_test_post_kids.pdf')
#
# # test with the triangle after resampling:
# g2=gplot.getSubplotPlotter( )
# params = [u'omegabh2',u'omegach2',u'theta',u'tau',u'logA',u'ns',u'sigma8']
# label=['LCDM', 'CPL']
# g2.add_legend(label,fontsize=10)
# g2.triangle_plot( [sample_planck_LCDM, sample_planck_CPL] ,params, filled=True, shaded=False)
# g2.export(results_dir+'2_test_post_planck.pdf')
#
# # test with the triangle after resampling:
# g2=gplot.getSubplotPlotter( )
# params = [u'omegabh2',u'omegach2',u'theta',u'tau',u'logA',u'ns',u'sigma8']
# label=['LCDM', 'CPL']
# g2.add_legend(label,fontsize=10)
# g2.triangle_plot( [sample_kidsplanck_LCDM, sample_kidsplanck_CPL] ,params, filled=True, shaded=False)
# g2.export(results_dir+'2_test_post_kidsplanck.pdf')

# get the bounds:
sample_1.rootdirname = os.path.join(results_dir, 'test_ph')
# sample_2.rootdirname = os.path.join(results_dir,'new_LCDM_kids')

all_samples=[sample_1]#, sample_2]
for sample in all_samples:
	sample.out_dir = results_dir
	# write covariance:
	sample.writeCovMatrix()
	# write correlation:
	sample.writeCorrelationMatrix()
	# write marginal statistics:
	sample.getMargeStats().saveAsText( sample.rootdirname + '.margestats' )
	# write likelihood statistics:
	sample.getLikeStats().saveAsText( sample.rootdirname + '.likestat'    )

#Plot only kids+Planck
g2 = gplot.getSinglePlotter(width_inch=4, ratio=0.7)
g2.plot_2d( sample_1, 'EFTw0', 'EFTwa', filled=True)#, lims=[  0.0,0.5,0.5,1.5])
# g2.add_legend(['old run','new run'], fontsize=10, legend_loc='best', colored_text=True)
g2.export(results_dir+'/test_ph.pdf')


exit()
