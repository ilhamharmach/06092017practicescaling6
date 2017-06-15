#! /usr/bin/env python

import numpy as NP
import pylab as PL
import FlowCytometryTools
from FlowCytometryTools import FCMeasurement
import glob as glob
import pickle
from scipy import stats

# FUNCTION: Loads data and converts to numpy array

# datafile=r'Data.001'

def load_file(File_name):
    Path = File_name
    sample = FCMeasurement(ID='Test Sample', datafile=Path)
    FSC=NP.array(sample.data[['FSC-H']]) # Forward scatter
    SSC=NP.array(sample.data[['SSC-H']]) # Size scatter
    GFP=NP.array(sample.data[['FL1-H']]) #GFP
    current_data=NP.array(sample.data[['FSC-H','SSC-H','FL1-H']]) # matrix that contains all of them
    sample_id = sample.meta[u'SAMPLE ID']
    return current_data,sample_id



# This function takes our current data file and erases some values that are above a threshold we fixed in
# order to get rid of all the supplementary signals the flow  cytometer takes into account
# and that we do not need

def gate_cells(cells,fsc_thr,ssc_thr,cut_max = True, cut_GFP_max = True):
    ind = NP.where(cells[:,0]>fsc_thr)
    cells2 = cells[ind[0],:]
    ind2 = NP.where(cells2[:,1] > ssc_thr)
    cells3 = cells2[ind2[0],:]
    if cut_max == True:
        ind3 = NP.where(cells3[:,1] < NP.max(cells3[:,1]))
        cells4 = cells3[ind3[0],:]
    else:
        cells4 = cells3
    if cut_GFP_max == True:
        ind4 = NP.where(cells4[:,2] < NP.max(cells4[:,2]))
        cells5 = cells4[ind4[0],:]
    else:
        cells5 = cells4
    return cells5




# Writing a function that will directly return all the sample IDs, the raw data and the filtered data for our files (so we won't be doing it with each file)

def filter_all_data(fsc_thr, ssc_thr, max_treat = True):
    all_data = {}
    strain_list = []
    file_list = glob.glob('Data*')
    for file in file_list:
        raw_data, name = load_file(file)
        strain_list.append(name)
    
    strain_list = NP.unique(NP.array(strain_list)) # Why?
    
    for file in file_list:
        raw_data, name = load_file(file)
        filtered_data = gate_cells(raw_data, fsc_thr, ssc_thr, cut_max = max_treat)
        all_data[name] = filtered_data  # putting it in a dictionary, and the data is callable from the name of the gene
    return all_data, strain_list # We return a dictionary that contains all the data for every file and a list that contains all the gene names that are present in the files.



# FUNCTION: Bin data by size
# Takes all data and bins across range of values.  Constant size bins enforced to permit
# background subtraction bin-wise.
    
# It bins the data meaning for x it takes the middle and for y it takes the mean , so for a range of points we only select one point , and it gives us a better understanding of the data , and a clearer vieuw of the curve.

def plot_size_bins(x, y, n_bin=100, plot_range = [0, 1023]):
	if plot_range == None:
		bins = NP.linspace(min(x),max(x),n_bin)
	else:
		bins = NP.linspace(plot_range[0],plot_range[1],n_bin)
	plot_data = NP.vstack((x,y)).transpose()
	binned_data = []
	for bin in range(len(bins)-1):
		temp_bin = plot_data[NP.where(plot_data[:,0] >= bins[bin])[0],:]
		temp_bin2 = temp_bin[NP.where(temp_bin[:,0] < bins[bin + 1]),:]
		if len(temp_bin2[0][:,0]) > 0:
			mean = NP.mean(temp_bin2[0][:,1])
			sd = NP.std(temp_bin2[0][:,1])
			binned_data.append([NP.mean([bins[bin], bins[bin + 1]]),mean, sd])
	return NP.array(binned_data)



# This function bins the data for the SSC and the GFP and puts the binned data in a new dictionary , so now we have the filtered data that we want

def bin_all(data):
    binned_data = {}
    for strain in data.keys():
        x = data[strain][:,1]
        y = data[strain][:,2]
        temp = plot_size_bins(x,y)
        binned_data[strain] = temp
    return binned_data






# These functions's goal is to substract the unnecessary things from our plot to be able to plot the filtered data
# For us ,the bg_strain is 'BY'
# Definition of the BG_strain?
# The strain list contains the names of the genes that we have in our files


def bg_sub_wrap(raw_data, binned_data, bg_strain, strain_list):
    if type(bg_strain) == str:
        bg_strain = [bg_strain] * len(strain_list) # why exactly
        if type(bg_strain) == list:
            bg_strain = bg_strain
        new_raw_data = {}
        new_binned_data = {}
        for sample in range(len(strain_list)):
            new_temp_binned, new_temp_raw = bg_sub(raw_data, binned_data, bg_strain, sample, strain_list)
            new_binned_data[strain_list[sample]] = new_temp_binned
            new_raw_data[strain_list[sample]] = new_temp_raw
    return new_binned_data, new_raw_data



#?

def bg_sub(raw_data, binned_data, bg_strain, sample, strain_list):
    strain = strain_list[sample] # One gene
    print strain
    bg = bg_strain[sample]
    print bg
    bg_bin_fit = bg_fit = NP.polyfit(raw_data[bg][:,1],raw_data[bg][:,2],1) #?
    bg_function = NP.poly1d(bg_fit) #?
    bg_sub_val = binned_data[strain][:,1] - bg_function(binned_data[strain][:,0])
    all_vals = NP.vstack([NP.transpose(binned_data[strain]), bg_sub_val]).transpose()
    raw_sub_val = NP.zeros(len(raw_data[strain][:,1]))
    bg_bin_fit = bg_fit = NP.polyfit(raw_data[bg][:,1],raw_data[bg][:,2],1)
    bg_function = NP.poly1d(bg_fit)
    sample_data = raw_data[strain]
    corrected_GFP = sample_data[:,2] - bg_function(sample_data[:,1])
    corrected_data = \
    NP.vstack((sample_data.transpose(),corrected_GFP)).transpose()
    all_vals_raw = corrected_data
    return all_vals, all_vals_raw



# Plot background subtracted graphs.
# Takes all strains in current data from binned data and plots them individually.  Saves
# plots as PDF.

def plot_save_data(proc_data, binned_data, file_type = 'pdf'):# What is proc_data?
    my_final_dictionary={}
    for sample in proc_data.keys():
        temp_data = proc_data[sample]
        temp_bin = binned_data[sample]
        PL.plot(proc_data[sample][:,1], proc_data[sample][:,3], '.', markersize = 12, \
        alpha = .05, label = "{0}".format(sample))
        PL.plot(temp_bin[:,0],temp_bin[:,3],'-k', linewidth = 2)
        a=proc_data[sample][:,1]
        b=proc_data[sample][:,3]
        slope, intercept =linear_model(a,b)
        my_final_dictionary[sample]={'proc_data':proc_data[sample],'binned data':temp_bin,'intercept':intercept,'scaling factor':intercept/NP.mean(proc_data[sample][:,3])}
        PL.plot(a, intercept+slope*a, 'r')
        PL.legend(loc = "upper left")
        PL.xlabel("SSC (Volume)")
        PL.ylabel("GFP Intensity (AU)")
        PL.ylim(0,1050)
        PL.xlim(0,1050)
        savepath = '{0}.{1}'.format(sample,file_type)
        PL.savefig(savepath)
        PL.close()
        pickle.dump(my_final_dictionary, open("all_data.p", "w"))


# Where to create the dictionary?

def analyze_data(fsc_thr, ssc_thr, bg_strain):
    print "Analyzing current folder..."
    print "Gating data..."
    data, strains = filter_all_data(fsc_thr, ssc_thr)
    print "Current folder data: {0}".format(strains)
    print "Binning data..."
    bin_data = bin_all(data)
    print "Performing background subtraction..."
    new_bin, new_raw = bg_sub_wrap(data, bin_data, bg_strain, strains)
    print "Plotting all data..."
    plot_save_data(new_raw, new_bin, file_type = "tif")
    print "Saving all data..."
    pickle.dump(new_raw, open("raw_data.p", "w"))
    pickle.dump(new_bin, open("binned_data.p", "w"))
    print "Process complete."






def linear_model(c,d):
    slope, intercept, r_value, p_value, std_err = stats.linregress(c,d)
    return slope, intercept





