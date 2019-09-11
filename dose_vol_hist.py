# function calculating cumulative frequencies
from scipy.stats import cumfreq
# for creating excel files
import xlsxwriter

### Producing data for a single dose-volume histogram
# Given 
#   d - a list (or 1-D NumPy array) of doses in voxels of a ROI,
#   numbins - number of bins dividing the range of doses.
# Returns a set of 2D-points defining a dose-volume histogram
#   * Fractions of volumes are dimensionless, i.e. from 0 to 1  
def dvh_calc(d, numbins=100):
    nvox=len(d) # number of voxels in the ROI
    ## Data of cumulative frequencies - CumfreqResult object
    freq_data=cumfreq(d,numbins=numbins)
    ## The lower bound of the 1st bin
    bin0=freq_data.lowerlimit
    ## Initializing the ouput array of 2-D points:
    #     if doses start above 0, add the 0 point
    out = [[0,1]] if bin0 > 0 else []
    out.append([bin0,1]) # add the point for the lower bound of the 1st bin
    ## Generating the rest of 2-D points
    out.extend(
        [ # n is the number of voxels with dose <= upper bound of current bin
          [
           # x-axis: dose
           bin0 + (i+1)*freq_data.binsize,
           # y-axis: fraction of voxels exceeding this dose
           (nvox-n)/nvox       
                  ]      
             for i,n in enumerate(freq_data.cumcount)
             ]
            )
    return out


### Producing Excel file with dose-volume histograms for one patient,
#       grouped by IMRT plans and by ROIs
# Given
#   fname - Excel file name
#   data - 3-D array describing dose distributions in several IMRT plans:
#       list[ for each plan: 
#           list[ for each ROI: 
#                   list or 1-D NumPy array of doses in the ROI ] ]
#        * the first 2 dimensions are rectangular [nr_of_plans x nr_of_ROIs],
#          the 3-rd dimension can be variable (different nrs. of voxels in ROIs)
#   plan_names - list of plans' display names
#   roi_names - list of display names for ROIs,
#   plan_abbrevs (optional) - list of strings abbreviating plans for naming worksheets
#       if None, plan_names are used instead,
#   numbins = 100 - number of bins dividing the range of doses,
#   roicharts_positions (optional) - if a sheet with ROI-specific charts is created,
#       list of positions of these charts as [cell row number, cell column number]
# Returns nothing
# Creates an Excel file:
#   * one worksheet per each IMRT plan, containing data for dose-volume histograms 
#           and the histograms plotted on one diagram
#   * if there is more than one IMRT plan, a worksheet named "Comparison",
#           containing for each volume its dose-volume histograms from all plans

def dvh2excel(fname, data, 
              plan_names, roi_names, plan_abbrevs=None, 
              numbins=100, roicharts_positions=None):
    if plan_abbrevs is None:
       plan_abbrevs=plan_names
    ## initializing an excel workbook
    workbook=xlsxwriter.Workbook(fname)
    # in case of multiple plans, a sheet with charts comparing plans for each ROI
    if len(data)>0:
        comparesheet=workbook.add_worksheet("Comparison")
        roicharts=[] # ROI-specific charts
        for i, roi_n in enumerate(roi_names):
            roicharts.append(
                workbook.add_chart({'type': 'scatter',"subtype":"smooth"}))
            roicharts[-1].set_title({"name":roi_n})
            roicharts[-1].set_x_axis({'major_gridlines': {'visible': True}})
    ## Filling the workbook with histograms
    # for each IMRT plan...
    for plan_abbr, plan_nm, plan_data in \
                zip(plan_abbrevs, plan_names, data): 
        # creating a worksheet and a diagram with histograms for this plan
        plansheet=workbook.add_worksheet(plan_abbr)
        planchart = workbook.add_chart({'type': 'scatter',"subtype":"smooth"})
        planchart.set_x_axis({'major_gridlines': {'visible': True}})
        # for each ROI...
        for i, roi_nm, roi_data in \
                zip(range(len(roi_names)), roi_names, plan_data):
            # data for the dose-volume histogram
            dvh_data=dvh_calc(roi_data, numbins=numbins)
            # name of the 2-column table for histogram data
            plansheet.write(0,i*2,roi_nm)
            # writing data of histogram points
            for j, point in enumerate(dvh_data):
                plansheet.write(j+1,i*2,point[0])
                plansheet.write(j+1,i*2+1,point[1])
            # adding series of the dvh to the plan-specific chart
            planchart.add_series({
            'categories': [plan_abbr, # worksheet name as part of excel reference
                           1, i*2, len(dvh_data), i*2], # reference to cell range
            'values':     [plan_abbr,
                           1, i*2+1, len(dvh_data), i*2+1],
            'name':       [plan_abbr, 0, i*2],
            'line': {'width':1.25}
                })
            # adding series of the dvh to the ROI-specific chart
            if len(data)>0:
                roicharts[i].add_series({
                'categories': [plan_abbr,
                               1, i*2, len(dvh_data), i*2],
                'values':     [plan_abbr,
                               1, i*2+1, len(dvh_data), i*2+1],
                'name':       plan_nm,
                'line': {'width':1.25}
                    })
        # placing the chart with plan-specific dose-volume histograms
        plansheet.insert_chart("B2",planchart)
    # placing ROI-specific charts on the "Comparison" sheet, 3 charts in a row
    if len(data)>0:
        # if not given, creating positions for charts as 3 charts in a row
        if roicharts_positions is None:
            roicharts_positions=[
                    [i*15,j*8] for i in range(int(len(roicharts)/3)+1)
                        for j in range(3)]
        for i,roi_ch in enumerate(roicharts):
            comparesheet.insert_chart(
                    roicharts_positions[i][0],
                    roicharts_positions[i][1],
                    roi_ch)
    workbook.close()