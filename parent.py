import os

import numpy as np
from scipy import sparse
from scipy.stats import cumfreq

import pandas as pd

import re # string regular expressions

import xlsxwriter
import dose_vol_hist as dvh

## List of bits of an integer starting from the least significant
def bitlist(n):
    return [1 if digit=='1' else 0 for digit in reversed(bin(n)[2:])]

## Read the m-file defining information about imrt data
def read_mfile(f):
    dic={} # init dictionary collecting result
 
    # name of the dataset
    dic["fname"]=f.readline()[:-1]
    # nr. of beams
    l=f.readline().strip("\n").split(maxsplit=1)
    if l[1]!="// liczba wiazek":
        print("Reading m-file ",f,
              ":\n   expected # // liczba wiazek,\n"
              "   got ",l)
        return
    dic["nbeams"]=int(l[0]) 
    # list of nrs of beamlets in beams
    dic["nbeamlist"]=[0 for i in range(dic["nbeams"])] 
    ## reading nrs. of beamlets for each beam
    while True:
        l=f.readline().strip("\n").split(maxsplit=1)
        try:
            dic["nbeamlist"][int(l[0])-1]=int(l[1])
        except:
            break # l contains next line-to-words after nrs. of beamlets
    dic["nvars"]=sum(dic["nbeamlist"])
    # reading nr. of voxels
    if l[1]!="// liczba vokseli":
        print("Reading m-file ",f,
              ":\n   expected # // liczba vokseli,\n"
              "   got ",l)
        return
    dic["nvox"]=int(l[0])
    # reading grid scaling
    l=f.readline().strip("\n").split(maxsplit=1)
    if l[1]!="// DoseGridScaling":
        print("Reading m-file ",f,
              ":\n   expected # // DoseGridScaling,\n"
              "   got ",l)
        return
    dic["DGScaling"]=float(l[0])
    # reading nr. of organs
    l=f.readline().strip("\n").split(maxsplit=1)
    if l[1]!="// liczba ROI":
        print("Reading m-file ",f,
              ":\n   expected # // liczba ROI,\n"
              "   got ",l)
        return
    dic["nroi"]=int(l[0])
    # reading organs list
    dic["roinames"]=["" for i in range(dic["nroi"])]
    dic["roicodes"]=[0 for i in range(dic["nroi"])]
    for i in range(dic["nroi"]):
        l=f.readline().strip("\n").split(maxsplit=1)
        dic["roicodes"][i]=len(bitlist(int(l[0])))-1
        dic["roinames"][i]=l[1]
    return dic

## Given a list of [int: int>0 <=> bit=1], return an int coded by these bits
def bits2int(l):
    return sum([2**i for i,b in enumerate(l) if b>0])

### Class patient describing an irradiation model
#   for a given patient and given equipment settings
#  # beamlets    
#  .nbeams - nr. of beams    
#  .nbeamlist - list of nrs of beamlets in all beams
#  .beams_offset - list of: for each beam, 1st beamlet index in the common list
#  # ROI    
#  .nroi - nr. of ROIs
#  .roicodes - for each ROI, nr. of the bit when coding (vox \in ROI) -> int
#  .roivox - for each ROI, list of voxel ids belonging to it        
#  .voxcodes - for each voxel, the code to which ROI it belongs 
#   (1=>1, 2=>2, 4=>3, 8=>4,...)
#       coding the ROI
class Patient:
    def __init__(self,file,dirnm=None):
    ## if file = handler, use dirnm to pass folder name
        if isinstance(file,str):
            f=open(file,"r")
            qclose=True # we will need to close the file
            self.dirnm=os.path.dirname(file)
        else:
            f=file
            qclose=False
            self.dirnm=dirnm
        if not(self.dirnm[-1]=="\\"):
            self.dirnm+="\\"        
    ## getting info from m-file
        for key, value in read_mfile(f).items():
            setattr(self,key,value)
        if qclose:
            f.close()
    ## calculating beams offset
        self.beams_offset=[sum(self.nbeamlist[0:i]) for i in range(self.nbeams)]
    ## read the d-matrices
        # column nr. to add to element indices when reading next parts
        icol=0 
        for i in range(self.nbeams):
            if i==0: # read and put the first part into the main matrix
                mm=pd.read_csv(
                    self.dirnm+"d_"+self.fname+"_"+str(i+1)+".txt",
                    sep=" ",skiprows=1,header=None).to_numpy()
            else: # read, adjust the part, then append to the main matrix
                mi=pd.read_csv(
                        self.dirnm+"d_"+self.fname+"_"+str(i+1)+".txt",
                        sep=" ",skiprows=1,header=None).to_numpy()
                mi[:,1]+=icol
                mm=np.vstack((mm,mi))
            icol+=self.nbeamlist[i]
        mm=mm.T
        self.d=sparse.csc_matrix(
            (mm[2]*self.DGScaling,(mm[0],mm[1])),
            shape=[self.nvox,icol]
            )
    ## read voxel-to-organ assignment
        vv=[int(a) for a in
            np.genfromtxt(self.dirnm+"v_"+self.fname+".txt").T]
        self.voxcodes=vv
        if len(vv)!=self.nvox:
            print("Voxel assignment length ",len(vv),
                  " vs. ",self.nvox," required!!!")
        self.roivox=[[] for i in range(self.nroi)]
        # codes of voxels -> voxel lists for ROIs 
        for i,v in enumerate(vv):
            for j,b in enumerate(bitlist(v)):
                if b==1:
                    self.roivox[j].append(i)
    ## reads a solution from txt results of Gurobi, returns a vector
    def read_gur_sol(self,file):
        # checking if fileis string 
        if isinstance(file, str):
            f=open(file)
            qclose=True # should close later
        else:
            f=file
            qclose=False
        # relative positions of x sub-vectors for all beams  
        addcols=[0]+\
            [sum(self.nbeamlist[:i]) for i in range(1,self.nbeams)]
        # initialize the solution vector
        xx=[0 for i in range(self.nvars)]
        for s in f:
            sl=s.split() # space-separated parts
            s1=re.split(r'(\d+)',sl[0]) # 1-st part spli by digits
            if s1[0]=="x" and s1[2]=="_": # fits variable name pattern
                xx[
                    addcols[int(s1[1])]+ # location of x for this beam
                    int(s1[3]) # component of x
                        ]=float(sl[1])
        if qclose:
            f.close()
        return np.array(xx)
    ## Given the vector/list of beamlet intensities,
    #  returns the vector of dose deopsition in voxels
    def dose(self,x):
        return self.d.dot(np.array(x).T)
    ## Given dose deposition vector d (numpy) and doi id (int), 
    #  next functions return dose distribution measures
    #! generalize to "r is int or list of ints"
    def roimax(self,d,r): # max. dose in voxels
        return np.max(d[self.roivox[r]])
    def roimin(self,d,r): # min. dose in voxels
        return np.min(d[self.roivox[r]])
    def roiavg(self,d,r): # average dose in voxels    
        return np.sum(d[self.roivox[r]])/len(self.roivox[r])
    def roidev(self,d,r,d0): # defiation (squares) from given dose d0 
        return np.sqrt(
                np.sum( (d[self.roivox[r]]-d0)**2 )/
                len(self.roivox[r])
                 )

if __name__=="__main__":

#    
    ## Reading patient data
    datadir="C:\\Mytemp\\IMRT\\1212485_start\\"
    p=Patient(datadir+"m_PARETO_3.txt")
    # English translations of ROI
    enroi=["Body","Brainstem","Br.stem 3mm","Spinal cord","Sp. cord 3mm",
       "Salvary L","Salvary R","Jaw","NT",'PTV1 (67.5)','PTV2 (60)','PTV3 (54)']    

    # Solutions data for the patient
    soldir="C:\\Mytemp\Dropbox (IBS PAN)\\IBS PAN Team Folder\\"+\
           "IMRT\\Mini-conf 2019\\DP presentation\\Calc\\"
    sol_fnames=[
        r"From JM\imrt_lp_BarConvTol=0.001.sol",
        r"results LS 24-05-2019\imrt_graniczenia_dla_funkcji_minmaxOARS.sol",
        r"results LS 24-05-2019\imrt_niejednorodnosc_PTV67_50.sol",
        r"results LS 24-05-2019\imrt_niejednorodnosc_PTV60.sol",
        r"results LS 24-05-2019\imrt_niejednorodnosc_PTV54.sol",
        r"results LS 24-05-2019\imrt_suma_niejednorodnosci_PTV.sol"            
            ]
    # names of the solutions for objective functions used
    sol_dispnames=[
            "Avg Salv. L&R",
            "Max in OARs",
            "Deviat. PTV1",
            "Deviat. PTV2",
            "Deviat. PTV3",
            "Dev. all PTVs"
            ]
    # list of ROIs to draw
    lroi=[9,10,11, #PTVs
          2,4,5,6,7,8
          ]
    # list of ROI-functions to display
    lroifun=[
        [9,"dev.",p.roidev,67.5], # deviation PTV1
        [10,"dev.",p.roidev,60], # deviation PTV2           
        [11,"dev.",p.roidev,54], # deviation PTV3
        [5,"avg.",p.roiavg], # average in salvary L
        [6,"avg.",p.roiavg], # average in salvary R
        [4,"max.",p.roimax], # max. in sp. cord + 3mm
        [2,"max.",p.roimax], # max. in br. stem + 3mm
        [7,"max.",p.roimax], # max. in jaw
        [8,"max.",p.roimax] ] # max. in NT

    doses=[]
    for si,sfname in enumerate(sol_fnames):
        xx=p.read_gur_sol(soldir+sfname)
        doses.append(
                [p.dose(xx)[p.roivox[iroi]]
                    for iroi in lroi]
                )
    dvh.dvh2excel(soldir+"out.xlsx",doses,sol_dispnames,[enroi[i] for i in lroi])
    

#    ## table of dose measures values for all solutions x selected ROIs
#    sol_fun_table=pd.DataFrame(
#            columns=sol_dispnames,
#            index=[p.roinames[l[0]]+" "+l[1] for l in lroifun]
#            )
#    for si,sfname in enumerate(sol_fnames):
#        xx=p.read_gur_sol(soldir+sfname)
#        ds=p.dose(xx)
#        for i,l in enumerate(lroifun):
#            sol_fun_table.iloc[i,si]=l[2](ds,l[0],*l[3:])
#    sol_fun_table.to_excel(soldir+"sol_table.xlsx")
#    workbook=xlsxwriter.Workbook(soldir+"out.xlsx")
#    # a sheet with charts for each ROI
#    chartsheet=workbook.add_worksheet("charts")
#    roicharts=[]
#    for i,iroi in enumerate(lroi):
#        roicharts.append(
#            workbook.add_chart({'type': 'scatter',"subtype":"smooth"}))
#        roicharts[-1].set_title({"name":enroi[iroi]})
#        roicharts[-1].set_x_axis({
#                'major_gridlines': {'visible': True},
#                        })
#    # sfname - filename of solution
#    for si,sfname in enumerate(sol_fnames):
#        ## read the solution and create worksheet
#        xx=p.read_gur_sol(soldir+sfname)
#        ds=p.dose(xx)
#        worksheet=workbook.add_worksheet(sol_dispnames[si])
#        chart = workbook.add_chart({'type': 'scatter',"subtype":"smooth"})
#        chart.set_x_axis({
#                'major_gridlines': {'visible': True},
#                        })
#        # iroi - nr. of ROI from the patient list 
#        for i,iroi in enumerate(lroi):
#            # preparing cumulative histogram
#            h=cumfreq(ds[p.roivox[iroi]],100)
#            b=h.lowerlimit  # 1st bin's lower bound
#            vtotal=len(p.roivox[iroi]) # nr. of voxels for scale
#            worksheet.write(0,i*2,enroi[iroi]) # ROI name
#            # writing all diagram elements
#            worksheet.write(1,i*2,0) # dose >=0 for 100%
#            worksheet.write(1,i*2+1,1)
#            worksheet.write(2,i*2,b) # dose <= min. dose for 100%
#            worksheet.write(2,i*2+1,1)            
#            for j,n in enumerate(h.cumcount):
#                b+=h.binsize
#                worksheet.write(j+3,i*2,b)
#                worksheet.write(j+3,i*2+1,(vtotal-n)/vtotal)
#            # adding series to the solution-specific chart
#            chart.add_series({
#            'categories': [sol_dispnames[si],
#                           1, i*2, len(h.cumcount)+2, i*2],
#            'values':     [sol_dispnames[si],
#                           1, i*2+1, len(h.cumcount)+2, i*2+1],
#            'name':       [sol_dispnames[si], 0, i*2],
#            'line': {'width':1.25}
#                })
#            # adding series to the roi-specific chart
#            roicharts[i].add_series({
#            'categories': [sol_dispnames[si],
#                           1, i*2, len(h.cumcount)+2, i*2],
#            'values':     [sol_dispnames[si],
#                           1, i*2+1, len(h.cumcount)+2, i*2+1],
#            'name':       sol_dispnames[si],
#            'line': {'width':1.25}
#                })
#            
#        worksheet.insert_chart("B2",chart)
#    # placing ROI-specific charts on the sheet
#    chartpos=[
#        [i*15,j*8] for i in range(int(len(lroi)/3)+1)        
#            for j in range(3)]
#    for i,ich in enumerate(roicharts):
#        chartsheet.insert_chart(
#                chartpos[i][0],
#                chartpos[i][1],
#                ich)
#    workbook.close()       
    
#### producing megabeamlets from organ structure
#    outdir="C:\\Mytemp\Dropbox (IBS PAN)\\IBS PAN Team Folder\\"+ \
#           "IMRT\\Experiments\\Megabeamlets\\"
#    b0=0 # nr. of the 1st beamlet in the considered beam
#    flog=open(outdir+"log.txt","w")
#    # creating megabeamlet files for each beam
#    for bi in range(p.nbeams):
#        # init: list of lists of beamlet nrs. (starting from 1 in each beam)
#        #       belonging to all possible combinations of ROIs
#        ll=[[] for i in range(2**p.nroi)]
#        # checking which ROIs does a beamlet intersect with
#        for i in range(p.nbeamlist[bi]):
#            # init: nrs of times the beamlet intersects with all rois
#            bitl=[0 for i in range(p.nroi)]
#            # go through indices of intersection with voxels 
#            # (nonzero elements of the beamlet-related column)
#            for j in p.d.getcol(b0+i).nonzero()[0]:
#                for bj,b in enumerate(bitlist(p.voxcodes[j])):
#                    if b>0:
#                        bitl[bj]=1
#            # put beamlet nr to the list corresponding to its ROIs combination
#            ll[bits2int(bitl)].append(i)
#        # writing megabeamlets lists to a file
#        nmega=0 # ounting nr. of obtained megabeamlets
#        with open(outdir+"beam_"+str(bi+1)+".txt","w") as f:
#            print("*** Writing for beam ",bi)
#            flog.write("*** Beam "+str(bi+1)+"\n")
#            for i,l in enumerate(ll):
#                if len(l)>0:
#                    # list of beamlet nrs to a file
#                    f.write(" ".join([str(a) for a in l])+"\n")
#                    # logging information
#                    sout=", ".join([
#                         p.roinames[j] for j,b in 
#                             enumerate(bitlist(i)) if b>0
#                            ])+": "+str(len(l))
#                    flog.write(sout+"\n")
#                    print(sout)
#                    nmega+=1
#            print("*** ",nmega," megabeamlets\n\n")
#        b0+=p.nbeamlist[bi]
#    flog.close()
    
