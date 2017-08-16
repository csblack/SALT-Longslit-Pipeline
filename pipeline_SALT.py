######********************************************
#Pipeline written by Christine Black.  v. 8/16/17

#This code uses PYRAF to go through SALT spectra and perform bias subtraction, trimming of any overscan regions, combinging and normalizing the flats, flat fielding the data, and cosmic ray removal

#This code expects to have all the spectra data (including the flat field images) in the same folder

#After each step in the reduction, the files will be appended with a letter denoting what has been done (b=bias,t=trim,f=flat,cr=cosmic ray removal)

#Flat normalization requires some interaction from the user

#In order for cosmic ray removal to work, ***this code requires that LA Cosmic (cosmics.py) be installed*** If LA Cosmic is not installed, comment out the CR Removal section.  Make sure the path (In the import section) is correct for the current computer

#From the directory with the data the command is 'python pipeline_SALT.py' to run the code

#The code assumes that you are using the SALT Longslit.  If you are not using these you will need to change the code accodringly.  Comments are located where these changes need to be made
######********************************************

######
#Importing Tasks
######

import glob, os
from pyraf import iraf
from iraf import images,noao,imred,ccdred,twodspec,longslit,apextract,onedspec
from iraf import hselect,bias,linebias,imarith,imdelete,imcopy,flatcombine,response,ccdproc,apall

# The following two lines are only needed as cosmic.py is not in this directory nor in the python path.  They would not be required if you copy cosmics.py in this directory.

###Path should be updated for your computer!!!
import sys
sys.path.append("/Users/Hogan/cosmics.py_0.4") # The directory that contains cosmic.py
import cosmics


######********************************************
#Reducing SALT Longslit spectra
######********************************************

######
#Puts data in list and file
######

file_list = glob.glob("*.fits")
dataraw_list=[]
text_file = open("dataraw.txt","wr")

for f in file_list:
    temp=f.split('.')
    g=temp[0]+'.fits'+'[SCI]'
    dataraw_list.append(g)


for f in range(0,len(dataraw_list)):
    if f < len(dataraw_list)-1:
        text_file.write(dataraw_list[f]+"\n")
    else:
        text_file.write(dataraw_list[f])
        
text_file.close()


######
#Clear CCDPROC
######
ccdproc.noproc="no"
ccdproc.fixpix="no"
ccdproc.overscan="no"
ccdproc.trim="no"
ccdproc.zerocor="no"
ccdproc.darkcor="no"
ccdproc.flatcor="no"
ccdproc.illumcor="no"
ccdproc.fringecor="no"
ccdproc.readcor="no"
ccdproc.scancor="no"


######
#Define Variables
######
data_list=[]
flat_list=[]
data_bt_list=[]
data_btf_list=[]
data_btfcr_list=[]
data_btfcr_comp_list=[]
data_btfcr_data_list=[]
data_forCR_list=[]
##
rdnoise = 23.0
gain = 1.0
##Change values if NOT using SALT Longslit#
mybias1="[20:1000,5:20]"
mytrim1="[1:1021,40:995]"

mybias2="[1080:2090,5:20]"
mytrim2="[1074:2097,40:995]"

mybias3="[2160:3140,5:20]"
mytrim3="[2146:3166,40:995]"

##
dispax = 1 # 1=horizontal image, 2=vertical image
##


######
#Subtract Bias
######
#This splits up each part of the chip into individual pieces and uses the overscan region below each part as the bias for each.  Then the parts are put back together
print 'Subtracting Bias (produces name.bt.fits) ...'

for f in dataraw_list:
    temp=f.split('.')

    linebias(f,str(temp[0])+'.1', bias=mybias1, trim=mytrim1 ,function="legendre",order=1, interactive="no")
    linebias(f,str(temp[0])+'.2', bias=mybias2, trim=mytrim2 ,function="legendre",order=1, interactive="no")
    linebias(f,str(temp[0])+'.3', bias=mybias3, trim=mytrim3 ,function="legendre",order=1, interactive="no")
    
    g=str(temp[0])+'.bt'+'.fits'
    
    #Adjust values if NOT using SALT Longslit#
    imarith(f+'[1:3166,40:995]','*',1.0 ,g)
    print 'Making final image...'
   
    f1=str(temp[0])+'.1'+'.fits'
    f2=str(temp[0])+'.2'+'.fits'
    f3=str(temp[0])+'.3'+'.fits'
    #Adjust values if NOT using SALT Longslit#
    imcopy(f1,g+'[1:1021,1:956]')
    imcopy(f2,g+'[1074:2097,1:956]')
    imcopy(f3,g+'[2146:3166,1:956]')

    print 'Deleting sections ...'
    imdelete(f1)
    imdelete(f2)
    imdelete(f3)
    
    data_bt_list.append(g)
    print 'Next Image'

#This selects all data that needs to be flat fielded and puts them in an array/file
for f in data_bt_list:
    if hselect(f,"$I","OBJECT!='FLAT'",Stdout=1):
        data_list.append(f)

text_file1 = open("data.bt.txt","wr")
text_file2 = open("data.btf.txt","wr")
for f in range(0,len(data_list)):
    temp=data_list[f].split('.')
    data_btf_list.append(str(temp[0])+'.'+str(temp[1])+'f'+'.fits')
    if f < len(data_list)-1:
        text_file1.write(data_list[f]+"\n")
        text_file2.write(str(temp[0])+'.'+str(temp[1])+'f'+'.fits'+"\n")
        
    else:
        text_file1.write(data_list[f])
        text_file2.write(str(temp[0])+'.'+str(temp[1])+'f'+'.fits')
        
text_file1.close()
text_file2.close()

######
#Combining/Normalizing Flats
######
print "Combining Flats (produces Flat) ..."
#This looks for all the flats and puts them in an array/file   
for f in data_bt_list:
    if hselect(f,"$I","OBJECT=='FLAT'",Stdout=1):
        flat_list.append(f)
print flat_list
text_file = open("flat_data.txt","wr")
for f in range(0,len(flat_list)):
    if f < len(flat_list)-1:
        text_file.write(flat_list[f]+"\n")
    else:
        text_file.write(flat_list[f])
        
text_file.close()

#This combines all the flats into one Flat.fits
flatcombine("@%s.txt" %'flat_data',output='Flat', combine="median", reject="avsigclip", gain=gain, rdnoise=rdnoise)

print "Normalize Flat (produces Flat.n) ..."
#Adjust values if NOT using SALT Longslit#
longslit.dispaxis = dispax

#This will run through the Flat normalization (this requires some interaction) and outputs Flat.n.fits
while True:
    response(calibration='Flat', normalization='Flat',response='Flat.n',interactive='yes',function='spline3')
    while True:
        ans=raw_input("Take a look at Flat.n.  Are you happy with it? ")
        if (ans == 'yes') or (ans == 'y'):
            break
        elif (ans == 'no') or (ans == 'n'):
            os.remove("Flat.n.fits")
            break
        elif (ans != 'no') or (ans != 'n') or (ans != 'yes') or (ans != 'y'):
            print "yes or no"
    if (ans == 'yes') or (ans == 'y'):
            break 


######
#Divides by Normalized Flat
######
#This divides all the data by the flat
print "Dividing Flat (produces name.btf.fits)"
#ccdproc("@%s.txt" %'data.bt_data',output="@%s.txt" %'data.btf_data',flatcor='yes',flat='Flat.n')
imarith("@%s.txt" %'data.bt', '/', "%s" %'Flat.n',"@%s.txt" %'data.btf')



######
#Remove Cosmic Rays using LA COSMIC
######
#This will run LA Cosmic on the Data and Standards (and not the arcs)

#Removes Comp Lamps from CR reduce
for f in data_btf_list:
    if hselect(f,"$I","OBJECT!='ARC'",Stdout=1):
        data_forCR_list.append(f)

text_file = open("data.btf_Forcr.txt","wr")
for f in range(0,len(data_forCR_list)):
    temp=data_forCR_list[f].split('.')
    if f < len(data_forCR_list)-1:
        text_file.write(str(temp[0])+'.'+str(temp[1])+'.fits'+"\n")
        
    else:
        text_file.write(str(temp[0])+'.'+str(temp[1])+'.fits')
        
text_file.close()

        
print "Running LA Cosmic (produces name.btf.cr.fits)"
#This is the code for LA Cosmic that has been copied here
text_file = open("data.btf.cr.txt","wr")
for f in range(0,len(data_forCR_list)):
    temp=data_forCR_list[f].split('.')
    if f < len(data_forCR_list)-1:
        text_file.write(str(temp[0])+'.'+str(temp[1])+'.cr'+'.fits'+"\n")
        
    else:
        text_file.write(str(temp[0])+'.'+str(temp[1])+'.cr'+'.fits')
        
text_file.close()


text_file = open("data.btf_Forcr.txt","r")
lines=text_file.readlines()

# Read the FITS :
for l in lines:
    array, header = cosmics.fromfits(l)
    # array is a 2D numpy array
    
    # Build the object :
    c = cosmics.cosmicsimage(array, gain=gain, readnoise=rdnoise, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
    # There are other options, check the manual...

    # Run the full artillery :
    c.run(maxiter = 4)
    temp=l.split('.')
    #data_btfcr_list.append(str(temp[0])+'.'+str(temp[1])+'.'+str(temp[2])+'.cr'+'.fits')
    # Write the cleaned image into a new FITS file, conserving the original header :
    cosmics.tofits(str(temp[0])+'.'+str(temp[1])+'.cr'+'.fits', c.cleanarray, header)

    # If you want the mask, here it is :
    cosmics.tofits("mask.fits", c.mask, header)
    # (c.mask is a boolean numpy array, that gets converted here to an integer array)





print "\n"+"Done!"

