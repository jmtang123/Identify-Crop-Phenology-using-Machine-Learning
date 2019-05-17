# This file has all functions: AddNdvi,removenn,intense,predata

# function to add NDVI band to image
def clip(image):
    import ee
    ee.Initialize()
    bound = ee.Geometry.Rectangle([89.985917,25.019917,90.799167,23.498889])
    return image.clip(bound)

def AddNdvi(image):
    red=image.select('sur_refl_b01')
    nir=image.select('sur_refl_b02')
    ndvi = (nir.subtract(red)).divide(nir.add(red)).rename('ndvi')
    return image.addBands(ndvi)

# turn the Pandas timestamp to day of year
def dayyear(list):
    dayn=[0]*len(list)
    for d in range(len(list)):
        dayn[d] = list[d].to_pydatetime().timetuple().tm_yday
    return dayn

# remove null and negative values from NDVI
def removenn(list):
    x=[0]*len(list)
    for i in range(len(list)):
        value=list[i]
        if value == None or value < 0:
            if i!= (len(list)-1):
                y=list[i+1]
                if y==None or y<0:
                    x[i]=list[i-1]
                    list[i]=list[i-1]
                else:
                    x[i]=(list[i-1]+list[i+1])/2
            else:
                x[i]=list[i-1]
        else:
            x[i]=list[i]
    return x

# detect the intensity and local maximum NDVI (highest value existing within 3 month)
def intense(list):
    intensity=0
    maxv=[]
    maxindex = []
    slist=[0]*2
    minv=[0]*2
    minindex=[0]*2
    n=len(list)
    for j in range(5,n-5):
        newlist=[list[j-5],list[j-4],list[j-3],list[j-2],list[j-1],list[j],list[j+1],list[j+2],list[j+3],list[j+4],list[j+5]]
        if list[j]==max(newlist) and list[j]>0.5:
            intensity=intensity+1
            print('a')
            maxv.append(list[j])
            maxindex.append(j)
    if intensity>3:
        print("There is something wrong in intensity")
    elif intensity==2:
        slist[0]=list[maxindex[0]:maxindex[1]+1]
        minv[0]=min(slist[0])
        minindex[0]=maxindex[0]+slist[0].index(min(slist[0])) # the minindex = first max + min location within small list
    elif intensity==3:
        slist[0] = list[maxindex[0]:maxindex[1] + 1]
        slist[1] = list[maxindex[1]:maxindex[2] + 1]
        minv[0] = min(slist[0])
        minv[1] = min(slist[1])
        minindex[0] = maxindex[0] + slist[0].index(min(slist[0]))
        minindex[1] = maxindex[1] + slist[1].index(min(slist[1]))
    return intensity,maxv,maxindex,minv,minindex

# return the ndvi data into the y as sigmoid function
def predata(list):
    mindata=min(list)
    maxdata=max(list)
    newlist=[0]*len(list)
    for i in range(len(list)):
        if list[i]==min(list):
            list[i]=list[i]+0.001
        elif list[i]==max(list):
            list[i]=list[i]-0.001
        newlist[i]=(list[i]-mindata)/(maxdata-mindata)
    return newlist

# sigmoid function used for the curve-fitting
def sigmoid(x, a, b):
    import numpy as np
    y = 1/(1+np.exp(a+b*x))
    return y

# curve-fitting
def curvefit(list1,list2):
    from scipy.optimize import curve_fit
    if list2[-1]>0.5:
        popt, pcov = curve_fit(sigmoid, list1, list2,p0=(20,-0.1))
    else:
        popt, pcov = curve_fit(sigmoid, list1, list2, p0=(-20, 0.1))
    return popt,pcov

# generate the ndvi date for visualization/validation
def validate(list1,x1,x2,list2):
    import numpy as np
    a=list2[0]
    b=list2[1]
    y=[0]*len(list1)
    for i in range(len(list1)):
        y[i]=x1+(x2-x1)*(1/(1+np.exp(a+b*list1[i])))
    return y

# calculate the change rate of curvature
# and then find the four days (greenup date GD, maturity date MD, senescence date SD, dormancy DD)
def curverate(list1,x1,x2,list2):
    import numpy as np
    a=list2[0]
    b=list2[1]
    c=x2-x1
    allday=range(list1[0],list1[-1]+1)
    ccr=[0]*len(allday)
    y=[0]*len(allday)
    greenD=0
    matureD=0
    senesD=0
    dormanD=0
    for i in range(len(allday)):
        z=np.exp(a+b*allday[i])
        y[i]=(c/(1+z))+x1
        s1=3*z*(1-z)*((1+z)**3)*(2*((1+z)**3)+(b**2)*(c**2)*z)
        s2=(((1+z)**4)+((b*c*z)**2))**2.5
        s3=((1+z)**2)*(1+2*z-5*(z**2))
        s4=(((1+z)**4)+((b*c*z)**2))**1.5
        ccr[i]=((b**3)*c*z*((s1/s2)-(s3/s4)))
    if b<0:
        top2index=sorted(np.argsort(ccr)[-2:])
    else:
        top2index = sorted(np.argsort(ccr)[:2])
    detectd1=[allday[i] for i in top2index][0]
    detectd2=[allday[i] for i in top2index][1]
    return y,ccr,detectd1,detectd2,allday
