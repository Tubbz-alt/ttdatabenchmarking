#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel7-gcc48-opt/bin/python

import psana;
import numpy as np;
#import matplotlib.pyplot as plt;
import TimeTool;
import argparse;

# AMO amo11816: 36 (490:495)
# XPP 12mu water runs: 106 (95:100) = 10%, 107 = 1%, 109 = 50%, 110 (95:100) = 100%
# XPP 20mu water runs: 99 (105:110) = .8%, 100 = .1%, 101 = 10%, 
# XPP diamond runs: 17 (205:225) = 30%, 20 = 50%, 21 = 100%, ?? this seems to fail 44 = 100%+slottedfoil10fs
# XPP SiN 1mu runs: 54 (205:225) = 5%, 58 = 50%, 59 = 100%, 64 (205:225) = 100%
# XPP SiN 2mu runs: 102 (105:110) = 1%, 103 = 10%, 104 (105:110) = 5%
limstr="Signal vert lims:\tr17-64 = (205:225)\tr99-r104 = (105:110)\tr106-r110 = (95:100)"
helpstr = 'Compute Matrix and shot sorting\tfor XPP interferometric method intensity scans.\t' + limstr;
parser = argparse.ArgumentParser(description=helpstr);#'Compute Matrix and shot sorting\n\tfor XPP interferometric method intensity scans');

parser.add_argument('-e','--exp',dest='expstr', type=str, help='experiment identifier', default='xppl3816');
parser.add_argument('-r','--run',dest='runstr', type=str, help='run number', default='21');
parser.add_argument('-d','--dirname',dest='dirname',type=str, help='directory name for data',default='../data/processed/');
parser.add_argument('-l','--llim',dest='llim',type=int,help='signal vrange lower limit (inclusive)',default=int(205));
parser.add_argument('-u','--ulim',dest='ulim',type=int,help='signal vrange upper limit (exclusive)',default=int(225));
parser.add_argument('-n','--skipshots',dest='skipshots',type=int,help='skip n shots',default=int(40));
parser.add_argument('-s','--skipsteps',dest='skipsteps',type=int,help='skip n steps',default=int(10));
parser.add_argument('-a','--attenuation',dest='atten',type=float,help='attenuation (actually transmission e.g. .1 xrays on sample = .1',default=float(1));
parser.add_argument('--subref',dest='subref',type=bool,help='bool to subtract reference',default=False);

args = parser.parse_args();

def gd2mj(gdets):
	# we are standardizing on the GDET 22 detector
	m43 = 1.10982; # from fitting
	m23 = 1.06575; # from fitting
	m53coefs = (1.0,414.849,-3344.54,11030.7,-16317.8,9119.7);
	m53pows = (1.0,5,6,7,8,9);
#g(x)=a*x+e*x**5+f*x**6+g*x**7+h*x**8+i*x**9
# a=1
#Final set of parameters            Asymptotic Standard Error
#=======================            ==========================
#e               = 414.849          +/- 151.6        (36.54%)
#f               = -3344.54         +/- 1230         (36.77%)
#g               = 11030.7          +/- 3723         (33.75%)
#h               = -16317.8         +/- 4986         (30.56%)
#i               = 9119.7           +/- 2493         (27.34%)
#
#correlation matrix of the fit parameters:
#                e      f      g      h      i      
#		e               1.000 
#		f              -0.999  1.000 
#		g               0.996 -0.999  1.000 
#		h              -0.991  0.996 -0.999  1.000 
#		i               0.986 -0.993  0.997 -0.999  1.000 
	vec=[];#np.array(dtype=float);
	for gd in gdets:
		vec.append(gd);
	mj = vec[3];
	mj += m23*vec[2];
	mj += m43*vec[4];
	x = vec[5]*np.ones(len(m53pows));
	mj += np.sum(m53coefs*np.power(x,m53pows));
	return mj/4.;

def tt_i2ps(i):
	ps = 0.100/(.025*(i-365)+34.5);
	sign = float(-1); # for XPP we need to subtract the TimeTool delay
	return sign*ps;



def i2lam(i):
        lset=550;
        nmPi=0.217;
        seterr=1.0051;
        return nmPi*i + seterr*lset - 110.072;


headstr = '#';
#'i' is the pixel index [ 0 .. 1023 ] and the wavelength [nm]
#runstr = args.runstr;
#expstr = str('xppl3816');
dsourcestr = 'exp=' + args.expstr + ':run=' + args.runstr;
print(dsourcestr)

# setting up TimeTool
bykick = int(162);
ttOptions = TimeTool.AnalyzeOptions(get_key='opal_1', eventcode_nobeam = bykick, sig_roi_x = '1 1022', sig_roi_y = '350 375',sb_roi_x = '1 1022',sb_roi_y = '150 175');
ttAnalyze = TimeTool.PyAnalyze(ttOptions);

ds = psana.DataSource(dsourcestr,module=ttAnalyze);
print(psana.DetNames('detectors'));


evr = psana.Detector('NoDetector.0:Evr.0')
det = psana.Detector('opal_0'); #'XppEndStation.0:Opal1000.0') # for XPP
TTdet = psana.Detector('opal_1'); #This is the TimeTool camera I think
cd = psana.Detector('ControlData')
EBdet = psana.Detector('EBeam'); #This I hope is the BLD data
GDdet = psana.Detector('FEEGasDetEnergy'); #This I hope is the FEE Gas Detector readings

# grabbing diode detectors for calibration #
inten_dets = (psana.Detector('NH2-SB1-IPM-01'), psana.Detector('XppMon_Pim0'), psana.Detector('XppMon_Pim1'), psana.Detector('XppSb2_Ipm'),psana.Detector('XppSb3_Ipm'),psana.Detector('XppEnds_Ipm0'));

num = 0.0;
y_init = 0;
y_final = 0;

vwin = (args.llim,args.ulim);
num = vwin[1]-vwin[0];

printsample = False;
ratio = .1;
#subref = False;

delayscale=float(1e12); #XPP



'''The third edition: take the average of each step and convert both axes to the right units. '''
for run in ds.runs():
	reference_img = np.zeros(1024,dtype=float);
	y = np.zeros(1024,dtype=float);
	nrefshots = int(0);
	for nstep,step in enumerate(run.steps()):
		printsample = nstep%20==0;
		if nstep%args.skipsteps==0:
			pvList = cd().pvControls();
			for pv in pvList:
				if y_init == 0:
					y_init = pv.value()
				y_final = pv.value()	
				print('Step', nstep, 'name/value',pv.name(),pv.value());
			for nevent,evt in enumerate(step.events()):
				ttResults = ttAnalyze.process(evt);
				ebResults = EBdet.get(evt);
				gdResults = GDdet.get(evt);

    				img = det.image(evt);  #img is not areal image, it is a matrix

				if (img is None):
					continue;

				if (printsample and (nevent+250)%500==0):
					filename = args.dirname + args.expstr + '_r' + args.runstr + '_image%d.step%d.dat' % (nevent,nstep);
					print('printing image ', filename);
					headstr='# sample image'
					np.savetxt(filename,img,fmt='%.6e',header=headstr);
					printsample = False;
					filename = args.dirname + args.expstr + '_r' + args.runstr + '_TTimage%d.step%d.dat' % (nevent,nstep);
					print('and image ' , filename);
					ttimg = TTdet.image(evt);
					headstr='# sample image TimeTool camera'
					np.savetxt(filename,ttimg,fmt='%.6e',header=headstr);
					if ttResults!=None:
						print(ttResults.position_time());

				ec = evr.eventCodes(evt)
				if bykick in ec: 
					y = np.sum(img[vwin[0]:vwin[1],:],axis=0)/num;
					if nrefshots==0:
						reference_img = y;#np.sum(img[vwin[0]:vwin[1],:],axis=0)/num;
					else:
						reference_img *= (1.-ratio);
						reference_img += (ratio)*y;#np.sum(img[vwin[0]:vwin[1],:],axis=0)/num;
					nrefshots += 1;
					refsum = sum(reference_img);
					
					try:
						Y = np.row_stack((Y,y));
						YS = np.row_stack((YS,(nstep,delayscale*y_final)));

					except NameError:
						Y=y;
						YS=(nstep,delayscale*y_final);

					continue;

				lineout = np.zeros(1024,dtype=float);
				if nevent%args.skipshots == 0:
					try:
						lineout = (np.sum(img[vwin[0]:vwin[1],:],axis=0)/num) ;
						tt_data = (ttResults.position_pixel(),ttResults.amplitude(),ttResults.position_time(),ttResults.position_fwhm());
						dd_data = (nstep,delayscale*y_final,delayscale*y_final+tt_i2ps(tt_data[0]));
						eb_data = (ebResults.ebeamL3Energy() , ebResults.ebeamCharge(), ebResults.ebeamEnergyBC1(), ebResults.ebeamEnergyBC2(), ebResults.ebeamLTU250(), ebResults.ebeamLTU450(), ebResults.ebeamLTUAngX(), ebResults.ebeamLTUAngY(), ebResults.ebeamLTUPosX(), ebResults.ebeamLTUPosY(), ebResults.ebeamUndAngX(), ebResults.ebeamUndAngY(), ebResults.ebeamUndPosX(), ebResults.ebeamUndPosY(), ebResults.ebeamPkCurrBC1(), ebResults.ebeamEnergyBC1(), ebResults.ebeamPkCurrBC2(), ebResults.ebeamEnergyBC2(), ebResults.ebeamDumpCharge());
						gd_data = ( gdResults.f_11_ENRC(), gdResults.f_12_ENRC(), gdResults.f_21_ENRC(), gdResults.f_22_ENRC(), gdResults.f_63_ENRC(), gdResults.f_64_ENRC() );
						ii_data = np.zeros(6*len(inten_dets),dtype=float) ;
						for i,intendet in enumerate(inten_dets):
							ii_data[i*6:i*6+4] = intendet.channel(evt).tolist();
							ii_data[i*6+4] = intendet.xpos(evt);
							ii_data[i*6+5] = intendet.ypos(evt);


					except:
						print("Something failed in filling event data values");
						continue;
					inten=args.atten*gd2mj(gd_data);
					sig_data = (inten,sum(lineout)/1024,(sum(lineout)-refsum)/1024);
					#print(sig_data);
					if args.subref:
						lineout = lineout - reference_img[:len(lineout)];
					try:

						R = np.row_stack((R,lineout));#it stacks from top to bottom
						D = np.row_stack((D,dd_data));
						T = np.row_stack((T,tt_data));
						E = np.row_stack((E,eb_data));
						G = np.row_stack((G,gd_data));
						S = np.row_stack((S,sig_data));
						I = np.row_stack((I,ii_data));
					
					except NameError:
						R = lineout;
						D = dd_data;
						T = tt_data;
						E = eb_data;
						G = gd_data;
						S = sig_data;
						I = ii_data;

lam = i2lam(np.arange(R.shape[1],dtype=float));


#for plot
#y_dim = int(np.shape(R)[0]);
#x_dim = int(np.shape(R)[1]);
#delay = delayscale*np.linspace(y_init,y_init,y_dim,dtype=float);

sorted_inds = np.argsort(D[:,1]);

if args.subref:
	runstr += "_refsub";



#filename = args.dirname;#'data_tmp/';
#filename+=expstr + '_r' + runstr + '_matrix_sorted.dat';
#headstr='# sorted matrix';
#np.savetxt(filename,R[sorted_inds,:],fmt='%.6e');


filename = args.dirname + args.expstr + '_r' + args.runstr + '_inds_sorted.dat';
headstr='index ordered to sort the raw matrix';
np.savetxt(filename,sorted_inds,fmt='%.6e',header=headstr);

filename = args.dirname + args.expstr + '_r' + args.runstr + '_justinrefs.dat';
headstr = 'matrix of bykick shots';
np.savetxt(filename,Y,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_justinrefsdelays.dat';
headstr = 'delays for the matrix of bykick shots';
np.savetxt(filename,YS,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_matrix.dat';
headstr = 'matrix of shots';
np.savetxt(filename,R,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_delays.dat';
headstr = 'nstep\tdelay\ttt_corrected delay\t';
np.savetxt(filename,D,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_signal.dat';
headstr = 'computed_mJ\tmeanlineout\tmean-refmean';
np.savetxt(filename,S,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_eb.dat';
headstr = 'ebeamL3Energy\tebeamCharge\tebeamEnergyBC1\tebeamEnergyBC2\tebeamLTU250\tebeamLTU450\tebeamLTUAngX\tebeamLTUAngY\tebeamLTUPosX\tebeamLTUPosY\tebeamUndAngX\tebeamUndAngY\tebeamUndPosX\tebeamUndPosY\tebeamPkCurrBC1\tebeamEnergyBC1\tebeamPkCurrBC2\tebeamEnergyBC2\tebeamDumpCharge()';
np.savetxt(filename,E,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_gd.dat';
headstr = 'f_11_ENRC\tf_12_ENRC\tf_21_ENRC\tf_22_ENRC\tf_63_ENRC\tf_64_ENRC';
np.savetxt(filename,G,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_tt.dat';
headstr = 'pos\tamp\ttime\tfwhm';
np.savetxt(filename,T,fmt='%.6e',header=headstr);
filename = args.dirname + args.expstr + '_r' + args.runstr + '_ii.dat';
headstr = 'each:ch0\tch1\tch2\tch3\tx\ty\t...'; # + '\t'.join(inten_dets);
np.savetxt(filename,I,fmt='%.6e',header=headstr);

filename = args.dirname + args.expstr + '_r' + args.runstr + '_wavelegths.dat';
headstr = 'wavelength';
np.savetxt(filename,lam,fmt='%.6e',header=headstr);

print('Done saving');

#plt.imshow(R,origin = 'lower',extent = [lam[0],lam[-1],delay[0],delay[-1]],aspect = 'auto')
#tick_locs_x = np.linspace(450,650,5)
#plt.xticks(tick_locs_x)

#tick_locs_y = np.linspace(delay[0],delay[-1],6)
#plt.yticks(tick_locs_y)

#plt.hot()
#plt.colorbar()
#plt.show()				
print('Done.');
