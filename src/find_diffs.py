#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel7-gcc48-opt/bin/python
###!/reg/neh/home/coffee/miniconda3/bin/python

import numpy as np
import argparse



parser = argparse.ArgumentParser(description='Compute Matrix of differnces for opal 1 calibration\n\tfor XPP interferometric method intensity scans');
parser.add_argument('-r','--run',dest='runstr', type=str, help='run number', default='110');
parser.add_argument('-e','--exp',dest='expstr', type=str, help='expname', default='xppl3816');
parser.add_argument('-d','--dirname',dest='dirname',type=str, help='directory name for data',default='data_tmp/');
parser.add_argument('-p','--precision',dest='precision',type=float, help='precision for the inner product',default=1e-5);


args = parser.parse_args();

def goodampl(x,y):
	m=0.085;
	X=np.sum(x,axis=1)/float(x.shape[1]);
	return np.where(y < m*X)[0];
def tt_i2ps(x):
	ps = 0.100/(.025*(x-365)+34.5);
	return ps;

datafile = args.dirname + args.expstr + '_r' + args.runstr + '_matrix.dat'
print('loading ',datafile);
data = np.loadtxt(datafile,dtype=float);
print('data.shape = ',data.shape);
datafile = args.dirname + args.expstr + '_r' + args.runstr + '_delays.dat'
dl_data = np.loadtxt(datafile,dtype=float);
datafile = args.dirname + args.expstr + '_r' + args.runstr + '_tt.dat'
tt_data = np.loadtxt(datafile,dtype=float);
datafile = args.dirname + args.expstr + '_r' + args.runstr + '_gd.dat'
gd_data = np.loadtxt(datafile,dtype=float);
datafile = args.dirname + args.expstr + '_r' + args.runstr + '_eb.dat'
eb_data = np.loadtxt(datafile,dtype=float);

goodinds = goodampl(gd_data[:,2:3],tt_data[:,1]);
data_keep = data[goodinds,:];
dl_keep = dl_data[goodinds,:];
tt_keep = tt_data[goodinds,:];
gd_keep = gd_data[goodinds,:];
eb_keep = eb_data[goodinds,:];
print('data_keep.shape = ',data_keep.shape);


cov = np.dot(data_keep,data_keep.T);
normvec=np.sqrt(np.diag(cov));
normmat=np.triu(cov/np.outer(normvec,normvec),k=1);
print('normmat.shape = ',normmat.shape);
threshold=(1.0-args.precision);
coords=np.argwhere(normmat>threshold);
print('coords.shape = ',coords.shape);
filename = args.expstr + '_r' + args.runstr + '_coords.out';
np.savetxt(filename,coords,fmt='%i');

print(tt_keep[coords[:10,1],0].reshape((10,1)));

outmat0 = np.zeros((coords.shape[0],data_keep.shape[1]),dtype=float);
outmat1 = np.zeros((coords.shape[0],data_keep.shape[1]),dtype=float);
outdiffsheader = "#diffTTind\tavgTTind\tdiffDL\tavgDL";
outdiffs = np.zeros((coords.shape[0],4),dtype=float);
#outdiffs[:,0]=tt_keep[coords[:,1].reshape((coords.shape[0],1)),0]-tt_keep[coords[:,0].reshape((coords.shape[0],1)),0];
outdiffs[:,0]=tt_keep[coords[:,1],0]-tt_keep[coords[:,0],0];
outdiffs[:,1]=(tt_keep[coords[:,1],0]+tt_keep[coords[:,0],0])/2.;
outdiffs[:,2]=dl_keep[coords[:,1],1]-dl_keep[coords[:,0],1];
outdiffs[:,3]=(dl_keep[coords[:,1],1]+dl_keep[coords[:,0],1])/2.;

outmat0=data_keep[coords[:,0],:];
outmat1=data_keep[coords[:,1],:];

diffsheadstr = "# presicion for match = 1 - %.2e" % args.precision;
filename = args.expstr + '_r' + args.runstr + '_diffs.out';
np.savetxt(filename,outdiffs,fmt='%.6f',header=diffsheadstr);
filename = args.expstr + '_r' + args.runstr + '_outmat0.out';
np.savetxt(filename,outmat0,fmt='%.6f');
filename = args.expstr + '_r' + args.runstr + '_outmat1.out';
np.savetxt(filename,outmat1,fmt='%.6f');

out=np.zeros((60,2),dtype=float)
out[:,1],bins = np.histogram(outdiffs[:,0]*(200./70),bins=60,range=(-15,15));
out[:,0] = bins[:-1];
filename = args.expstr + '_r' + args.runstr + '_hist0.out';
np.savetxt(filename,out,fmt='%i');

out[:,1],bins = np.histogram(outdiffs[:,0],bins=60,range=(20,50));
out[:,0] = bins[:-1];
filename = args.expstr + '_r' + args.runstr + '_hist1.out';
np.savetxt(filename,out,fmt='%i');

out[:,1],bins = np.histogram(outdiffs[:,0],bins=60,range=(55,85));
out[:,0] = bins[:-1];
filename = args.expstr + '_r' + args.runstr + '_hist2.out';
np.savetxt(filename,out,fmt='%i');

out[:,1],bins = np.histogram(outdiffs[:,0],bins=60,range=(75,105));
out[:,0] = bins[:-1];
filename = args.expstr + '_r' + args.runstr + '_hist3.out';
np.savetxt(filename,out,fmt='%i');

out[:,1],bins = np.histogram(outdiffs[:,0],bins=60,range=(125,155));
out[:,0] = bins[:-1];
filename = args.expstr + '_r' + args.runstr + '_hist4.out';
np.savetxt(filename,out,fmt='%i');

raise SystemExit(0);
