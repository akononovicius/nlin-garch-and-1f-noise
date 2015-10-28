import java.io.*;
import java.math.*;
import java.math.BigDecimal;
import java.util.Random;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class garch {
	public static void getCommandLine(String[] str) {
		if(!(str.length>0)) {
			System.out.println("Initiating with default settings...");
			return ;
		}
		int i=0;
		try {
			for(i=0;i<str.length;i++) {
				if(str[i].equalsIgnoreCase("--core") || str[i].equalsIgnoreCase("-c")) {
					i++;
					commonVariables.core=Integer.parseInt(str[i]);
				} else if(str[i].equalsIgnoreCase("--realizations") || str[i].equalsIgnoreCase("-r")) {
					i++;
					commonVariables.realizations=Integer.parseInt(str[i]);
				} else if(str[i].equalsIgnoreCase("--z") || str[i].equalsIgnoreCase("-z")) {
					commonVariables.var=1;
				} else if(str[i].equalsIgnoreCase("--jobs") || str[i].equalsIgnoreCase("-j")) {
					i++;
					commonVariables.core=Integer.parseInt(str[i]);
					i++;
					commonVariables.realizations=Integer.parseInt(str[i]);
				} else if(str[i].equalsIgnoreCase("--realizationPoints") || str[i].equalsIgnoreCase("-rp")) {
					i++;
					int pts=Integer.parseInt(str[i]);
					pts=(int)Math.pow(2,Math.ceil(Math.log(pts)/Math.log(2)));
					commonVariables.realizationPoints=pts;
				} else if(str[i].equalsIgnoreCase("--outPoints") || str[i].equalsIgnoreCase("-op")) {
					i++;
					commonVariables.outputPoints=Integer.parseInt(str[i]);
					commonVariables.multiOutputPoints=commonVariables.outputPoints*30;
				} else if(str[i].equalsIgnoreCase("--gatherPoints") || str[i].equalsIgnoreCase("-gp")) {
					i++;
					commonVariables.gatherPoints=Integer.parseInt(str[i]);
					if(commonVariables.gatherPoints<=0) commonVariables.gatherPoints=1;
				} else if(str[i].equalsIgnoreCase("-dt")) {
					i++;
					commonVariables.dt=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-dt60")) {
					i++;
					commonVariables.dt60=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-st")) {
					i++;
					commonVariables.dt=Double.parseDouble(str[i]);
					commonVariables.dt60=commonVariables.dt;
				} else if(str[i].equalsIgnoreCase("--outputTemplate") || str[i].equalsIgnoreCase("-ot")) {
					i++;
					commonVariables.outputTemplate=str[i];
				} else if(str[i].equalsIgnoreCase("--outputTime") || str[i].equalsIgnoreCase("--outTime") || str[i].equalsIgnoreCase("-outTime")) {
					commonVariables.outputTime=true;
				} else if(str[i].equalsIgnoreCase("--useAbs")) {
					commonVariables.useAbs=true;
				} else if(str[i].equalsIgnoreCase("--noLimits")) {
					commonVariables.useLimits=false;
				} else if(str[i].equalsIgnoreCase("--nonlinear") || str[i].equalsIgnoreCase("-nlin")) {
					commonVariables.linear=false;
				} else if(str[i].equalsIgnoreCase("--nonLinearFeedback") || str[i].equalsIgnoreCase("-nlinfb")) {
					commonVariables.nonLinearFeedback=true;
				} else if(str[i].equalsIgnoreCase("--absoluteNoise") || str[i].equalsIgnoreCase("-absnoise")) {
					commonVariables.absoluteNoise=true;
				} else if(str[i].equalsIgnoreCase("--limits") || str[i].equalsIgnoreCase("-lim")) {
					i++;
					commonVariables.min=Double.parseDouble(str[i]);
					i++;
					commonVariables.max=Double.parseDouble(str[i]);
					commonVariables.useLimits=true;
				} else if(str[i].equalsIgnoreCase("--pdfLimits") || str[i].equalsIgnoreCase("-pdflim")) {
					i++;
					commonVariables.pdfMin=commonFunctions.LogBase10(Double.parseDouble(str[i]));
					i++;
					commonVariables.pdfMax=commonFunctions.LogBase10(Double.parseDouble(str[i]));
				} else if(str[i].equalsIgnoreCase("--pdfLogLimits") || str[i].equalsIgnoreCase("-pdfllim")) {
					i++;
					commonVariables.pdfMin=Double.parseDouble(str[i]);
					i++;
					commonVariables.pdfMax=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-x0")) {
					i++;
					commonVariables.x0=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-eta")) {
					i++;
					commonVariables.eta=Double.parseDouble(str[i]);
				} else if(str[i].equalsIgnoreCase("-seed")) {
					i++;
					commonVariables.seed=Long.parseLong(str[i]);
				} else if(str[i].equalsIgnoreCase("-pars")) {
					i++;
					commonVariables.a=Double.parseDouble(str[i]);
					i++;
					commonVariables.b=Double.parseDouble(str[i]);
					i++;
					commonVariables.c=Double.parseDouble(str[i]);
				} else {
					System.out.println("Error in input!");
					System.out.println("-> "+str[i]+" <-");
					System.exit(0);
				}
				if(!(i<str.length-1)) System.out.println("Commandline arguments were fully processed.");
			}
		} catch (Exception e) {
			System.out.println("Error in input!");
			System.out.println("-> "+str[i]+" <-");
			System.exit(0);
		}
	}
	public static void main(String [] args) {
		getCommandLine(args);
		launcher gijos=new launcher();
	}
}

class launcher {
	private SimpleDateFormat sdf=new SimpleDateFormat("HH:mm:ss.SSS yyyy-MM-dd");
	private long startTime=0;//time (ms) when program was started
	private long timerPause=30000;//minimal pause (ms) between reports
	private long lastReport=0;//time (ms) of the last report
	private double[] pdf=new double[commonVariables.multiOutputPoints];
	private double[] spec=new double[commonVariables.realizationPoints];
	private int[] perc=null;//percentage
	private gija[] equations=null;//equation objects
	private int completed=0;//completed threads
	private int reported=0;//number of reports
	public launcher() {
		startTime=(Calendar.getInstance()).getTimeInMillis();
		System.out.println("Thread launcher has started!");
		int realZingsnis=Math.max(commonVariables.realizations/commonVariables.core,1);
		perc=new int[commonVariables.core];
		equations=new gija[commonVariables.core];
		if(commonVariables.seed<0) commonVariables.seed=System.currentTimeMillis();
		for(int i=0;i<commonVariables.core;i++) {
			perc[i]=0;
			equations[i]=new gija(commonVariables.seed,commonVariables.a,commonVariables.b,commonVariables.c,commonVariables.eta,commonVariables.dt,commonVariables.pdfMin,commonVariables.pdfMax,commonVariables.min,commonVariables.max,commonVariables.useLimits,commonVariables.x0,commonVariables.useAbs,commonVariables.nonLinearFeedback,commonVariables.absoluteNoise,commonVariables.gatherPoints,commonVariables.realizationPoints,realZingsnis,commonVariables.multiOutputPoints,commonVariables.outputTime,commonVariables.var,i,this);
			System.out.println("Thread nr. "+i+" was launched. Sent "+realZingsnis+" jobs.");
		}
	}
	//thread should pass histogram array, spectra array and its number
	synchronized public void put(double[] p, double[] s, int nr) {
		completed++;
		for(int i=0;i<pdf.length;i++) pdf[i]+=p[i];
		for(int i=0;i<spec.length;i++) spec[i]+=s[i];
		perc[nr]=100;
		if(completed>=commonVariables.core) {
			System.out.println("All threads has reported their results. Preparing results for final output.");
			int first=-1;
			int last=-1;
			for(int i=0;i<pdf.length;i++) {
				if((pdf[i]>0)&&(first==-1)&&(i>0)) first=i;
				if(pdf[i]>0) last=i;
				pdf[i]/=((double)commonVariables.core);
			}
			if(first>0) first--;
			if(last<pdf.length-1) last++;
			for(int i=0;i<spec.length;i++) spec[i]/=((double)commonVariables.core);
			double xstep=(commonVariables.pdfMax-commonVariables.pdfMin)/((double)commonVariables.multiOutputPoints);
			outputarr(commonFunctions.logPdfModification(pdf,commonVariables.pdfMin+first*xstep,commonVariables.pdfMin+last*xstep,commonVariables.pdfMin,xstep,commonVariables.outputPoints),commonVariables.outputTemplate+commonVariables.preffix+"dist",6,false);
			outputarr(commonFunctions.specModification(spec,commonVariables.gatherPoints*60*commonVariables.dt/commonVariables.dt60,commonVariables.outputPoints,false),commonVariables.outputTemplate+commonVariables.preffix+"spec",6,false);
			System.out.println("Done. Exiting.");
			System.exit(0);
		}
	}
	//function for a thread to report its progess (number, percentage)
	synchronized public void report(int nr, int vperc) {
		perc[nr]=vperc;
		Calendar cal=Calendar.getInstance();
		long nowMs=cal.getTimeInMillis();
		if(nowMs-lastReport>timerPause) {
			lastReport=nowMs;
			System.out.println(sdf.format(cal.getTime()));
			nowMs=cal.getTimeInMillis()-startTime;
			System.out.println("  Output template: "+commonVariables.outputTemplate);
			System.out.println("  Seed: "+commonVariables.seed);
			System.out.println("  Runtime: "+nowMs+" ms");
			int tmin=101;
			for(int i=0;i<perc.length;i++) {
				if(tmin>perc[i]) tmin=perc[i];
			}
			System.out.println("  Completion: "+tmin+"%");
			for(int i=0;i<perc.length;i++) {
				System.out.println("    ["+i+" - "+perc[i]+"%]");
			}
			long eta=(long)Math.floor((100.0/((double)tmin)-1.0)*nowMs);
			cal.setTimeInMillis(startTime+nowMs+eta);
			System.out.println("  ETA: "+sdf.format(cal.getTime())+" (+"+eta+" ms)");
			System.out.flush();
		}
	}
	//output 2D array
	public static void outputarr(double[][] arr, String name, int afterComma, boolean aboveZero) {
		BufferedWriter out=null;
		String fileName=name;
		try {
			File file=new File(fileName);
			file.createNewFile();
			out=new BufferedWriter(new FileWriter(fileName));
			for(int i=0;i<arr.length;i++) {
				if(((aboveZero)&&(arr[i][1]>0))||(!aboveZero)) {
					for(int j=0;j<arr[i].length-1;j++) out.write(round(arr[i][j],afterComma)+" ");
					out.write(round(arr[i][arr[i].length-1],afterComma)+"");
					out.newLine();
				}
			} 
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	//rounding
	public static double round(double d, int i) {
		int j = i;
		BigDecimal bigdecimal = new BigDecimal(d);
		bigdecimal = bigdecimal.setScale(j,bigdecimal.ROUND_HALF_UP);
		d = bigdecimal.doubleValue();
		return d;
	}
}

class gija implements Runnable {
	private static double matlog10=Math.log(10);
	private static double matlog2=Math.log(2);
	private double dt=1;
	private double x0=1;
	private double pdfMin=0.01;
	private double pdfMax=100;
	private double[] x=null;
	private double[] pdf=null;
	private double[] spec=null;
	private int realizationPoints=32768;
	private int realizations=1;
	private generalCarcass equation=null;
	private launcher parrent=null;
	private long reportNext=32768/100;
	private int proc=1;
	private int nr=1;
	private boolean outputTime=true;
	private int gatherPoints=1;
	public gija() {}
	public gija(long seed, double va, double vb, double vc, double veta, double vdt, double vpmin, double vpmax, double vfmin, double vfmax, boolean vf, double vx, boolean vabs, boolean lnfb, boolean vabsn, int vgp, int vrp, int vr, int vo, boolean outT, int var, int vnr, launcher vpar) {
		if(commonVariables.linear) {
			equation=new linearGarch(va,vb,vc,vfmin,vfmax,vf,vnr,vabs,seed);
		} else {
			equation=new nonLinearGarch(va,vb,vc,vfmin,vfmax,vf,vnr,veta,vabs,lnfb,vabsn,seed);
		}
		equation.var=var;
		pdfMin=vpmin;
		pdfMax=vpmax;
		dt=vdt;
		x0=Math.min(Math.max(vx,vfmin),vfmax);
		gatherPoints=vgp;
		realizationPoints=vrp;
		realizations=vr;
		reportNext=(int)Math.floor(proc*realizations*(((double)realizationPoints)/100.0));
		x=new double[realizationPoints];
		pdf=new double[vo];
		for(int i=0;i<pdf.length;i++) pdf[i]=0;
		spec=new double[realizationPoints];
		for(int i=0;i<spec.length;i++) spec[i]=0;
		parrent=vpar;
		nr=vnr;
		outputTime=outT;
		new Thread(this,"Thread nr. "+nr).start();
	}
	public void run() {
		double pdfStep=(pdfMax-pdfMin)/((double)pdf.length);
		if(Double.isNaN(x0)) {
			equation.setZ(0.5);
		} else equation.setZ(x0);
		for(int j=0;j<realizations;j++) {
			double[] pdft=new double[pdf.length];
			for(int i=0;i<pdf.length;i++) pdft[i]=0;
			long alreadyMade=j*realizationPoints;
			for(int i=0;i<realizationPoints;i++) {
				if(gatherPoints>1) {
					double grez=0;
					for(int g=0;g<gatherPoints;g++) {
						double eqrez=equation.step(dt);
						grez+=eqrez;
						eqrez=Math.log(eqrez)/matlog10;
						if((eqrez>=pdfMin)&&(eqrez<=pdfMax)) {
							int tmp=(int)Math.floor((eqrez-pdfMin)/pdfStep);
							if((tmp>=0)&&(tmp<pdft.length)) pdft[tmp]++;
						}						
					}
					x[i]=grez/((double)gatherPoints);
				} else {
					x[i]=equation.step(dt);
					double eqrez=Math.log(x[i])/matlog10;
					if((eqrez>=pdfMin)&&(eqrez<=pdfMax)) {
						int tmp=(int)Math.floor((eqrez-pdfMin)/pdfStep);
						if((tmp>=0)&&(tmp<pdft.length)) pdft[tmp]++;
					}
				}
				if((i+alreadyMade)>=reportNext) {
					parrent.report(nr,proc);
					proc++;
					reportNext=(int)Math.floor(proc*realizations*(((double)realizationPoints)/100.0));
				}
			}
			for(int i=0;i<pdf.length;i++) pdf[i]+=(pdft[i]/((double)realizationPoints));
			pdft=null;
			double[] spect=makeSpectra(preformFFT(realToComplex(x,true)));
			for(int i=0;i<=spec.length/2;i++) spec[i]+=spect[i];
		}
		if(realizations>1) {
			for(int i=0;i<pdf.length;i++) pdf[i]/=((double)realizations);
			for(int i=0;i<=spec.length/2;i++) spec[i]/=((double)realizations);
		}
		if(outputTime) outputarr(x,commonVariables.outputTemplate+commonVariables.preffix+nr+".time",3);
		parrent.put(pdf,spec,nr);
	}
	public static void outputarr(double[] arr, String name, int poKablelio) {
		BufferedWriter out=null;
		String fileName=name;
		try {
			File file=new File(fileName);
			file.createNewFile();
			out=new BufferedWriter(new FileWriter(fileName));
			for(int i=0;i<arr.length;i++) {
				out.write(commonFunctions.round(arr[i],poKablelio)+"");
				out.newLine();
			} 
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static double[][] realToComplex(double[] arrx, boolean subtractMean) {
		double[][] arr=new double[arrx.length][2];
		double mean=0;
		double tmean=0;
		for(int i=0;i<arr.length;i++) {
			arr[i][0]=arrx[i];
			arr[i][1]=0;
			if(subtractMean) {
				tmean+=arr[i][0];
				if(tmean>Double.MAX_VALUE/100.0) {
					mean+=(tmean/((double)arr.length));
					tmean=0;
				}
			}
		}
		if(subtractMean) {
			if(tmean!=0) mean+=(tmean/((double)arr.length));
			for(int i=0;i<arr.length;i++) arr[i][0]-=mean;
		}
		return arr;
	}
	public static double[][] preformFFT(double[][] arr) {
		int n=arr.length;
		int nm1=n-1;
		int nd2=(int)(n/2);
		int m=cint(Math.log(n)/matlog2);
		int j=nd2;
		for(int i=1;i<nm1;i++) {//bit reversal
			if(i<=j) {
				double tr=arr[j][0];
				double ti=arr[j][1];
				arr[j][0]=arr[i][0];
				arr[j][1]=arr[i][1];
				arr[i][0]=tr;
				arr[i][1]=ti;
			}
			int k=nd2;
			while(k<=j) {
				j-=k;
				k=(int)(k/2);
			}
			j+=k;
		}		
		for(int l=1;l<m+1;l++) {
			int le=cint(Math.pow(2,l));
			int le2=(int)(le/2);
			double ur=1;
			double ui=0;
			double tr=0;
			double ti=0;
			double sr=Math.cos(Math.PI/((double)le2));
			double si=-Math.sin(Math.PI/((double)le2));
			for(j=1;j<le2+1;j++) {
				int jm1=j-1;
				for(int i=jm1;i<=nm1;i+=le) {
					int ip=i+le2;
					tr=arr[ip][0]*ur-arr[ip][1]*ui;
					ti=arr[ip][1]*ur+arr[ip][0]*ui;
					arr[ip][0]=arr[i][0]-tr;
					arr[ip][1]=arr[i][1]-ti;
					arr[i][0]+=tr;
					arr[i][1]+=ti;
				}
				tr=ur;
				ur=tr*sr-ui*si;
				ui=tr*si+ui*sr;
			}
		}
		return arr;
	}
	public static double[] makeSpectra(double[][] arr) {
		double[] rez=new double[(int)Math.floor(arr.length/2)+1];
		for(int i=0;i<rez.length;i++) {
			rez[i]=(arr[i][0]*arr[i][0]+arr[i][1]*arr[i][1]);
		}
		return rez;
	}
	private static int cint(double expr) {
		double dif=Math.abs((expr-(int)expr));
		if(dif==0.5) {
			if(((int)expr)%2==0) {
				return (int)expr;
			} else {
				return (int)expr+1;
			}
		} else if(dif>0.5) {
			return (int)expr+1;
		} else {
			return (int)expr;
		}
	}
}

class generalCarcass {
	protected Random gen=new Random();
	protected int numeris=-1;
	public double a=1;
	public double b=1;
	public double c=1;
	public int var=0;
	protected double min=0.01;
	protected double max=100;
	protected boolean useLimits=true;
	protected boolean useAbs=false; 
	protected double kappa=1;
	protected double lastOmega=0.5;
	protected double lastSigma2=0.5;
	public generalCarcass() {}
	public double step(double dt) {return 0;}
	public void setZ(double val) {lastSigma2=val;lastOmega=gen.nextGaussian();}
}

class linearGarch extends generalCarcass {
	public linearGarch(double va, double vb, double vc, double vfmin, double vfmax, boolean vf, int vnr, boolean vabs, long seed) {
		a=va;
		b=vb;
		c=vc;
		min=vfmin;
		max=vfmax;
		useLimits=vf;
		useAbs=vabs;
		numeris=vnr;
		if(vnr==0) commonVariables.preffix="lin.";
		if(seed<0) gen=new Random(System.currentTimeMillis()+numeris*13);
		else gen=new Random(seed+numeris*13);
		makeOmega();
	}
	public double step(double dt) {
		lastSigma2=a+lastSigma2*(c+b*lastOmega*lastOmega);
		if(useLimits) lastSigma2=Math.max(Math.min(lastSigma2,max),min);
		if(useAbs) lastSigma2=Math.abs(lastSigma2);
		makeOmega();
		if(var==1) return Math.sqrt(lastSigma2)*lastOmega;
		return lastSigma2;
	}
	private void makeOmega() {
		lastOmega=gen.nextGaussian();
	}
}

class nonLinearGarch extends generalCarcass {
	private int eta=1;
	private boolean nonLinearFeedback=false;
	private boolean absoluteNoise=false;
	public nonLinearGarch(double va, double vb, double vc, double vfmin, double vfmax, boolean vf, int vnr, double veta, boolean vabs, boolean nlfb, boolean vabsn, long seed) {
		a=va;
		b=vb;
		c=vc;
		min=vfmin;
		max=vfmax;
		useLimits=vf;
		numeris=vnr;
		if(vnr==0) commonVariables.preffix="nlin.";
		useAbs=vabs;
		nonLinearFeedback=nlfb;
		absoluteNoise=vabsn;
		eta=(int)Math.floor(2.0*veta);//actual eta values are different from the internally used eta values
		if(seed<0) gen=new Random(System.currentTimeMillis()+numeris*13);
		else gen=new Random(seed+numeris*13);
		makeOmega();
	}
	public double step(double dt) {
		double zeta=Math.sqrt(lastSigma2)*lastOmega;//eta=0.5
		if(eta!=1) {//actually eta!=0.5
			//faster rasing to power in increments by 0.5
			double rez=1;
			for(int i=0;i<eta;i++) rez*=zeta;
			zeta=rez;
		}
		if(nonLinearFeedback) {
			double feedback=Math.sqrt(lastSigma2);//eta=0.5
			if(eta!=1) {//eta!=0.5
				double rez=1;
				for(int i=0;i<eta;i++) rez*=feedback;
				feedback=rez;
			}
			lastSigma2=a+b*zeta+lastSigma2-c*feedback;
		} else {
			lastSigma2=a+b*zeta+lastSigma2*c;
		}
		if(useLimits) lastSigma2=Math.max(Math.min(lastSigma2,max),min);
		if(useAbs) lastSigma2=Math.abs(lastSigma2);
		makeOmega();
		if(var==1) return Math.sqrt(lastSigma2)*lastOmega;
		return lastSigma2;
	}
	private void makeOmega() {
		lastOmega=gen.nextGaussian();
		if(absoluteNoise) lastOmega=Math.abs(lastOmega);
	}
}

class commonVariables {
	public static long seed=-1;
	public static int var=0;//which variable is outputed? (0: sigma_t^2; 1: z_t)
	public static int core=1;//number of processor cores (paralel execution if multicore processor)
	public static int realizations=1;//total number of realizations
	public static int realizationPoints=32768;//number of points in realization
	public static int gatherPoints=10;//number of points to fold into
	public static int outputPoints=640;//number of points to output in PDF or PSD
	public static int multiOutputPoints=10000;//internal number of points for estimating PDF
	public static String preffix="lin.";//output template (file prefix)
	public static String outputTemplate="rez.";//output template (file prefix)
	public static double pdfMin=0.001;//lower bound for PDF
	public static double pdfMax=0.999;//upper bound for PDF
	public static boolean outputTime=false;//output time series?
	public static double dt=1;//time interval between two points in realization
	public static double dt60=60;//how much time coincides with minute?
	public static double min=0.001;//limit from below
	public static double max=0.999;//limit from above
	public static boolean useLimits=true;//impose limits on the solution?
	public static boolean useAbs=false;//take absolute value of the solution?
	public static double x0=0;//initial condition
	public static double a=1;//model parameters
	public static double b=1;//model parameters
	public static double c=1;//model parameters
	public static double eta=1;//model parameters (mu = 2 eta)
	public static boolean linear=true;
	public static boolean absoluteNoise=false;
	public static boolean nonLinearFeedback=false;
}

class commonFunctions {
	private static double matlog10=Math.log(10);
	//Faster lg calculation
	public static double LogBase10(double a) {
		return Math.log(a)/matlog10;
	}
	//advanced rounding function
	public static double round(double d, int i) {
		int j=i;
		BigDecimal bigdecimal=new BigDecimal(d);
		bigdecimal=bigdecimal.setScale(j,bigdecimal.ROUND_HALF_UP);
		d=bigdecimal.doubleValue();
		return d;
	}
	//linear spectra to logarithmic spectra
	public static double[][] specModification(double[] spec, double timeTick, int outPoints, boolean smoothen) {
		double normavimas=commonFunctions.LogBase10(2*timeTick/spec.length);
		double scale=commonFunctions.LogBase10(spec.length)+commonFunctions.LogBase10(timeTick);
		double llim=0;
		double rlim=commonFunctions.LogBase10(spec.length/2.0);
		double lstep=(rlim-llim)/((double)outPoints);
		double clim=llim+lstep;
		int i=1;
		int intervale=0;
		double[][] rez=new double[outPoints+1][2];
		int used=0;
		double total=0;
		double oldX=0;
		while(clim<=rlim) {
			while(commonFunctions.LogBase10(i)<clim) {
				total+=spec[i];
				i++;
				intervale++;
			}
			if(total>0) {
				if(used==0) {
					oldX=Math.pow(10,clim-scale);
					rez[used][0]=commonFunctions.LogBase10(oldX/2.0);
					rez[used][1]=commonFunctions.LogBase10(total/((double)intervale));//oldX);
				} else {
					double newX=Math.pow(10,clim-scale);
					rez[used][0]=commonFunctions.LogBase10((newX+oldX)/2.0);
					rez[used][1]=commonFunctions.LogBase10(total/((double)intervale));//(newX-oldX));
					oldX=newX;
				}
				rez[used][1]+=normavimas;
				used++;
			}
			intervale=0;
			total=0;
			clim+=lstep;
		}
		double[][] rez2=new double[used][2];
		for(int ii=0;ii<used;ii++) {
			rez2[ii][0]=rez[ii][0];
			rez2[ii][1]=rez[ii][1];
		}
		if(smoothen) {
			/*movavg time window 3*/
			for(int ii=0;ii<used;ii++) {
				if((ii>0)&&(ii<used-1)) {
					rez2[ii][1]=(rez2[ii-1][1]+rez2[ii][1]+rez2[ii+1][1])/3.0;
				}
			}
		}
		return rez2;
	}
	//linear pdf to logarithmic pdf
	public static double[][] logPdfModification(double[] pdf, double llim, double rlim, double xlim, double xstep, int outPoints) {
		double[][] rez=new double[outPoints][2];
		for(int i=0;i<rez.length;i++) {
			rez[i][0]=0;
			rez[i][1]=0;
		}
		int moved=0;
		double curlim=xlim;
		while(curlim<llim) {
			curlim+=xstep;
			moved++;
		}
		double lstep=(rlim-llim)/((double)(outPoints-1));
		int used=0;
		while((llim<=rlim)&&(used<outPoints)) {
			double integral=0;
			llim+=lstep;
			while((curlim<llim)&&(moved<pdf.length)) {
				curlim+=xstep;
				integral+=pdf[moved];
				moved++;
			}
			if(integral>0) {
				rez[used][0]=llim-0.5*lstep;
				if(used>0) rez[used][1]=commonFunctions.LogBase10(integral/(Math.pow(10,rez[used][0])-Math.pow(10,rez[used-1][0])));
				else rez[used][1]=commonFunctions.LogBase10(integral/(Math.pow(10,rez[used][0])-Math.pow(10,rez[used][0]-lstep)));
				used++;
			}
		}
		if(used<outPoints) {
			double[][] rez2=new double[used][2];
			for(int ii=0;ii<used;ii++) {
				rez2[ii][0]=rez[ii][0];
				rez2[ii][1]=rez[ii][1];
			}
			rez=new double[used][2];
			for(int ii=0;ii<used;ii++) {
				rez[ii][0]=rez2[ii][0];
				rez[ii][1]=rez2[ii][1];
			}
			rez2=null;
		}
		return rez;
	}
}